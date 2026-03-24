"""Support for I/O of the AMP/ASCO/CAP 2017 GKS Study Statements format to define submissions.

Currently only supports Therapeutic, Diagnostic, and Prognostic Assertions. This
assumes you are using MetaKB (https://github.com/cancervariants/metakb) to generate GKS
JSON files.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38

"""

from types import MappingProxyType
import typing
from logzero import logger
from ga4gh.core.models import MappableConcept, iriReference
from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.va_spec.aac_2017 import (
    VariantTherapeuticResponseStudyStatement,
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
)
from ga4gh.va_spec.base import Contribution, Document, EvidenceLine
from ga4gh.vrs.models import Syntax, Expression

from clinvar_api.models import (
    AffectedStatus,
    Assembly,
    CitationDb,
    CollectionMethod,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionCitation,
    SubmissionContainer,
    SubmissionClinicalImpactSubmission,
)
from clinvar_api.models.sub_payload import (
    SomaticClinicalImpactClassification,
    SubmissionObservedInSomatic,
)
from clinvar_api.msg.sub_payload import (
    SomaticClinicalImpactClassificationDescription,
)

from clinvar_this.io.gks_json.base import BatchMetadata, GksJsonTransformer


# Mapping from GKS classification to clinical impact classification
_IMPACT_CLASS_MAPPING = MappingProxyType(
    {
        "Tier I": SomaticClinicalImpactClassificationDescription.STRONG,
        "Tier II": SomaticClinicalImpactClassificationDescription.POTENTIAL,
        "Tier III": SomaticClinicalImpactClassificationDescription.UNKNOWN,
        "Tier IV": SomaticClinicalImpactClassificationDescription.BENIGN_LIKELY_BENIGN,
    }
)


class Aac2017GksJsonTransformer(GksJsonTransformer):
    """Class for transforming AMP/ASCO/CAP 2017 GKS Study Statements into submission format"""

    @staticmethod
    def _get_variant_hgvs(
        variant: CategoricalVariant,
    ) -> str | None:
        """Get HGVS expression for a variant

        :param variant: Variant record
        :return: cDNA RefSeq HGVS expression or genomic RefSeq HGVS expression for a
            variant, if provided.
        """

        def get_hgvs(
            expressions: typing.List[Expression],
            syntax: typing.Union[Syntax.HGVS_C, Syntax.HGVS_G],
        ) -> str | None:
            """Get a HGVS expression from a list of expressions for a given syntax

            :param expressions: List of representations specified by nomenclature or
                syntax for a variant
            :param syntax: The syntax to find an expression for
            :return: HGVS expression for a given syntax, if found
            """
            expr_prefix = "NM" if syntax == Syntax.HGVS_C else "NC"
            hgvs = [
                expr.value
                for expr in expressions
                if expr.syntax == syntax and expr.value.startswith(expr_prefix)
            ]
            if hgvs:
                return sorted(hgvs)[0]

            return None

        expressions = None
        if getattr(variant, "constraints", None):
            expressions = variant.constraints[0].root.allele.expressions
        else:
            # Not able to normalize
            try:
                expressions_ext = next(
                    ext for ext in variant.extensions if ext.name == "expressions"
                ).value
                expressions = [Expression(**ext) for ext in expressions_ext]
            except (StopIteration, TypeError):
                return None

        return get_hgvs(expressions, Syntax.HGVS_C) or get_hgvs(
            expressions, Syntax.HGVS_G
        )

    @staticmethod
    def _get_observed_in(
        allele_origin_qualifier: MappableConcept,
        collection_method: CollectionMethod,
        affected_status: AffectedStatus,
    ) -> list[SubmissionObservedInSomatic] | None:
        """Get observed in value

        `collection_method` and `affected_status` are hard coded

        :param allele_origin_qualifier: Allele origin qualifier
        :param collection_method: The specific type of method used for the study statement
        :return: The observed in value, allele origin mapping exists
        """
        allele_origin = allele_origin_qualifier.name
        if not allele_origin:
            return None

        return [
            SubmissionObservedInSomatic(
                allele_origin=allele_origin,
                collection_method=collection_method,
                affected_status=affected_status,
            )
        ]

    @staticmethod
    def _get_citations(
        evidence_lines: typing.List[EvidenceLine],
    ) -> typing.List[SubmissionCitation]:
        """Get list of PubMed citations and supporting evidence for an assertion

        :param evidence_lines: A list of evidence-based arguments that may support or
            refute the validity of a proposition for a study statement
        :return: List of PubMed citations and supporting evidence URLs for a study
            statement
        """
        citations = []

        reported_in_documents: list[Document] = []
        for evidence_line in evidence_lines:
            if reported_in := evidence_line.reportedIn or []:
                reported_in_documents.extend(reported_in)

            if evidence_items := evidence_line.hasEvidenceItems:
                for evidence_item in evidence_items:
                    if reported_in := evidence_item.reportedIn or []:
                        reported_in_documents.extend(reported_in)

        try:
            if all((document.pmid for document in reported_in_documents)):
                # Use DB ID
                for document in reported_in_documents:
                    submission_citation = SubmissionCitation(
                        db=CitationDb.PUBMED, id=document.pmid
                    )
                    if submission_citation not in citations:
                        citations.append(submission_citation)
                return citations
        except AttributeError:
            for document in reported_in_documents:
                if isinstance(document, Document):
                    if document_pmid := document.pmid:
                        citations.append(
                            SubmissionCitation(
                                url=f"https://pubmed.ncbi.nlm.nih.gov/{document_pmid}"
                            )
                        )
                elif isinstance(document, iriReference):
                    citations.append(SubmissionCitation(url=document.root))

        return citations

    @staticmethod
    def _get_date_last_evaluated(contributions: list[Contribution]) -> str | None:
        """The date that the classification was last evaluated by the submitter

        :param contributions: List of contributions
        :return: The last record in the list of contributions
        """
        try:
            contribution: Contribution = contributions[-1]
        except IndexError:
            return None

        return contribution.date.strftime("%Y-%m-%d")

    def _get_clinical_impact_submission(
        self,
        record: VariantTherapeuticResponseStudyStatement
        | VariantDiagnosticStudyStatement
        | VariantPrognosticStudyStatement,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None
    ) -> SubmissionClinicalImpactSubmission:
        """Get clinical impact submission for a therapeutic, diagnostic, or prognostic assertion

        Assertions with Substitutes therapies will be separated by semicolons
        and the Study Statement's description will be updated to include this note.

        :param record: The therapeutic, diagnostic, or prognostic assertion
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found. If not found,
            ``variant_coords`` must be provided. This takes priority over
            ``variant_coords``.
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: The clinical impact submission for a therapeutic, diagnostic, or
            prognostic assertion
        """
        proposition = record.proposition

        if hasattr(proposition, "objectTherapeutic"):
            therapeutic = proposition.objectTherapeutic.root
            drug_for_therapeutic_assertion = self.get_drugs(therapeutic)
        else:
            drug_for_therapeutic_assertion = None

        return SubmissionClinicalImpactSubmission(
            record_status=RecordStatus.NOVEL,
            local_id=proposition.subjectVariant.id or proposition.subjectVariant.name,
            submitted_assembly=submitted_assembly,
            local_key=record.id,
            observed_in=observed_in,
            condition_set=self.get_condition_set(proposition),
            variant_set=self.get_variant_set(proposition, variant_hgvs=variant_hgvs),
            clinical_impact_classification=SomaticClinicalImpactClassification(
                clinical_impact_classification_description=_IMPACT_CLASS_MAPPING[
                    record.classification.primaryCoding.code.root
                ],
                assertion_type_for_clinical_impact=self.gks_predicate_to_assertion[
                    proposition.predicate
                ],
                comment=self.get_comment(record),
                citation=self._get_citations(record.hasEvidenceLines),
                drug_for_therapeutic_assertion=drug_for_therapeutic_assertion,
                date_last_evaluated=self._get_date_last_evaluated(record.contributions),
            ),
        )

    def records_to_submission_container(
        self,
        study_statements: typing.List[
            VariantTherapeuticResponseStudyStatement
            | VariantDiagnosticStudyStatement
            | VariantPrognosticStudyStatement
        ],
        batch_metadata: BatchMetadata,
    ) -> SubmissionContainer:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param study_statements: List of Therapeutic Response, Diagnostic,
            or Prognostic Assertions represented as GKS Variant Therapeutic Response,
            Diagnostic, or Prognostic Study Statements
        :param batch_metadata: Batch-wide settings
            The properties will be assigned to all variants/samples in the batch.
        :return: Submission container data structures
        """
        clinical_impact_submissions: typing.List[
            SubmissionClinicalImpactSubmission
        ] = []

        for study_statement in study_statements:
            variant = study_statement.proposition.subjectVariant
            variant_hgvs = self._get_variant_hgvs(variant)
            if not variant_hgvs:
                logger.warning("No hgvs found for statement ID: %s", study_statement.id)
                continue

            if stmt_method_type := study_statement.specifiedBy.methodType:
                method_type = CollectionMethod(stmt_method_type)
            else:
                method_type = batch_metadata.collection_method

            observed_in = self._get_observed_in(
                study_statement.proposition.alleleOriginQualifier,
                method_type,
                batch_metadata.affected_status,
            )
            if not observed_in:
                logger.warning(
                    "No observed in found for statement ID: %s", study_statement.id
                )
                continue

            clinical_impact_submission: SubmissionClinicalImpactSubmission = (
                self._get_clinical_impact_submission(
                    study_statement,
                    observed_in,
                    variant_hgvs=variant_hgvs,
                    submitted_assembly=batch_metadata.submitted_assembly
                )
            )
            clinical_impact_submissions.append(clinical_impact_submission)

        return SubmissionContainer(
            assertion_criteria=SubmissionAssertionCriteria(
                db=CitationDb.PUBMED,
                id="27993330",  # AMP/ASCO/CAP
            ),
            clinical_impact_submission=clinical_impact_submissions,
        )
