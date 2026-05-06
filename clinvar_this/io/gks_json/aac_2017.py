"""Support for I/O of the AMP/ASCO/CAP 2017 GKS formatted data to define submissions.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38

"""

from types import MappingProxyType
import typing
from logzero import logger, logfile
from ga4gh.core.models import MappableConcept, iriReference
from ga4gh.cat_vrs.models import CategoricalVariant, DefiningAlleleConstraint
from ga4gh.va_spec.aac_2017 import (
    VariantClinicalSignificanceStatement,
    AmpAscoCapClassificationCode,
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
    AlleleOrigin,
    SomaticClinicalImpactClassificationDescription,
)

from clinvar_this.io.gks_json.base import BatchMetadata, GksJsonTransformer

logfile("aac_2017.log")


# Mapping from GKS classification code to ClinVar clinical impact classification
_IMPACT_CLASS_MAPPING = MappingProxyType(
    {
        AmpAscoCapClassificationCode.TIER_1: SomaticClinicalImpactClassificationDescription.STRONG,
        AmpAscoCapClassificationCode.TIER_2: SomaticClinicalImpactClassificationDescription.POTENTIAL,
        AmpAscoCapClassificationCode.TIER_3: SomaticClinicalImpactClassificationDescription.UNKNOWN,
        AmpAscoCapClassificationCode.TIER_4: SomaticClinicalImpactClassificationDescription.BENIGN_LIKELY_BENIGN,
    }
)


class Aac2017GksJsonTransformer(GksJsonTransformer):
    """Class for transforming AMP/ASCO/CAP 2017 GKS formatted data to define submissions"""

    @staticmethod
    def _get_variant_hgvs(
        variant: CategoricalVariant,
    ) -> str | None:
        """Retrieve a HGVS expression for a variant

        Checks the first constraint for an expression. Only support
        DefiningAlleleConstraints at the moment.

        If no constraints found, then checks the expressions extension.

        Order matters: the first matching expression is returned. cDNA RefSeq HGVS
        expressions are prioritized over genomic RefSeq HGVS expressions.

        :param variant: Categorical Variant
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
        if getattr(variant, "constraints", None) and variant.constraints:
            constraint = variant.constraints[0]
            if isinstance(constraint.root, DefiningAlleleConstraint):
                expressions = constraint.root.allele.expressions
        else:
            # Case where a VRS Allele is unable to be represented, can store
            # expressions as an extension named 'expressions'
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

        :param allele_origin_qualifier: Allele origin qualifier
        :param collection_method: The specific type of method used for the statement
        :param affected_status: Whether or not the individual(s) in each observation
            were affected by the condition for the interpretation
        :return: The observed in value
        """
        allele_origin: str | None = allele_origin_qualifier.name
        if not allele_origin:
            return None

        return [
            SubmissionObservedInSomatic(
                allele_origin=AlleleOrigin(allele_origin),
                collection_method=collection_method,
                affected_status=affected_status,
            )
        ]

    @staticmethod
    def _get_citations(
        evidence_lines: list[EvidenceLine],
    ) -> list[SubmissionCitation]:
        """Extract unique citations from evidence lines and related evidence items.

        Citations may be sourced from:
        - `citations` extensions attached to evidence lines
        - PubMed IDs (`pmid`) on reported documents
        - IRI reference URLs attached to reported documents

        If all reported documents contain PMIDs, citations are submitted using
        PubMed database identifiers. Otherwise, citations fall back to URL-based
        references.

        :param evidence_lines: Evidence lines that may contain supporting citation
            information.
        :return: A deduplicated list of submission citations.
        """

        citations: list[SubmissionCitation] = []
        reported_in_documents: list[Document | iriReference] = []

        def add_citation(citation: SubmissionCitation) -> None:
            """Append a citation if it has not already been added."""
            if citation not in citations:
                citations.append(citation)

        for evidence_line in evidence_lines:
            # Temporary support for citations stored in extensions
            for ext in evidence_line.extensions or []:
                if ext.name == "citations" and isinstance(ext.value, list):
                    for citation_url in ext.value:
                        add_citation(SubmissionCitation(url=citation_url))

            reported_in_documents.extend(evidence_line.reportedIn or [])

            for evidence_item in evidence_line.hasEvidenceItems or []:
                reported_in_documents.extend(evidence_item.reportedIn or [])

        documents_have_all_pmids = all(
            isinstance(document, Document) and document.pmid
            for document in reported_in_documents
        )

        if documents_have_all_pmids:
            for document in reported_in_documents:
                add_citation(
                    SubmissionCitation(
                        db=CitationDb.PUBMED,
                        id=document.pmid,
                    )
                )

            return citations

        for document in reported_in_documents:
            if isinstance(document, Document):
                if document.pmid:
                    add_citation(
                        SubmissionCitation(
                            url=f"https://pubmed.ncbi.nlm.nih.gov/{document.pmid}"
                        )
                    )
            elif isinstance(document, iriReference):
                add_citation(SubmissionCitation(url=document.root))

        return citations

    @staticmethod
    def _get_date_last_evaluated(contributions: list[Contribution]) -> str | None:
        """The date that the classification was last evaluated by the submitter

        Expects `contributions` to be ordered chronologically, with the most
        recent contribution as the final item in the list.

        :param contributions: Ordered list of contributions
        :return: The formatted date (`YYYY-MM-DD`) of the most recent contribution,
            or `None` if no contributions are present.
        """
        try:
            contribution: Contribution = contributions[-1]
        except IndexError:
            return None

        return contribution.date.strftime("%Y-%m-%d")

    def _get_clinical_impact_submission(
        self,
        record: VariantClinicalSignificanceStatement,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ) -> SubmissionClinicalImpactSubmission:
        """Get clinical impact novel submission for a therapeutic, diagnostic, or prognostic assertion

        Assertions with Substitutes therapies will be separated by semicolons
        and the Statement's description will be updated to include this note.

        Local ID will use the proposition's variant ID or name.

        Local Key will use the `record`'s ID.

        :param record: The therapeutic, diagnostic, or prognostic assertion
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: The clinical impact submission for a therapeutic, diagnostic, or
            prognostic assertion
        """
        proposition = record.proposition
        target_proposition = record.hasEvidenceLines[0].targetProposition
        if hasattr(target_proposition, "objectTherapeutic"):
            therapeutic = target_proposition.objectTherapeutic.root
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
                    target_proposition.predicate
                ],
                comment=self.get_comment(record),
                citation=self._get_citations(record.hasEvidenceLines),
                drug_for_therapeutic_assertion=drug_for_therapeutic_assertion,
                date_last_evaluated=self._get_date_last_evaluated(
                    record.contributions or []
                ),
            ),
        )

    def records_to_submission_container(
        self,
        statements: typing.List[VariantClinicalSignificanceStatement],
        batch_metadata: BatchMetadata,
    ) -> SubmissionContainer:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param statements: List of Therapeutic Response, Diagnostic, or Prognostic
            statements
        :param batch_metadata: Batch-wide settings
            The properties will be assigned to all variants/samples in the batch.
        :return: Submission container data structures
        """
        clinical_impact_submissions: typing.List[
            SubmissionClinicalImpactSubmission
        ] = []

        for statement in statements:
            variant = statement.proposition.subjectVariant
            variant_hgvs = self._get_variant_hgvs(variant)
            if not variant_hgvs:
                logger.warning(
                    "Skipping statement. No hgvs found for statement ID: %s",
                    statement.id,
                )
                continue

            if stmt_method_type := statement.specifiedBy.methodType:
                method_type = CollectionMethod(stmt_method_type)
            else:
                method_type = batch_metadata.collection_method

            observed_in = self._get_observed_in(
                statement.proposition.alleleOriginQualifier,
                method_type,
                batch_metadata.affected_status,
            )
            if not observed_in:
                logger.warning(
                    "Skipping statement. No observed in found for statement ID: %s",
                    statement.id,
                )
                continue

            clinical_impact_submission: SubmissionClinicalImpactSubmission = (
                self._get_clinical_impact_submission(
                    statement,
                    observed_in,
                    variant_hgvs=variant_hgvs,
                    submitted_assembly=batch_metadata.submitted_assembly,
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
