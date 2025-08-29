"""Support for I/O of the AMP/ASCO/CAP 2017 GKS Study Statements format to define submissions.

Currently only supports Therapeutic, Diagnostic, and Prognostic Assertions. This
assumes you are using MetaKB (https://github.com/cancervariants/metakb) to generate GKS
JSON files.
"""

from types import MappingProxyType
import typing

from ga4gh.core.models import MappableConcept
from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.va_spec.aac_2017 import (
    VariantTherapeuticResponseStudyStatement,
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
)
from ga4gh.va_spec.base import Contribution, EvidenceLine
from ga4gh.vrs.models import Syntax, Expression

from clinvar_api.models import (
    AffectedStatus,
    AlleleOrigin,
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

from clinvar_this.io.gks_json.base import GksJsonTransformer


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
            except StopIteration:
                return

        return get_hgvs(expressions, Syntax.HGVS_C) or get_hgvs(
            expressions, Syntax.HGVS_G
        )

    @staticmethod
    def _get_observed_in(
        allele_origin_qualifier: MappableConcept,
        method_type: str
    ) -> list[SubmissionObservedInSomatic] | None:
        """Get observed in value

        `collection_method` and `affected_status` are hard coded

        :param allele_origin_qualifier: Allele origin qualifier
        :param method_type: The specific type of method used for the study statement
        :return: The observed in value, allele origin mapping exists
        """
        allele_origin = allele_origin_qualifier.name
        if not allele_origin:
            return None

        try:
            collection_method = CollectionMethod(method_type)
        except ValueError:
            return None

        return [
            SubmissionObservedInSomatic(
                allele_origin=allele_origin,
                collection_method=collection_method,
                affected_status=AffectedStatus.UNKNOWN,
            )
        ]

    @staticmethod
    def _get_citations(
        mp_id: str,
        evidence_lines: typing.List[EvidenceLine],
    ) -> typing.List[SubmissionCitation]:
        """Get list of PubMed citations and supporting evidence for an assertion

        :param mp_id: CIViC molecular profile ID for assertion
        :param evidence_lines: A list of evidence-based arguments that may support or
            refute the validity of a proposition for a CIViC assertion
        :return: List of PubMed citations and supporting CIViC evidence URLs for an
            assertion
        """
        citations = [SubmissionCitation(url=f"https://identifiers.org/{mp_id}")]
        for evidence_line in evidence_lines:
            for evidence_item in getattr(evidence_line, "hasEvidenceItems", []):
                citations.append(
                    SubmissionCitation(
                        url=f"https://civicdb.org/links/evidence/{evidence_item.id.split('civic.eid:')[-1]}"
                    )
                )

                docs = evidence_item.reportedIn or []
                for doc in docs:
                    if getattr(doc, "pmid", None):
                        citations.append(
                            SubmissionCitation(
                                url=f"https://pubmed.ncbi.nlm.nih.gov/{doc.pmid}"
                            )
                        )

        return citations

    @staticmethod
    def _get_date_last_evaluated(contributions: list[Contribution]) -> str | None:
        """The date that the classification was last evaluated by the submitter

        :param contributions: List of contributions
        :return: The date (yyyy-mm-dd) the classification was last evaluated by the
            submitter
        """
        try:
            contribution: Contribution = contributions[0]
        except IndexError:
            return

        return contribution.date.strftime("%Y-%m-%d")

    def _get_clinical_impact_submission(
        self,
        record: VariantTherapeuticResponseStudyStatement
        | VariantDiagnosticStudyStatement
        | VariantPrognosticStudyStatement,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
    ) -> SubmissionClinicalImpactSubmission:
        """Get clinical impact submission for a therapeutic, diagnostic, or prognostic assertion

        Assertions with Substitutes therapies will be separated by semicolons
        and the Study Statement's description will be updated to include this note.

        :param record: The therapeutic, diagnostic, or prognostic assertion
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found. If not found,
            ``variant_coords`` must be provided. This takes priority over
            ``variant_coords``.
        :param variant_coords: The chromosome coordinates for a variant, if found. If
            not found, ``variant_hgvs`` must be provided
        :return: The clinical impact submission for a therapeutic, diagnostic, or
            prognostic assertion
        """
        proposition = record.proposition

        if hasattr(proposition, "objectTherapeutic"):
            therapeutic = proposition.objectTherapeutic.root
            drug_for_therapeutic_assertion = self.get_drugs(therapeutic)
        else:
            drug_for_therapeutic_assertion = None
        mp_id = proposition.subjectVariant.id

        return SubmissionClinicalImpactSubmission(
            record_status=RecordStatus.NOVEL,
            local_id=mp_id,
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
                citation=self._get_citations(mp_id, record.hasEvidenceLines),
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
    ) -> typing.List[SubmissionContainer]:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param study_statements: List of Therapeutic Response, Diagnostic,
            or Prognostic Assertions represented as GKS Variant Therapeutic Response,
            Diagnostic, or Prognostic Study Statements
        :return: A list of submission container data structures
        """
        clinical_impact_submissions: typing.List[SubmissionContainer] = []

        for study_statement in study_statements:
            variant = study_statement.proposition.subjectVariant
            variant_hgvs = self._get_variant_hgvs(variant)
            if not variant_hgvs:
                continue

            observed_in = self._get_observed_in(
                study_statement.proposition.alleleOriginQualifier,
                study_statement.specifiedBy.methodType
            )
            if not observed_in:
                continue

            clinical_impact_submission = self._get_clinical_impact_submission(
                study_statement,
                observed_in,
                variant_hgvs=variant_hgvs,
            )
            clinical_impact_submissions.append(clinical_impact_submission)

        return SubmissionContainer(
            assertion_criteria=SubmissionAssertionCriteria(
                db=CitationDb.PUBMED,
                id="27993330",  # AMP/ASCO/CAP
            ),
            clinical_impact_submission=clinical_impact_submissions,
        )
