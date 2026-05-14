"""Support for I/O of the AMP/ASCO/CAP 2017 GKS formatted data to define Clinical Impact submissions.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38

"""

from types import MappingProxyType
from logzero import logfile
from ga4gh.va_spec.aac_2017 import (
    VariantClinicalSignificanceStatement,
    AmpAscoCapClassificationCode,
)
from ga4gh.va_spec.base import (
    DiagnosticPredicate,
    PrognosticPredicate,
    TherapeuticResponsePredicate,
)

from clinvar_api.models import (
    Assembly,
    CitationDb,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionClinicalImpactSubmission,
)
from clinvar_api.models.sub_payload import (
    SomaticClinicalImpactClassification,
    SubmissionObservedInSomatic,
)
from clinvar_api.msg.sub_payload import (
    SomaticClinicalImpactAssertionType,
    SomaticClinicalImpactClassificationDescription,
)

from clinvar_this.io.gks_json.base import GksJsonTransformer

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


class ClinicalImpactTransformer(
    GksJsonTransformer[VariantClinicalSignificanceStatement]
):
    """Class for transforming AMP/ASCO/CAP 2017 GKS formatted data to define Clinical Impact submissions"""

    submission_container_attribute = "clinical_impact_submission"
    assertion_criteria = SubmissionAssertionCriteria(
        db=CitationDb.PUBMED,
        id="27993330",  # AMP/ASCO/CAP Guidelines, 2017
    )
    gks_statement_cls = VariantClinicalSignificanceStatement

    # Mapping from GKS classification code to ClinVar clinical impact classification
    impact_class_mapping = MappingProxyType(
        {
            AmpAscoCapClassificationCode.TIER_1: SomaticClinicalImpactClassificationDescription.STRONG,
            AmpAscoCapClassificationCode.TIER_2: SomaticClinicalImpactClassificationDescription.POTENTIAL,
            AmpAscoCapClassificationCode.TIER_3: SomaticClinicalImpactClassificationDescription.UNKNOWN,
            AmpAscoCapClassificationCode.TIER_4: SomaticClinicalImpactClassificationDescription.BENIGN_LIKELY_BENIGN,
        }
    )

    # Mapping from GKS predicate type to ClinVar assertion type for clinical impact
    gks_predicate_to_assertion = {
        TherapeuticResponsePredicate.RESISTANCE: SomaticClinicalImpactAssertionType.THERAPEUTIC_RESISTANCE,
        TherapeuticResponsePredicate.SENSITIVITY: SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
        DiagnosticPredicate.EXCLUSIVE: SomaticClinicalImpactAssertionType.DIAGNOSTIC_EXCLUDES_DIAGNOSIS,
        DiagnosticPredicate.INCLUSIVE: SomaticClinicalImpactAssertionType.DIAGNOSTIC_SUPPORTS_DIAGNOSIS,
        PrognosticPredicate.BETTER_OUTCOME: SomaticClinicalImpactAssertionType.PROGNOSTIC_BETTER_OUTCOME,
        PrognosticPredicate.WORSE_OUTCOME: SomaticClinicalImpactAssertionType.PROGNOSTIC_POOR_OUTCOME,
    }

    def _get_submission(
        self,
        statement: VariantClinicalSignificanceStatement,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ) -> SubmissionClinicalImpactSubmission:
        """Transform a GKS clinical significance statement into a ClinVar novel clinical impact submission.

        These statements support AMP/ASCO/CAP therapeutic, diagnostic, and prognostic assertions.

        Assertions with Substitutes therapies will be separated by semicolons
        and the Statement's description will be updated to include this note.

        Local ID will use the proposition's variant ID or name.

        Local Key will use the `record`'s ID.

        :param statement: GKS statement (therapeutic, diagnostic, or prognostic) to transform
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: The clinical impact submission corresponding to a GKS Clinical
            Significance statement
        """
        proposition = statement.proposition
        target_proposition = statement.hasEvidenceLines[0].targetProposition
        if hasattr(target_proposition, "objectTherapeutic"):
            therapeutic = target_proposition.objectTherapeutic.root
            drug_for_therapeutic_assertion = self._get_drugs(therapeutic)
        else:
            drug_for_therapeutic_assertion = None

        return SubmissionClinicalImpactSubmission(
            record_status=RecordStatus.NOVEL,
            local_id=proposition.subjectVariant.id or proposition.subjectVariant.name,
            submitted_assembly=submitted_assembly,
            local_key=statement.id,
            observed_in=observed_in,
            condition_set=self._get_condition_set(proposition),
            variant_set=self._get_variant_set(proposition, variant_hgvs=variant_hgvs),
            clinical_impact_classification=SomaticClinicalImpactClassification(
                clinical_impact_classification_description=_IMPACT_CLASS_MAPPING[
                    statement.classification.primaryCoding.code.root
                ],
                assertion_type_for_clinical_impact=self.gks_predicate_to_assertion[
                    target_proposition.predicate
                ],
                comment=self._get_comment(statement),
                citation=self._get_citations(statement.hasEvidenceLines),
                drug_for_therapeutic_assertion=drug_for_therapeutic_assertion,
                date_last_evaluated=self._get_date_last_evaluated(
                    statement.contributions or []
                ),
            ),
        )
