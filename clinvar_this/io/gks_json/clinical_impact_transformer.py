"""Support for I/O of the AMP/ASCO/CAP 2017 GKS formatted data to define Clinical Impact submissions.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38

"""

from types import MappingProxyType

from ga4gh.core.models import MappableConcept, iriReference
from ga4gh.va_spec.aac_2017 import (
    AmpAscoCapClassificationCode,
    VariantClinicalSignificanceStatement,
)
from ga4gh.va_spec.base import (
    DiagnosticPredicate,
    EvidenceLine,
    PrognosticPredicate,
    TherapeuticResponsePredicate,
    TherapyGroup,
    VariantDiagnosticProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from logzero import logfile
from pydantic.dataclasses import dataclass

from clinvar_api.models import (
    Assembly,
    CitationDb,
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


@dataclass(slots=True)
class DrugContext:
    """Container for drug context information"""

    assertion_type: SomaticClinicalImpactAssertionType | None = None
    therapeutic: TherapyGroup | MappableConcept | None = None
    drug: str | None = None


class ClinicalImpactTransformer(GksJsonTransformer[VariantClinicalSignificanceStatement]):
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

    def _get_drug_for_therapeutic_assertion(
        self, therapeutic: MappableConcept | TherapyGroup
    ) -> str | None:
        """Get the name for drug(s)

        :param therapeutic: Therapeutic record. Assumes ``name`` is provided in
            ``MappableConcept`` objects.
        :return: The formatted name for a therapeutic record. Multiple therapies for
            combination and substitution will be separated by semicolons.
        """
        if isinstance(therapeutic, MappableConcept):
            return therapeutic.name

        return ";".join(sorted([t.name for t in therapeutic.therapies if t.name]))

    def _get_drug_context_from_evidence_lines(
        self, evidence_lines: list[EvidenceLine]
    ) -> DrugContext:
        """Retrieve drug context information found in evidence lines

        :param evidence_lines: Evidence lines containing potential drug context
            information
        :return: Drug context information, if found in evidence lines
        """
        for el in evidence_lines:
            if isinstance(el, iriReference):
                continue

            target_proposition = el.targetProposition

            if not isinstance(
                target_proposition,
                (
                    VariantDiagnosticProposition,
                    VariantPrognosticProposition,
                    VariantTherapeuticResponseProposition,
                ),
            ):
                continue

            assertion_type = self.gks_predicate_to_assertion[target_proposition.predicate]

            therapeutic = None
            drug = None

            if isinstance(target_proposition, VariantTherapeuticResponseProposition):
                therapeutic_root = target_proposition.objectTherapeutic.root

                if isinstance(therapeutic_root, (TherapyGroup, MappableConcept)):
                    therapeutic = therapeutic_root
                    drug = self._get_drug_for_therapeutic_assertion(therapeutic)

            return DrugContext(assertion_type=assertion_type, therapeutic=therapeutic, drug=drug)

        return DrugContext()

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
        evidence_lines = self._get_evidence_lines(statement.hasEvidenceLines)
        drug_context = self._get_drug_context_from_evidence_lines(evidence_lines)

        return SubmissionClinicalImpactSubmission(
            **self._build_shared_submission_kwargs(
                statement=statement,
                submitted_assembly=submitted_assembly,
                observed_in=observed_in,
                variant_hgvs=variant_hgvs,
            ),
            clinical_impact_classification=SomaticClinicalImpactClassification(
                **self._build_shared_classification_kwargs(
                    statement.description,
                    drug_context.therapeutic,
                    evidence_lines,
                    statement.contributions,
                ),
                clinical_impact_classification_description=_IMPACT_CLASS_MAPPING[
                    statement.classification.primaryCoding.code.root
                ],
                assertion_type_for_clinical_impact=drug_context.assertion_type,
                drug_for_therapeutic_assertion=drug_context.drug,
            ),
        )
