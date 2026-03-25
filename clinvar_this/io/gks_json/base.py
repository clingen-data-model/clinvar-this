"""Support for I/O of the GKS JSON format to define submissions.

Currently only supports Therapeutic, Diagnostic, and Prognostic Assertions.
"""

from abc import ABC, abstractmethod
import json
import typing

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import (
    VariantDiagnosticStudyStatement,
    VariantPrognosticStudyStatement,
    VariantTherapeuticResponseStudyStatement,
)
from ga4gh.va_spec.base import (
    DiagnosticPredicate,
    VariantDiagnosticProposition,
    VariantPrognosticProposition,
    PrognosticPredicate,
    MembershipOperator,
    TherapeuticResponsePredicate,
    TherapyGroup,
    VariantTherapeuticResponseProposition,
)
from ga4gh.vrs.models import MolecularVariation
from pydantic import BaseModel, ConfigDict, ValidationError


from clinvar_api.models import (
    SubmissionCondition,
    SubmissionContainer,
    SubmissionVariant,
    SubmissionVariantSet,
    CollectionMethod,
    AffectedStatus,
    Assembly
)
from clinvar_api.models.sub_payload import (
    SubmissionConditionSetSomatic,
    SubmissionVariantGene,
)
from clinvar_api.msg.sub_payload import (
    ConditionDb,
    SomaticClinicalImpactAssertionType,
    SomaticClinicalImpactClassificationDescription,
)
from clinvar_this import exceptions
from clinvar_this.io.base import TransformIO


class BatchMetadata(BaseModel):
    """Batch-wide settings for GKS JSON import.

    The properties will be assigned to all variants/samples in the batch.
    """

    model_config = ConfigDict(frozen=True)

    affected_status: AffectedStatus = AffectedStatus.UNKNOWN
    collection_method: CollectionMethod = CollectionMethod.NOT_PROVIDED
    submitted_assembly: Assembly = Assembly.GRCH38


def batch_metadata_from_mapping(
    keys_values: typing.Iterable[str],
) -> BatchMetadata:
    """Convert configuration from ``KEY=VALUE`` strings to ``BatchMetadata``

    If values are not provided, then will use defaults
    """
    field_types = {
        name: value for (name, value) in typing.get_type_hints(BatchMetadata).items()
    }
    kwargs = {}
    for key_value in keys_values:
        if "=" not in key_value:
            raise exceptions.ArgumentsError(f"Invalid key/value pair in {key_value}")
        key, value = key_value.split("=")
        if key in field_types:
            try:
                kwargs[key] = field_types[key](value)
            except ValueError as e:
                raise exceptions.ArgumentsError(e)

    return BatchMetadata(**kwargs)


class GksJsonTransformer(TransformIO, ABC):
    """Class for transforming GKS JSON input data from various formats into submission format"""

    # Mapping from GKS predicate type to assertion type for clinical impact
    gks_predicate_to_assertion = {
        TherapeuticResponsePredicate.RESISTANCE: SomaticClinicalImpactAssertionType.THERAPEUTIC_RESISTANCE,
        TherapeuticResponsePredicate.SENSITIVITY: SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
        DiagnosticPredicate.EXCLUSIVE: SomaticClinicalImpactAssertionType.DIAGNOSTIC_EXCLUDES_DIAGNOSIS,
        DiagnosticPredicate.INCLUSIVE: SomaticClinicalImpactAssertionType.DIAGNOSTIC_SUPPORTS_DIAGNOSIS,
        PrognosticPredicate.BETTER_OUTCOME: SomaticClinicalImpactAssertionType.PROGNOSTIC_BETTER_OUTCOME,
        PrognosticPredicate.WORSE_OUTCOME: SomaticClinicalImpactAssertionType.PROGNOSTIC_POOR_OUTCOME,
    }

    @abstractmethod
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

        :param study_statements: List of GKS study statements
        :return: Submission container data structures
        """

    @staticmethod
    def _read_file(
        inputf: typing.TextIO,
    ) -> typing.List[
        VariantTherapeuticResponseStudyStatement
        | VariantDiagnosticStudyStatement
        | VariantPrognosticStudyStatement
    ]:
        """Get list of Variant Therapeutic Response, Diagnostic, or Prognostic Study Statements from a file

        For now, prognostic study statements are NOT supported

        :param inputf: Text file-like object containing input GKS Statement data
        :raises exceptions.InvalidFormat: If there was an error decoding JSON
        :return: A list of Variant Therapeutic Response, Diagnostic, or Prognostic Study Statements
        """
        study_statements: typing.List[
            VariantTherapeuticResponseStudyStatement
            | VariantDiagnosticStudyStatement
            | VariantPrognosticStudyStatement
        ] = []

        try:
            data = json.load(inputf)
        except json.JSONDecodeError as e:
            err_msg = "Error decoding GKS JSON"
            raise exceptions.InvalidFormat(err_msg) from e

        if "gks_records" not in data:
            err_msg = "Invalid GKS JSON: missing required key `gks_records` (must be a list of statements)"
            raise KeyError(err_msg)

        statements = data["gks_records"]

        supported_study_stmts: typing.List[callable] = [
            VariantTherapeuticResponseStudyStatement,
            VariantDiagnosticStudyStatement,
            VariantPrognosticStudyStatement,
        ]

        for statement in statements:
            for study_stmt in supported_study_stmts:
                try:
                    study_statements.append(study_stmt(**statement))
                except ValidationError:
                    pass

        return study_statements

    @staticmethod
    def get_variant_aliases(
        variant: MolecularVariation | CategoricalVariant,
    ) -> typing.List[str] | None:
        """Get the aliases for a variant

        :param variant: Variant record
        :return: Aliases for a variant, if found
        """
        return variant.aliases

    @staticmethod
    def get_drugs(therapeutic: TherapyGroup | MappableConcept) -> str:
        """Get the name for drug(s)

        :param therapeutic: Therapeutic record. Assumes ``name`` is provided in
            ``MappableConcept`` objects.
        :return: The formatted name for a therapeutic record. Multiple therapies for
            combination and substitution will be separated by semicolons.
        """
        if isinstance(therapeutic, MappableConcept):
            return therapeutic.name

        return ";".join(sorted([t.name for t in therapeutic.therapies]))

    @staticmethod
    def get_comment(
        record: VariantTherapeuticResponseStudyStatement
        | VariantDiagnosticStudyStatement
        | VariantPrognosticStudyStatement,
    ) -> str | None:
        """Get comment from a study statement

        :param record: GKS study statement
            Assumes ``name`` is provided in ``MappableConcept`` objects.
        :return: Comment for a given study statement.
            If the therapeutic is a substitute group, the original comment will be
            updated to make note that these are in substitution (deviating from clinvar
            api schema which notes that the therapies are in combination)
        """
        comment = record.description
        if isinstance(record, VariantTherapeuticResponseStudyStatement):
            therapeutic = record.proposition.objectTherapeutic.root
            if (
                isinstance(therapeutic, TherapyGroup)
                and therapeutic.membershipOperator == MembershipOperator.OR
            ):
                comment = f"{record.description or ''} NOTE: These therapies are in substitution.".strip()
        return comment

    @staticmethod
    def get_condition_set(
        proposition: VariantTherapeuticResponseProposition
        | VariantDiagnosticProposition
        | VariantPrognosticProposition,
    ) -> SubmissionConditionSetSomatic:
        """Get condition from a proposition

        We will prioritize sending DB/IDs over names

        :param proposition: Proposition for a given study statement.
            Assumes a single condition is used for the study statement and that
            ``name`` is provided in ``MappableConcept`` objects.
        :return: The condition for which the variant is interpreted
        """
        condition = (
            proposition.conditionQualifier
            if isinstance(proposition, VariantTherapeuticResponseProposition)
            else proposition.objectCondition
        )

        condition_db = None
        condition_id = None
        condition_name = None

        if condition_primary_coding := condition.root.primaryCoding:
            condition_primary_coding_id = condition_primary_coding.id
            condition_primary_coding_code = condition_primary_coding.code.root

            if condition_primary_coding_id.startswith(f"{ConditionDb.MONDO.value}_"):
                condition_db = ConditionDb.MONDO
                condition_id = condition_primary_coding_code
            elif condition_primary_coding_id.startswith("MIM:"):
                condition_db = ConditionDb.OMIM
                condition_id = condition_primary_coding_code
            elif condition_primary_coding_id.startswith(
                f"{ConditionDb.MEDGEN.value.lower()}:"
            ):
                condition_db = ConditionDb.MEDGEN
                condition_id = condition_primary_coding_code
            elif condition_primary_coding_id.startswith(
                f"{ConditionDb.ORPHANET.value.lower()}:"
            ):
                condition_db = ConditionDb.ORPHANET
                condition_id = f"ORPHA{condition_primary_coding_code}"
            elif condition_primary_coding_id.startswith(
                f"{ConditionDb.MESH.value.lower()}:"
            ):
                condition_db = ConditionDb.MESH
                condition_id = condition_primary_coding_code
            elif condition_primary_coding_code.startswith(f"{ConditionDb.HP.value}:"):
                condition_db = ConditionDb.HP
                condition_id = condition_primary_coding_code

            if not condition_db and not condition_id:
                condition_name = condition_primary_coding.name

        if not condition_db and not condition_id:
            condition_name = condition.root.name

        return SubmissionConditionSetSomatic(
            condition=[
                SubmissionCondition(
                    name=condition_name, db=condition_db, id=condition_id
                )
            ]
        )

    @staticmethod
    def get_clinical_impact_classification_description(
        record: VariantTherapeuticResponseStudyStatement
        | VariantDiagnosticStudyStatement
        | VariantPrognosticStudyStatement,
    ) -> SomaticClinicalImpactClassificationDescription:
        """Get AMP/ASCO/CAP classification

        :param record: GKS study statement
            Assumes ``classification`` uses ``primaryCode``
        :return: _description_
        """
        return SomaticClinicalImpactClassificationDescription(
            record.classification.primaryCode.root
        )

    def get_variant_set(
        self,
        proposition: VariantTherapeuticResponseProposition
        | VariantDiagnosticProposition
        | VariantPrognosticProposition,
        variant_hgvs: str | None = None,
    ) -> SubmissionVariantSet:
        """Get variant set

        This assumes only a single submission variant.

        :param proposition: Proposition for a given study statement.
        :param variant_hgvs: The HGVS expression for a variant, if found.
        :return: Variant set for a proposition
        """
        return SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs=variant_hgvs,
                    gene=[
                        SubmissionVariantGene(
                            symbol=proposition.geneContextQualifier.name
                        )
                    ],
                    alternate_designations=self.get_variant_aliases(
                        proposition.subjectVariant
                    ),
                )
            ]
        )

    def get_assertion_type_for_clinical_impact(
        self,
        proposition: VariantTherapeuticResponseProposition
        | VariantDiagnosticProposition
        | VariantPrognosticProposition,
    ) -> SomaticClinicalImpactAssertionType:
        """Get assertion type for clinical impact for a given proposition

        :param proposition: Proposition for a given study statement.
        :return: Assertion type for clinical impact
        """
        return self.gks_predicate_to_assertion[proposition.predicate]
