"""Support for I/O of the GKS JSON format to define submissions.

Currently only supports Therapeutic, Diagnostic, and Prognostic Assertions that follow
AMP/ASCO/CAP 2017 guidelines.
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
    Assembly,
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

    # Mapping from GKS predicate type to ClinVar assertion type for clinical impact
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
        statements: typing.List[
            VariantTherapeuticResponseStudyStatement
            | VariantDiagnosticStudyStatement
            | VariantPrognosticStudyStatement
        ],
        batch_metadata: BatchMetadata,
    ) -> SubmissionContainer:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param statements: List of GKS statements
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
        """Get list of Variant Therapeutic Response, Diagnostic, or Prognostic Statements
        from a file

        Expects `gks_records` key to contain list of GKS formatted statements

        :param inputf: Text file-like object containing input GKS Statement data
        :raise exceptions.InvalidFormat: If there was an error decoding JSON
        :raise KeyError: If JSON is missing `gks_records` key
        :return: A list of Variant Therapeutic Response, Diagnostic, or Prognostic
            Statements
        """
        statements: typing.List[
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

        gks_records = data["gks_records"]

        supported_stmts: typing.List[callable] = [
            VariantTherapeuticResponseStudyStatement,
            VariantDiagnosticStudyStatement,
            VariantPrognosticStudyStatement,
        ]

        for gks_record in gks_records:
            for supported_stmt in supported_stmts:
                try:
                    statements.append(supported_stmt(**gks_record))
                except ValidationError:
                    pass

        return statements

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
        """Get comment from a statement

        :param record: GKS statement
            Assumes ``name`` is provided in ``MappableConcept`` objects.
        :return: Comment for a given statement.
            If the therapeutic is a substitute group, the original comment will be
            updated to make note that these are in substitution (deviating from ClinVar
            API schema which notes that the therapies are in combination)
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
        proposition: (
            VariantTherapeuticResponseProposition
            | VariantDiagnosticProposition
            | VariantPrognosticProposition
        ),
    ) -> SubmissionConditionSetSomatic:
        """Build a somatic condition set from a proposition.

        Database identifiers are preferred over condition names when supported
        ontology mappings are available.

        Assumes:
        - A single condition is associated with the proposition
        - `MappableConcept.name` is populated when no supported mapping exists

        :param proposition: Proposition for a given statement
        :return: Condition set for the interpreted variant.
        """

        condition = (
            proposition.conditionQualifier
            if isinstance(proposition, VariantTherapeuticResponseProposition)
            else proposition.objectCondition
        )

        condition_db: ConditionDb | None = None
        condition_id: str | None = None
        condition_name: str | None = None

        if primary_coding := condition.root.primaryCoding:
            coding_id = primary_coding.id
            coding_code = primary_coding.code.root

            # Mapping of supported condition database prefixes to:
            # (coding ID prefix, submission database enum, transformed submission identifier)
            db_mappings: tuple[
                tuple[str, ConditionDb, str],
                ...,
            ] = (
                (f"{ConditionDb.MONDO.value}_", ConditionDb.MONDO, coding_code),
                ("MIM:", ConditionDb.OMIM, coding_code),
                (
                    f"{ConditionDb.MEDGEN.value.lower()}:",
                    ConditionDb.MEDGEN,
                    coding_code,
                ),
                (
                    f"{ConditionDb.ORPHANET.value.lower()}:",
                    ConditionDb.ORPHANET,
                    f"ORPHA{coding_code}",
                ),
                (
                    f"{ConditionDb.MESH.value.lower()}:",
                    ConditionDb.MESH,
                    coding_code,
                ),
            )

            for coding_id_prefix, db, mapped_id in db_mappings:
                if coding_id.startswith(coding_id_prefix):
                    condition_db = db
                    condition_id = mapped_id
                    break

            if coding_code.startswith(f"{ConditionDb.HP.value}:"):
                condition_db = ConditionDb.HP
                condition_id = coding_code

            if not condition_db and not condition_id:
                condition_name = primary_coding.name

        if not condition_db and not condition_id:
            condition_name = condition.root.name

        return SubmissionConditionSetSomatic(
            condition=[
                SubmissionCondition(
                    name=condition_name,
                    db=condition_db,
                    id=condition_id,
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

        :param record: GKS statement
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

        :param proposition: Proposition for a given statement.
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

        :param proposition: Proposition for a given statement.
        :return: Assertion type for clinical impact
        """
        return self.gks_predicate_to_assertion[proposition.predicate]
