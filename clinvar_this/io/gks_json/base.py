"""Support for I/O of the GKS JSON format to define submissions.

Currently only supports Therapeutic Assertions.
"""

from abc import ABC, abstractmethod
import json
import typing

from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import VariantTherapeuticResponseStudyStatement
from ga4gh.va_spec.base import (
    TherapeuticResponsePredicate,
    TherapyGroup,
    VariantTherapeuticResponseProposition,
)
from ga4gh.vrs.models import MolecularVariation
from pydantic import ValidationError


from clinvar_api.models import (
    SubmissionChromosomeCoordinates,
    SubmissionCondition,
    SubmissionContainer,
    SubmissionVariant,
    SubmissionVariantSet,
)
from clinvar_api.models.sub_payload import (
    SubmissionConditionSetSomatic,
    SubmissionVariantGene,
)
from clinvar_api.msg.sub_payload import (
    SomaticClinicalImpactAssertionType,
    SomaticClinicalImpactClassificationDescription,
)
from clinvar_this import exceptions
from clinvar_this.io.base import TransformIO


class GksJsonTransformer(TransformIO, ABC):
    """Class for transforming GKS JSON input data from various formats into submission format"""

    # Mapping from GKS predicate type to assertion type for clinical impact
    gks_predicate_to_assertion = {
        TherapeuticResponsePredicate.RESISTANCE: SomaticClinicalImpactAssertionType.THERAPEUTIC_RESISTANCE,
        TherapeuticResponsePredicate.SENSITIVITY: SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
    }

    @abstractmethod
    def records_to_submission_container(
        self,
        study_statements: typing.List[VariantTherapeuticResponseStudyStatement],
    ) -> typing.List[SubmissionContainer]:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param study_statements: List of GKS study statements
        :return: A list of submission container data structures
        """

    @staticmethod
    def _read_file(
        inputf: typing.TextIO,
    ) -> typing.List[VariantTherapeuticResponseStudyStatement]:
        """Get list of Variant Therapeutic Response Study Statements from a file

        For now, diagnostic and prognostic study statements are NOT supported

        :param inputf: Text file-like object containing input GKS Statement data
        :raises exceptions.InvalidFormat: If there was an error decoding JSON
        :return: A list of Variant Therapeutic Response Study Statements
        """
        study_statements: typing.List[VariantTherapeuticResponseStudyStatement] = []

        try:
            statements = json.load(inputf)
        except json.JSONDecodeError as e:
            err_msg = "Error decoding JSON"
            raise exceptions.InvalidFormat(err_msg) from e

        for statement in statements:
            try:
                study_statements.append(
                    VariantTherapeuticResponseStudyStatement(**statement)
                )
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

        return ";".join([t.name for t in therapeutic.therapies])

    @staticmethod
    def get_comment(record: VariantTherapeuticResponseStudyStatement) -> str | None:
        """Get comment from a study statement

        :param record: GKS study statement
            Assumes ``name`` is provided in ``MappableConcept`` objects.
        :return: Comment for a given study statement.
            If the therapeutic is a substitute group, the original comment will be
            updated to make note that these are in substitution (deviating from clinvar
            api schema which notes that the therapies are in combination)
        """
        therapeutic = record.proposition.objectTherapeutic.root
        if (
            isinstance(therapeutic, TherapyGroup)
            and therapeutic.groupType.name == "TherapeuticSubstituteGroup"
        ):
            comment = f"{record.description or ''} NOTE: These therapies are in substitution.".strip()
        else:
            comment = record.description
        return comment

    @staticmethod
    def get_condition_set(
        proposition: VariantTherapeuticResponseProposition,
    ) -> SubmissionConditionSetSomatic:
        """Get condition from a proposition

        :param proposition: Proposition for a given study statement.
            Assumes a single condition is used for the study statement and that
            ``name`` is provided in ``MappableConcept`` objects.
        :return: The condition for which the variant is interpreted
        """
        return SubmissionConditionSetSomatic(
            condition=[
                SubmissionCondition(
                    name=proposition.conditionQualifier.root.name,
                )
            ]
        )

    @staticmethod
    def get_clinical_impact_classification_description(
        record: VariantTherapeuticResponseStudyStatement,
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
        proposition: VariantTherapeuticResponseProposition,
        variant_hgvs: str | None = None,
        variant_coords: SubmissionChromosomeCoordinates | None = None,
    ) -> SubmissionVariantSet:
        """Get variant set

        This assumes only a single submission variant.

        :param proposition: Proposition for a given study statement.
        :param variant_hgvs: The HGVS expression for a variant, if found. If not found,
            ``variant_coords`` must be provided. This takes priority over
            ``variant_coords``.
        :param variant_coords: The chromosome coordinates for a variant, if found. If
            not found, ``variant_hgvs`` must be provided
        :return: Variant set for a proposition
        """
        return SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs=variant_hgvs,
                    chromosome_coordinates=variant_coords
                    if variant_hgvs is None
                    else None,
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
        self, proposition: VariantTherapeuticResponseProposition
    ) -> SomaticClinicalImpactAssertionType:
        """Get assertion type for clinical impact for a given proposition

        :param proposition: Proposition for a given study statement.
        :return: Assertion type for clinical impact
        """
        return self.gks_predicate_to_assertion[proposition.predicate]
