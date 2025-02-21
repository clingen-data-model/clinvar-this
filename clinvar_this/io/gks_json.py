"""Support for I/O of the GKS JSON format to define submissions.

Currently only supports CIViC Therapeutic Assertions. This assumes you are using
MetaKB (https://github.com/cancervariants/metakb) to generate GKS JSON files.
"""

import datetime
import json
import typing

from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import VariantTherapeuticResponseStudyStatement
from ga4gh.va_spec.base import TherapeuticResponsePredicate, EvidenceLine, TherapyGroup
from ga4gh.vrs.models import Syntax, Expression, Variation
from pydantic import ConfigDict, ValidationError
import requests

from clinvar_api.models import (
    AffectedStatus,
    AlleleOrigin,
    Assembly,
    Chromosome,
    CitationDb,
    CollectionMethod,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionChromosomeCoordinates,
    SubmissionCitation,
    SubmissionCondition,
    SubmissionContainer,
    SubmissionClinicalImpactSubmission,
    SubmissionVariant,
    SubmissionVariantSet,
)
from clinvar_api.models.sub_payload import (
    SomaticClinicalImpactClassification,
    SubmissionConditionSetSomatic,
    SubmissionObservedInSomatic,
    SubmissionVariantGene,
)
from clinvar_api.msg.sub_payload import (
    SomaticClinicalImpactAssertionType,
    SomaticClinicalImpactClassificationDescription,
)
from clinvar_this import exceptions
from clinvar_this.io.base import TransformIO, BatchMetadata


class BatchMetadataGksJson(BatchMetadata):
    """Batch-wide settings for GKS JSON import.

    The properties will be assigned to all variants/samples in the batch.
    """

    model_config = ConfigDict(frozen=True)
    affected_status: typing.Optional[AffectedStatus] = None


class GksJsonTransformer(TransformIO):
    """Class for transforming GKS JSON input data from various formats into submission format"""

    batch_metadata_defaults = BatchMetadataGksJson(
        collection_method=CollectionMethod.CURATION,
        allele_origin=AlleleOrigin.SOMATIC,
        affected_status=AffectedStatus.UNKNOWN,
    )

    # Mapping from GKS classification to ClinVar
    gks_class_to_impact_class = {
        "Tier I": SomaticClinicalImpactClassificationDescription.STRONG,
        "Tier II": SomaticClinicalImpactClassificationDescription.POTENTIAL,
        "Tier III": SomaticClinicalImpactClassificationDescription.UNKNOWN,
        "Tier IV": SomaticClinicalImpactClassificationDescription.BENIGN_LIKELY_BENIGN,
    }

    gks_predicate_to_assertion = {
        TherapeuticResponsePredicate.RESISTANCE: SomaticClinicalImpactAssertionType.THERAPEUTIC_RESISTANCE,
        TherapeuticResponsePredicate.SENSITIVITY: SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
    }

    @staticmethod
    def _read_file(
        inputf: typing.TextIO,
    ) -> typing.List[VariantTherapeuticResponseStudyStatement]:
        """Get list of CIViC Variant Therapeutic Response Study Statements from a file

        For now, diagnostic and prognostic CIViC assertions are not supported

        :param inputf: Text file-like object containing input GKS Statement data
        :raises exceptions.InvalidFormat: If there was an error decoding JSON
        :return: A list of CIViC Assertions represented as Variant Therapeutic Response
            Study Statements
        """
        tr_assertions: typing.List[dict] = []

        try:
            civic_assertions = json.load(inputf)
        except json.JSONDecodeError as e:
            err_msg = "Error decoding JSON"
            raise exceptions.InvalidFormat(err_msg) from e

        for civic_assertion in civic_assertions:
            try:
                tr_assertions.append(
                    VariantTherapeuticResponseStudyStatement(**civic_assertion)
                )
            except ValidationError:
                pass

        return tr_assertions

    @staticmethod
    def _get_variant_hgvs(variant: Variation) -> str | None:
        """Get HGVS expression for a variant

        :param variant: Variant record from CIViC
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

        if getattr(variant, "constraints", None):
            expressions = variant.constraints[0].root.allele.expressions
        else:
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
    def _get_variant_coords(
        variant: Variation,
    ) -> SubmissionChromosomeCoordinates | None:
        """Get the chromosome coordinates for a variant

        Looks in the variant's extensions to find an extension with the name
        'CIViC representative coordinate'

        :param variant: Variant record from CIViC
        :return: Chromosome coordinates for a variant, if found
        """
        try:
            rep_coords_ext = next(
                ext
                for ext in variant.extensions
                if ext.name == "CIViC representative coordinate"
            ).value
        except StopIteration:
            return

        return SubmissionChromosomeCoordinates(
            assembly=Assembly(rep_coords_ext["reference_build"]),
            chromosome=Chromosome(rep_coords_ext["chromosome"]),
            start=rep_coords_ext["start"],
            stop=rep_coords_ext["stop"],
        )

    @staticmethod
    def _get_citations(
        mp_id: str,
        evidence_lines: typing.List[EvidenceLine],
    ) -> typing.List[SubmissionCitation]:
        """Get list of PubMed citations and supporting CIViC evidence for an assertion

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
                            SubmissionCitation(db=CitationDb.PUBMED, id=str(doc.pmid))
                        )

        return citations

    @staticmethod
    def _get_drugs(therapeutic: TherapyGroup | MappableConcept) -> str:
        """Get the name for drug(s)

        :param therapeutic: CIViC therapeutic record
        :return: The formatted name for a CIViC therapeutic record. Multiple therapies
            for combination and substitution will be separated by semicolons.
        """
        if isinstance(therapeutic, MappableConcept):
            return therapeutic.name

        return ";".join([t.name for t in therapeutic.therapies])

    @staticmethod
    def _get_submitted_date(record_id: int) -> str | None:
        """Get the date (yyyy-mm-dd) for the submission for an assertion

        Since CIViCPy does not include event information, use the CIVIC GraphQL API

        :param record_id: The ID for the CIViC assertion
        :return: The date (yyyy-mm-dd) for the submission for an assertion
        """
        query = f"""
            {{
                assertion(id: {record_id}) {{
                    submissionEvent {{
                        createdAt
                    }}
                }}
            }}
        """

        resp = requests.post(
            "https://civicdb.org/api/graphql",
            json={"query": query},
            headers={"Content-Type": "application/json"},
        )

        if resp.status_code == 200:
            resp = resp.json()
            submitted_date = (
                resp.get("data", {})
                .get("assertion", {})
                .get("submissionEvent", {})
                .get("createdAt", "")
            )
            submitted_date = datetime.datetime.strptime(
                submitted_date, "%Y-%m-%dT%H:%M:%SZ"
            )
            submitted_date = submitted_date.strftime("%Y-%m-%d")
            return submitted_date
        return

    def _get_clinical_impact_submission(
        self,
        record: VariantTherapeuticResponseStudyStatement,
        variant_hgvs: str | None = None,
        variant_coords: SubmissionChromosomeCoordinates | None = None,
    ) -> SubmissionClinicalImpactSubmission:
        """Get clinical impact submission for a CIViC therapeutic assertion

        Assertions with Substitutes therapies will be separated by semicolons
        and the CIViC assertion's description will be updated to include this note.

        Allele origin, collection method, and affected status will always use the
        ``batch_metadata_defaults``

        :param record: The CIViC therapeutic assertion
        :param variant_hgvs: The HGVS expression for a variant, if found. If not found,
            ``variant_coords`` must be provided. This takes priority over
            ``variant_coords``.
        :param variant_coords: The chromosome coordinates for a variant, if found. If
            not found, ``variant_hgvs`` must be provided
        :return: The clinical impact submission for a CIViC therapeutic assertion
        """
        proposition = record.proposition
        therapeutic = proposition.objectTherapeutic.root

        if (
            isinstance(therapeutic, TherapyGroup)
            and therapeutic.groupType.name == "TherapeuticSubstituteGroup"
        ):
            description = (
                f"{record.description} NOTE: These therapies are in substitution."
            )
        else:
            description = record.description

        mp_id = proposition.subjectVariant.id

        return SubmissionClinicalImpactSubmission(
            record_status=RecordStatus.NOVEL,
            local_id=mp_id,
            local_key=record.id,
            observed_in=[
                SubmissionObservedInSomatic(
                    allele_origin=self.batch_metadata_defaults.allele_origin,
                    collection_method=self.batch_metadata_defaults.collection_method,
                    affected_status=self.batch_metadata_defaults.affected_status,
                )
            ],
            condition_set=SubmissionConditionSetSomatic(
                condition=[
                    SubmissionCondition(
                        name=proposition.conditionQualifier.root.name,
                    )
                ]
            ),
            variant_set=SubmissionVariantSet(
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
                    )
                ]
            ),
            clinical_impact_classification=SomaticClinicalImpactClassification(
                clinical_impact_classification_description=self.gks_class_to_impact_class[
                    record.classification.primaryCode.root
                ],
                assertion_type_for_clinical_impact=self.gks_predicate_to_assertion[
                    proposition.predicate
                ],
                comment=description,
                citation=self._get_citations(mp_id, record.hasEvidenceLines),
                drug_for_therapeutic_assertion=self._get_drugs(therapeutic),
                date_last_evaluated=self._get_submitted_date(
                    int(record.id.split("civic.aid:")[-1])
                ),
            ),
        )

    def records_to_submission_container(
        self,
        civic_tr_assertions: typing.List[VariantTherapeuticResponseStudyStatement],
    ) -> typing.List[SubmissionContainer]:
        """Transform GKS records to submission container data structures

        Will only submit using clinical impact submissions

        :param civic_tr_assertions: List of CIViC Therapeutic Assertions represented
            as GKS Variant Therapeutic Response Study Statements
        :return: A list of submission container data structures
        """
        clinical_impact_submissions: typing.List[SubmissionContainer] = []

        for civic_tr_assertion in civic_tr_assertions:
            variant = civic_tr_assertion.proposition.subjectVariant
            variant_hgvs = self._get_variant_hgvs(variant)
            variant_coords = None
            if not variant_hgvs:
                variant_coords = self._get_variant_coords(variant)
                if not variant_coords:
                    continue

            clinical_impact_submission = self._get_clinical_impact_submission(
                civic_tr_assertion,
                variant_hgvs=variant_hgvs,
                variant_coords=variant_coords,
            )
            clinical_impact_submissions.append(clinical_impact_submission)

        return SubmissionContainer(
            assertion_criteria=SubmissionAssertionCriteria(
                db=CitationDb.PUBMED,
                id="27993330",  # AMP/ASCO/CAP
            ),
            clinical_impact_submission=clinical_impact_submissions,
        )
