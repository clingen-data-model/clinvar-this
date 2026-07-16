"""Support for I/O of the GKS JSON format to define ClinVar submissions.

Currently only supports:
* Therapeutic, Diagnostic, and Prognostic Statements that follow
AMP/ASCO/CAP 2017 guidelines. These will map to `clinical_impact_submission`.
* Oncogenic Statements that follow ClinGen/CGC/VICC 2022 guidelines. These will map to `oncogenicity_submission`.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38
"""

import json
import re
from abc import ABC, abstractmethod
from typing import Any, Generic, Iterable, Literal, TextIO, TypeVar, get_type_hints

from ga4gh.cat_vrs.models import CategoricalVariant, DefiningAlleleConstraint
from ga4gh.core.models import MappableConcept, iriReference
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
from ga4gh.va_spec.base import (
    Contribution,
    Document,
    EvidenceLine,
    MembershipOperator,
    Statement,
    StudyResult,
    TherapyGroup,
    VariantDiagnosticProposition,
    VariantOncogenicityProposition,
    VariantPrognosticProposition,
    VariantTherapeuticResponseProposition,
)
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from ga4gh.vrs.models import Allele, Expression, MolecularVariation, Syntax
from logzero import logfile, logger
from pydantic import BaseModel, ConfigDict

from clinvar_api.models import (
    AffectedStatus,
    AlleleOrigin,
    Assembly,
    CitationDb,
    CollectionMethod,
    ConditionDb,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionCitation,
    SubmissionClinicalImpactSubmission,
    SubmissionCondition,
    SubmissionConditionSetSomatic,
    SubmissionContainer,
    SubmissionObservedInSomatic,
    SubmissionOncogenicitySubmission,
    SubmissionVariant,
    SubmissionVariantGene,
    SubmissionVariantSet,
)
from clinvar_this import exceptions
from clinvar_this.io.base import TransformIO

logfile("gks_json.log")

# Supported GKS Statement types for ClinVar Submission
GksStatementT = TypeVar(
    "GksStatementT",
    VariantClinicalSignificanceStatement,
    VariantOncogenicityStatement,
)

# Name of ClinVar Submission Container attribute
SubmissionContainerAttribute = Literal[
    "clinical_impact_submission",
    "oncogenicity_submission",
]

# Supported ClinVar Submission types
SubmissionT = TypeVar(
    "SubmissionT",
    SubmissionClinicalImpactSubmission,
    SubmissionOncogenicitySubmission,
)


CLINVAR_ACCESSION_RE = re.compile(r"^[A-Z]{3}\d{9}(\.\d+)?$")


class BatchMetadata(BaseModel):
    """Batch-wide settings for GKS JSON import.

    The properties will be assigned to all variants/samples in the batch.
    """

    model_config = ConfigDict(frozen=True)

    affected_status: AffectedStatus = AffectedStatus.UNKNOWN
    collection_method: CollectionMethod = CollectionMethod.NOT_PROVIDED
    submitted_assembly: Assembly = Assembly.GRCH38


def batch_metadata_from_mapping(
    keys_values: Iterable[str],
) -> BatchMetadata:
    """Convert configuration from ``KEY=VALUE`` strings to ``BatchMetadata``

    If values are not provided, then will use defaults
    """
    field_types = {
        name: value for (name, value) in get_type_hints(BatchMetadata).items()
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


class GksJsonTransformer(TransformIO, ABC, Generic[GksStatementT]):
    """Abstract base class for transforming GKS JSON input data from various formats into submission format"""

    submission_container_attribute: SubmissionContainerAttribute
    assertion_criteria: SubmissionAssertionCriteria
    gks_statement_cls: type[GksStatementT]

    @classmethod
    def _read_file(
        cls,
        inputf: TextIO,
    ) -> list[GksStatementT]:
        """Read GKS statements from a JSON file object.

        Expects the JSON file to contain a ``gks_records`` key with a list of GKS
        statement records compatible with the transformer's configured
        ``gks_statement_cls``.

        :param inputf: Text file object containing GKS JSON data
        :raise exceptions.InvalidFormat: If the input is not valid JSON or does not
            contain the required ``gks_records`` key
        :raise ValueError: If one or more GKS records does not validate to the
            configured GKS statement class
        :return: Parsed list of GKS statements
        """
        statements: list[GksStatementT] = []

        try:
            data = json.load(inputf)
        except json.JSONDecodeError as e:
            err_msg = "Error decoding GKS JSON"
            raise exceptions.InvalidFormat(err_msg) from e

        if "gks_records" not in data:
            msg = "Invalid GKS JSON: missing required key `gks_records` (must be a list of statements)"
            raise exceptions.InvalidFormat(msg)

        gks_records = data["gks_records"]

        for i, gks_record in enumerate(gks_records):
            try:
                statements.append(cls.gks_statement_cls(**gks_record))
            except Exception as e:
                msg = f"Failed to validate GKS statement at index {i} into {cls.gks_statement_cls.__name__}: statement.id={gks_record.get('id')}"
                raise ValueError(msg) from e

        return statements

    @staticmethod
    def _get_local_id(variant: CategoricalVariant | Allele) -> str | None:
        """Get local ID

        :param variant: Variant associated to statement
        :return: Variant ID or name, if exists
        """
        return variant.id or variant.name

    @staticmethod
    def _get_variant_hgvs(
        variant: Allele | CategoricalVariant,
    ) -> str | None:
        """Retrieve a HGVS expression for a variant

        For Categorical Variants, checks the first constraint for an expression. Only
        support extracting from DefiningAlleleConstraints at the moment.
        If no constraints found, then checks the expressions extension. These are cases
        where an HGVS expression is unable to be representing using VRS.

        For VRS Alleles, checks the `expressions` field.

        Order matters: the first matching expression is returned. cDNA RefSeq HGVS
        expressions are prioritized over genomic RefSeq HGVS expressions.

        :param variant: Variant associated to statement
        :return: cDNA RefSeq HGVS expression or genomic RefSeq HGVS expression for a
            variant, if provided.
        """

        def get_hgvs(
            expressions: list[Expression] | None,
            syntax: Literal[Syntax.HGVS_C, Syntax.HGVS_G],
        ) -> str | None:
            """Get a HGVS expression from a list of expressions for a given syntax

            :param expressions: List of representations specified by nomenclature or
                syntax for a variant
            :param syntax: The syntax to find an expression for
            :return: HGVS expression for a given syntax, if found
            """
            if not expressions:
                return None

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

        if isinstance(variant, CategoricalVariant):
            if getattr(variant, "constraints", None) and variant.constraints:
                constraint = variant.constraints[0]
                if isinstance(constraint.root, DefiningAlleleConstraint) and isinstance(
                    constraint.root.allele, Allele
                ):
                    expressions = constraint.root.allele.expressions
            else:
                # Case where a VRS Allele is unable to be represented, can store
                # expressions as an extension named 'expressions'
                try:
                    expressions_ext = next(
                        ext
                        for ext in variant.extensions or []
                        if ext.name == "expressions"
                    ).value
                    if isinstance(expressions_ext, list):
                        expressions = [Expression(**ext) for ext in expressions_ext]
                except (StopIteration, TypeError):
                    return None
        else:
            expressions = variant.expressions

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
    def _get_comment(
        description: str | None,
        therapeutic: TherapyGroup | MappableConcept | None,
    ) -> str | None:
        """Get comment from a statement

        :param description: Description for GKS statement
        :param therapeutic: Therapeutic for GKS statement, if one exists
        :return: Comment for a given statement.
            If the therapeutic is a substitute group, the original comment will be
            updated to make note that these are in substitution (deviating from ClinVar
            API schema which notes that the therapies are in combination)
        """
        if therapeutic and (
            isinstance(therapeutic, TherapyGroup)
            and therapeutic.membershipOperator == MembershipOperator.OR
        ):
            description = f"{description or ''} NOTE: These therapies are in substitution.".strip()
        return description

    @staticmethod
    def _get_condition_set(
        proposition: (
            VariantTherapeuticResponseProposition
            | VariantDiagnosticProposition
            | VariantPrognosticProposition
            | VariantOncogenicityProposition
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

        if isinstance(proposition, VariantTherapeuticResponseProposition):
            condition = proposition.conditionQualifier
        elif isinstance(proposition, VariantOncogenicityProposition):
            condition = proposition.objectTumorType
        else:
            condition = proposition.objectCondition

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
    def _get_citations(
        reported_in: list[Document | iriReference] | None,
        evidence_lines: list[EvidenceLine],
    ) -> list[SubmissionCitation]:
        """Extract unique citations from statement documents, evidence lines,
        and related evidence items.

        Citations may be sourced from:
        - Statement documents (reported_in field)
        - PubMed IDs (`pmid`) on reported documents
        - DOI IDs (`doi`) on reported documents
        - IRI reference URLs attached to reported documents

        Each citation is submitted using the most specific reference available:
        PMID, DOI, or URL.

        :param reported_in: Documents in which the statement is reported
        :param evidence_lines: Evidence lines that may contain supporting citation
            information.
        :return: A deduplicated list of submission citations.
        """

        def add_citation(citation: SubmissionCitation) -> None:
            """Append a citation if it has not already been added."""
            if citation not in citations:
                citations.append(citation)

        def add_reference(
            reference: str | iriReference | Document | None,
        ) -> None:
            """Convert a reference into a SubmissionCitation."""
            if reference is None:
                return

            if isinstance(reference, str):
                add_citation(SubmissionCitation(url=reference))
                return

            if isinstance(reference, iriReference):
                add_citation(SubmissionCitation(url=reference.root))
                return

            if isinstance(reference, Document):
                if reference.pmid:
                    add_citation(
                        SubmissionCitation(
                            db=CitationDb.PUBMED,
                            id=reference.pmid,
                        )
                    )
                elif reference.doi:
                    add_citation(
                        SubmissionCitation(
                            db=CitationDb.DOI,
                            id=reference.doi,
                        )
                    )

        citations: list[SubmissionCitation] = []
        reported_in_documents: list[Document | iriReference] = reported_in or []

        for evidence_line in evidence_lines:
            reported_in_documents.extend(evidence_line.reportedIn or [])

            for evidence_item in evidence_line.hasEvidenceItems or []:
                reported_in = None
                if isinstance(evidence_item, Statement | EvidenceLine):
                    reported_in = evidence_item.reportedIn
                elif isinstance(evidence_item, StudyResult):
                    reported_in = evidence_item.root.reportedIn

                reported_in_documents.extend(reported_in or [])

        for document in reported_in_documents:
            add_reference(document)

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

        if contribution.date:
            return contribution.date.strftime("%Y-%m-%d")

    def _build_shared_submission_kwargs(
        self,
        statement: GksStatementT,
        observed_in: list[SubmissionObservedInSomatic],
        variant: CategoricalVariant | Allele,
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ) -> dict[str, Any]:
        """Build submission kwargs that are used in both clinical impact and
        oncogenicity submissions

        :param statement: GKS statement
        :param observed_in: List of distinct observations
        :param variant: Variant associated to statement
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: Dictionary containing shared submission kwargs to be passed to the
            Submission container
        """
        proposition = statement.proposition

        clinvar_accession = next(
            (
                str(extension.value)
                for extension in statement.extensions or []
                if extension.name == "clinvar_accession"
            ),
            None,
        )

        if clinvar_accession and not CLINVAR_ACCESSION_RE.match(clinvar_accession):
            logger.warning(
                "Statement ID %s ClinVar accession %s does not match ClinVar regex",
                statement.id,
                clinvar_accession,
            )
            clinvar_accession = None

        record_status = (
            RecordStatus.NOVEL if clinvar_accession is None else RecordStatus.UPDATE
        )

        return {
            "clinvar_accession": clinvar_accession,
            "record_status": record_status,
            "local_id": self._get_local_id(variant),
            "submitted_assembly": submitted_assembly,
            "local_key": statement.id,
            "observed_in": observed_in,
            "condition_set": self._get_condition_set(proposition),
            "variant_set": self._get_variant_set(
                proposition.geneContextQualifier,
                variant,
                variant_hgvs=variant_hgvs,
            ),
        }

    def _build_shared_classification_kwargs(
        self,
        reported_in: list[Document | iriReference] | None,
        description: str | None,
        therapeutic: TherapyGroup | MappableConcept | None,
        evidence_lines: list[EvidenceLine] | None,
        contributions: list[Contribution] | None,
    ):
        """Build classification kwargs that are used in both clinical impact and
        oncogenicity classifications

        :param reported_in: Documents in which the statement is reported
        :param description: Description for GKS statement
        :param therapeutic: Therapeutic for GKS statement
        :param evidence_lines: Evidence lines associated to GKS statement
        :param contributions: Contributions for GKS statement
        :return: Dictionary containing shared classification kwargs to be passed to the
            classification container
        """
        return {
            "comment": self._get_comment(description, therapeutic),
            "citation": self._get_citations(reported_in, evidence_lines or []),
            "date_last_evaluated": self._get_date_last_evaluated(contributions or []),
        }

    def _get_variant_set(
        self,
        gene_context: MappableConcept,
        variant: CategoricalVariant | Allele,
        variant_hgvs: str | None = None,
    ) -> SubmissionVariantSet:
        """Get variant set

        This assumes only a single submission variant.

        :param gene_context: Gene associated to statement.
        :param variant: Variant associated to statement.
        :param variant_hgvs: The HGVS expression for a variant, if found.
        :return: Variant set for a proposition
        """
        return SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs=variant_hgvs,
                    gene=[SubmissionVariantGene(symbol=gene_context.name)],
                    alternate_designations=variant.aliases,
                )
            ]
        )

    def _get_method_type(
        self,
        statement: GksStatementT,
        batch_metadata: BatchMetadata,
    ) -> CollectionMethod:
        """Get the collection method associated with a GKS statement.

        Attempts to derive the collection method from the statement's
        ``specifiedBy.methodType`` field. If the method type is missing or
        cannot be converted into a supported ``CollectionMethod``, the
        batch-level default collection method is used instead.

        :param statement: GKS statement containing optional method type information
        :param batch_metadata: Batch-wide metadata containing fallback collection method settings
        :return: Collection method associated with the statement
        """
        try:
            if stmt_method_type := statement.specifiedBy.methodType:
                return CollectionMethod(stmt_method_type)
        except ValueError:
            pass

        return batch_metadata.collection_method

    @staticmethod
    def _get_evidence_lines(els: list[EvidenceLine | iriReference] | None):
        els = els or []
        return [el for el in els if isinstance(el, EvidenceLine)]

    @abstractmethod
    def _get_submission(
        self,
        statement: GksStatementT,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ) -> SubmissionT:
        """Transform a GKS statement into a ClinVar novel submission

        Local ID will use the proposition's variant ID or name.

        Local Key will use the `record`'s ID.

        If `clinvar_accession` extension exists in `statement`, then this variant will
        have record status as `update` rather than `novel`.

        :param statement: GKS statement instance to transform
        :param observed_in: List of distinct ClinVar somatic observations associated
            with the statement
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: ClinVar submission model corresponding to the transformer type
        """

    def records_to_submission_container(
        self,
        statements: list[GksStatementT],
        batch_metadata: BatchMetadata,
    ) -> SubmissionContainer:
        """Transform GKS statements into a ClinVar submission container.

        Converts all provided GKS statements into ClinVar submission models
        and packages them into a ``SubmissionContainer``.

        Statements that cannot be transformed due to missing required data
        (for example HGVS or ``ObservedIn`` information) are skipped.

        :param statements: List of parsed GKS statements to transform
        :param batch_metadata: Batch-wide metadata applied to all generated submissions
        :return: ClinVar submission container containing transformed submissions
        """

        submissions = []

        for statement in statements:
            statement_id = statement.id

            variant: Allele | CategoricalVariant
            prop_variant = statement.proposition.subjectVariant

            if isinstance(prop_variant, MolecularVariation):
                if not isinstance(prop_variant.root, Allele):
                    logger.warning(
                        "Skipping statement. Molecular Variation is not an Allele for statement ID: %s",
                        statement_id,
                    )
                    continue

                variant = prop_variant.root
            elif isinstance(prop_variant, CategoricalVariant):
                variant = prop_variant
            else:
                logger.warning(
                    "Skipping statement. Variant is not a Categorical Variant or Allele for statement ID: %s",
                    statement_id,
                )
                continue

            variant_hgvs = self._get_variant_hgvs(variant)

            if not variant_hgvs:
                logger.warning(
                    "Skipping statement. No HGVS found for statement ID: %s",
                    statement_id,
                )
                continue

            method_type = self._get_method_type(statement, batch_metadata)

            observed_in = self._get_observed_in(
                statement.proposition.alleleOriginQualifier,
                method_type,
                batch_metadata.affected_status,
            )

            if not observed_in:
                logger.warning(
                    "Skipping statement. No observed_in found for statement ID: %s",
                    statement_id,
                )
                continue

            submissions.append(
                self._get_submission(
                    statement,
                    observed_in,
                    variant,
                    variant_hgvs=variant_hgvs,
                    submitted_assembly=batch_metadata.submitted_assembly,
                )
            )

        return SubmissionContainer(
            assertion_criteria=self.assertion_criteria,
            **{self.submission_container_attribute: submissions},
        )
