"""Support for I/O of the GKS JSON format to define ClinVar submissions.

Currently only supports:
* Therapeutic, Diagnostic, and Prognostic Statements that follow
AMP/ASCO/CAP 2017 guidelines. These will map to `clinical_impact_submission`.
* Oncogenic Statements that follow ClinGen/CGC/VICC 2022 guidelines. These will map to `oncogenicity_submission`.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38
"""

from abc import ABC, abstractmethod
import json
import typing

from ga4gh.cat_vrs.models import CategoricalVariant, DefiningAlleleConstraint
from ga4gh.core.models import MappableConcept
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
from ga4gh.va_spec.base import (
    VariantDiagnosticProposition,
    VariantPrognosticProposition,
    MembershipOperator,
    TherapyGroup,
    VariantTherapeuticResponseProposition,
    VariantOncogenicityProposition,
)
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from ga4gh.vrs.models import Expression, Syntax
from pydantic import BaseModel, ConfigDict


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
from clinvar_api.msg.sub_payload import ConditionDb, RecordStatus
from clinvar_this import exceptions
from clinvar_this.io.base import TransformIO
from ga4gh.core.models import iriReference
from ga4gh.va_spec.base import (
    Contribution,
    Document,
    EvidenceLine,
)

from clinvar_api.models import (
    CitationDb,
    SubmissionCitation,
)
from clinvar_api.models.sub_payload import (
    SubmissionObservedInSomatic,
)
from clinvar_api.msg.sub_payload import (
    AlleleOrigin,
)
from logzero import logger, logfile

from clinvar_api.models import (
    SubmissionAssertionCriteria,
    SubmissionClinicalImpactSubmission,
    SubmissionOncogenicitySubmission,
)

logfile("gks_json.log")

# Supported GKS Statement types for ClinVar Submission
GksStatementT = typing.TypeVar(
    "GksStatementT",
    VariantClinicalSignificanceStatement,
    VariantOncogenicityStatement,
)

# Name of ClinVar Submission Container attribute
SubmissionContainerAttribute = typing.Literal[
    "clinical_impact_submission",
    "oncogenicity_submission",
]

# Supported ClinVar Submission types
SubmissionT = typing.TypeVar(
    "SubmissionT",
    SubmissionClinicalImpactSubmission,
    SubmissionOncogenicitySubmission,
)


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


class GksJsonTransformer(TransformIO, ABC, typing.Generic[GksStatementT]):
    """Class for transforming GKS JSON input data from various formats into submission format"""

    submission_container_attribute: SubmissionContainerAttribute
    assertion_criteria: SubmissionAssertionCriteria
    gks_statement_cls: type[GksStatementT]

    @classmethod
    def _read_file(
        cls,
        inputf: typing.TextIO,
    ) -> list[GksStatementT]:
        """Read GKS statements from a JSON file object.

        Expects the JSON file to contain a ``gks_records`` key with a
        list of GKS statement records compatible with the transformer's
        configured ``gks_statement_cls``.

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
    def _get_local_id(variant: CategoricalVariant):
        return variant.id or variant.name

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
    def _get_drugs(therapeutic: TherapyGroup | MappableConcept) -> str:
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

    def _build_shared_submission_kwargs(
        self,
        statement: GksStatementT,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ):
        """Build submission kwargs that are used in both clinical impact and
        oncogenicity submissions

        :param statement: GKS statement
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: _description_
        """
        proposition = statement.proposition
        variant = proposition.subjectVariant

        return {
            "record_status": RecordStatus.NOVEL,
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

    def _build_shared_classification_kwargs(self, description: str | None, therapeutic: TherapyGroup | MappableConcept | None, evidence_lines: list[EvidenceLine] | None, contributions: list[Contribution] | None):
        return {
            "comment": self._get_comment(description, therapeutic),
            "citation": self._get_citations(evidence_lines or []),
            "date_last_evaluated": self._get_date_last_evaluated(
                contributions or []
            ),
        }

    def _get_variant_set(
        self,
        gene_context: MappableConcept,
        variant: CategoricalVariant,
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
        statements: typing.List[GksStatementT],
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
            variant = statement.proposition.subjectVariant

            if not isinstance(variant, CategoricalVariant):
                continue

            variant_hgvs = self._get_variant_hgvs(variant)

            if not variant_hgvs:
                logger.warning(
                    "Skipping statement. No HGVS found for statement ID: %s",
                    statement.id,
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
                    statement.id,
                )
                continue

            submissions.append(
                self._get_submission(
                    statement,
                    observed_in,
                    variant_hgvs=variant_hgvs,
                    submitted_assembly=batch_metadata.submitted_assembly,
                )
            )

        return SubmissionContainer(
            assertion_criteria=self.assertion_criteria,
            **{self.submission_container_attribute: submissions},
        )
