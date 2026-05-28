"""Support for I/O of the ClinGen/CGC/VICC 2022 GKS formatted data to define Oncogenicity submissions.

Example usage:
$ clinvar-this batch import path_to_gks_json -m affected_status=yes -m "collection_method=clinical testing" -m submitted_assembly=GRCh38

"""

from clinvar_this.io.gks_json.base import GksJsonTransformer
from clinvar_api.models import (
    Assembly,
    CitationDb,
    SubmissionAssertionCriteria,
    SubmissionOncogenicitySubmission,
    SomaticOncogenicityClassification,
)
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
from clinvar_api.models.sub_payload import (
    SubmissionObservedInSomatic,
)


class OncogenicityTransformer(GksJsonTransformer[VariantOncogenicityStatement]):
    """Class for transforming ClinGen/CGC/VICC 2022 GKS formatted data to define Oncogenicity submissions"""

    submission_container_attribute = "oncogenicity_submission"
    assertion_criteria = SubmissionAssertionCriteria(
        db=CitationDb.PUBMED,
        id="36063163",  # ClinGen/CGC/VICC Guidelines for Oncogenicity, 2022
    )
    gks_statement_cls = VariantOncogenicityStatement

    def _get_submission(
        self,
        statement: VariantOncogenicityStatement,
        observed_in: list[SubmissionObservedInSomatic],
        variant_hgvs: str | None = None,
        submitted_assembly: Assembly | None = None,
    ) -> SubmissionOncogenicitySubmission:
        """Transform a GKS oncogenicity statement into a ClinVar novel oncogenicity submission.

        These statements support ClinGen/CGC/VICC oncogenicity assertions.

        Local ID will use the proposition's variant ID or name.

        Local Key will use the `record`'s ID.

        :param statement: GKS statement (oncogenicity) to transform
        :param observed_in: List of distinct observations
        :param variant_hgvs: The HGVS expression for a variant, if found
        :param submitted_assembly: The genome assembly used to call the variant.
            Required if `variant_hgvs` is non-null
        :return: The oncogenicity submission corresponding to a GKS Oncogenicity
            statement
        """
        evidence_lines = self._get_evidence_lines(statement.hasEvidenceLines)

        return SubmissionOncogenicitySubmission(
            **self._build_shared_submission_kwargs(
                statement=statement,
                submitted_assembly=submitted_assembly,
                observed_in=observed_in,
                variant_hgvs=variant_hgvs,
            ),
            oncogenicity_classification=SomaticOncogenicityClassification(
                **self._build_shared_classification_kwargs(
                    statement.description,
                    None,
                    evidence_lines,
                    statement.contributions,
                ),
                oncogenicity_classification_description=statement.classification.primaryCoding.code.root.capitalize(),
            ),
        )
