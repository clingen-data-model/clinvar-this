"""Module for testing GKS Oncogenicity Transformer"""

import json

from conftest import DATA_DIR
from deepdiff import DeepDiff
from ga4gh.va_spec.ccv_2022 import VariantOncogenicityStatement
import pytest

from clinvar_api.models import (
    AffectedStatus,
    AlleleOrigin,
    CitationDb,
    CollectionMethod,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionCitation,
    SubmissionCondition,
    SubmissionContainer,
    SubmissionOncogenicitySubmission,
    SubmissionVariant,
    SubmissionVariantSet,
)
from clinvar_api.models.sub_payload import (
    OncogenicityClassificationDescription,
    SomaticOncogenicityClassification,
    SubmissionConditionSetSomatic,
    SubmissionObservedInSomatic,
    SubmissionVariantGene,
)
from clinvar_api.msg.sub_payload import Assembly
from clinvar_this.io.gks_json.base import BatchMetadata
from clinvar_this.io.gks_json.oncogenicity_transformer import OncogenicityTransformer


@pytest.fixture(scope="module")
def oncogenicity_transformer():
    """Create test fixture for ClinicalImpactTransformer"""
    return OncogenicityTransformer()


@pytest.fixture(scope="module")
def civic_metadata():
    """Create test fixture for CIViC metadata to use"""
    return BatchMetadata(
        affected_status=AffectedStatus.UNKNOWN,
        collection_method=CollectionMethod.CURATION,
        submitted_assembly=Assembly.GRCH37,
    )


@pytest.fixture(scope="module")
def oncogenicity_gks_json_data():
    """Create test fixture for Oncogenicity GKS JSON data"""
    with (DATA_DIR / "oncogenicity_civic_assertions.json").open() as f:
        return json.load(f)["gks_records"]


@pytest.fixture(scope="module")
def civic_aid202(oncogenicity_gks_json_data):
    """Create test fixture for CIViC AID 202"""
    return VariantOncogenicityStatement(**oncogenicity_gks_json_data[0])


@pytest.fixture(scope="module")
def oncogenicity_assertion_criteria():
    """Create test fixture for AMP/ASCO/CAP 2017 assertion criteria"""
    return SubmissionAssertionCriteria(
        db=CitationDb.PUBMED,
        id="36063163",
    )


@pytest.fixture(scope="module")
def civic_aid202_submission():
    """Create test fixture for CIViC AID 202 submission"""
    return SubmissionOncogenicitySubmission(
        record_status=RecordStatus.NOVEL,
        local_id="civic.mpid:113",
        submitted_assembly=Assembly.GRCH37,
        local_key="civic.aid:202",
        observed_in=[
            SubmissionObservedInSomatic(
                affected_status=AffectedStatus.UNKNOWN,
                allele_origin=AlleleOrigin.SOMATIC,
                collection_method=CollectionMethod.CURATION,
            )
        ],
        condition_set=SubmissionConditionSetSomatic(
            condition=[
                SubmissionCondition(
                    name="Medullary Thyroid Carcinoma",
                )
            ]
        ),
        variant_set=SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs="NM_020975.4:c.2753T>C",
                    gene=[SubmissionVariantGene(symbol="RET")],
                    alternate_designations=["MET918THR"],
                )
            ]
        ),
        oncogenicity_classification=SomaticOncogenicityClassification(
            oncogenicity_classification_description=OncogenicityClassificationDescription.LIKELY_ONCOGENIC,
            comment="Published sequencing studies have shown that RET mutations are very common in medullary thryoid carcinoma (MTC) and M918T is the most common specific variant, especially in the MEN2B clinical subtype of familial disease (civic.EID:78) but also in sporadic cases(civic.EID:12800). M918T mutations may predict worse outcomes (civic.EID:74). Biochemical and functional characterization demonstrates that the M918T mutation leads to functional activation of RET relative to wild-type through multiple complementary mechanisms, including increased ATP affinity (>10-fold) and complex stability, reduced conformational rigidity, and the promotion of ligand-independent dimerization and autophosphorylation (civic.EID:12805). Exogenous expression has been shown to induce transformation of Ba/F3 cells (civic.EID:11723), and drive colony formation in NIH3T3 cells (civic.EID:12709, OS2). RET M918T occurs in the region of the tyrosine kinase domain which is associated with multiple endocrine neoplasia type 2 B (OM1). RET M918T is predicted to be deleterious (CHASMplus score 0.314 > VECS gene-specific cutoff of 0.22, OP1). Eleven instances of the variant occur in cancerhotspots.org (V2): 6 Thyroid, 4 Adrenal Gland, 1 Breast (OP3). The variant is absent in gnomAD database (v4.1.0, OP4). Together these criteria indicate that M918T is likely oncogenic, with a score of 9.",
            citation=[
                SubmissionCitation(url="https://civicdb.org/links/evidence/74"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/12800"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/78"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/12711"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/12805"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/11723"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/12709"),
            ],
            date_last_evaluated="2026-04-16",
        ),
    )


@pytest.fixture(scope="module")
def civic_oncogenicity_submissions(
    civic_aid202_submission, oncogenicity_assertion_criteria
):
    """Create test fixture for CIViC oncogenicity submissions"""
    return SubmissionContainer(
        assertion_criteria=oncogenicity_assertion_criteria,
        oncogenicity_submission=[civic_aid202_submission],
    )


def test_read_file(
    oncogenicity_transformer,
    civic_aid202,
):
    """Ensure that read_file method works correctly"""
    path = DATA_DIR / "oncogenicity_civic_assertions.json"
    expected_assertions = [
        civic_aid202,
    ]

    with path.open("rt") as inputf:
        actual = oncogenicity_transformer.read_file(file=inputf)
    assert actual == expected_assertions

    actual = oncogenicity_transformer.read_file(path=path)
    assert actual == expected_assertions


def test_records_to_submission_container(
    oncogenicity_transformer,
    civic_metadata,
    civic_aid202,
    civic_oncogenicity_submissions,
):
    """Ensure that records_to_submission_container works correctly for oncogenicity statements"""
    actual = oncogenicity_transformer.records_to_submission_container(
        [civic_aid202], civic_metadata
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_oncogenicity_submissions.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}
