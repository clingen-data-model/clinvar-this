"""Module for testing GKS Clinical Impact Transformer"""

import json
import re

from conftest import DATA_DIR
from deepdiff import DeepDiff
from ga4gh.va_spec.aac_2017 import VariantClinicalSignificanceStatement
import pytest

from clinvar_api.models import (
    AffectedStatus,
    AlleleOrigin,
    CitationDb,
    CollectionMethod,
    RecordStatus,
    SubmissionAssertionCriteria,
    SubmissionCitation,
    SubmissionClinicalImpactSubmission,
    SubmissionCondition,
    SubmissionContainer,
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
    Assembly,
    SomaticClinicalImpactAssertionType,
    SomaticClinicalImpactClassificationDescription,
)
from clinvar_this import exceptions
from clinvar_this.io.gks_json.base import BatchMetadata
from clinvar_this.io.gks_json.clinical_impact_transformer import (
    ClinicalImpactTransformer,
)


@pytest.fixture(scope="module")
def clinical_impact_transformer():
    """Create test fixture for ClinicalImpactTransformer"""
    return ClinicalImpactTransformer()


@pytest.fixture(scope="module")
def civic_metadata():
    """Create test fixture for CIViC metadata to use"""
    return BatchMetadata(
        affected_status=AffectedStatus.UNKNOWN,
        collection_method=CollectionMethod.CURATION,
        submitted_assembly=Assembly.GRCH37,
    )


@pytest.fixture(scope="module")
def clinical_impact_gks_json_data():
    """Create test fixture for AAC 2017 GKS JSON data"""
    with (DATA_DIR / "clinical_impact_civic_assertions.json").open() as f:
        return json.load(f)["gks_records"]


@pytest.fixture(scope="module")
def civic_aid6(clinical_impact_gks_json_data):
    """Create test fixture for CIViC AID6"""
    return VariantClinicalSignificanceStatement(**clinical_impact_gks_json_data[0])


@pytest.fixture(scope="module")
def civic_aid7(clinical_impact_gks_json_data):
    """Create test fixture for CIViC AID7"""
    return VariantClinicalSignificanceStatement(**clinical_impact_gks_json_data[1])


@pytest.fixture(scope="module")
def civic_aid20(clinical_impact_gks_json_data):
    """Create test fixture for CIViC AID20"""
    return VariantClinicalSignificanceStatement(**clinical_impact_gks_json_data[2])


@pytest.fixture(scope="module")
def civic_aid9(clinical_impact_gks_json_data):
    """Create test fixture for CIViC AID9"""
    return VariantClinicalSignificanceStatement(**clinical_impact_gks_json_data[3])


@pytest.fixture(scope="module")
def civic_aid200(clinical_impact_gks_json_data):
    """Create test fixture for CIViC AID200"""
    return VariantClinicalSignificanceStatement(**clinical_impact_gks_json_data[4])


@pytest.fixture(scope="module")
def amp_asco_cap_assertion_criteria():
    """Create test fixture for AMP/ASCO/CAP 2017 assertion criteria"""
    return SubmissionAssertionCriteria(
        db=CitationDb.PUBMED,
        id="27993330",  # AMP/ASCO/CAP
    )


@pytest.fixture(scope="module")
def civic_aid7_submission():
    """Create test fixture for CIViC AID7 submission"""
    return SubmissionClinicalImpactSubmission(
        record_status=RecordStatus.NOVEL,
        local_id="civic.mpid:12",
        submitted_assembly=Assembly.GRCH37,
        local_key="civic.aid:7",
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
                    name="Melanoma",
                )
            ]
        ),
        variant_set=SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs="NM_004333.4:c.1799T>A",
                    gene=[SubmissionVariantGene(symbol="BRAF")],
                    alternate_designations=["VAL600GLU", "V640E", "VAL640GLU"],
                )
            ]
        ),
        clinical_impact_classification=SomaticClinicalImpactClassification(
            clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
            assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
            comment="Combination treatment of BRAF inhibitor dabrafenib and MEK inhibitor trametinib is recommended for adjuvant treatment of stage III or recurrent melanoma with BRAF V600E mutation detected by the approved THxID kit, as well as first line treatment for metastatic melanoma. The treatments are FDA approved based on studies including the Phase III COMBI-V, COMBI-D and COMBI-AD Trials. Combination therapy is now recommended above BRAF inhibitor monotherapy. Cutaneous squamous-cell carcinoma and keratoacanthoma occur at lower rates with combination therapy than with BRAF inhibitor alone.",
            citation=[
                SubmissionCitation(url="https://civicdb.org/links/evidence/3758"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/25399551"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6178"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/28891408"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6940"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/23020132"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6938"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/25265492"),
            ],
            drug_for_therapeutic_assertion="Dabrafenib;Trametinib",
        ),
    )


@pytest.fixture(scope="module")
def civic_tr_submissions(civic_aid7_submission, amp_asco_cap_assertion_criteria):
    """Create test fixture for CIViC AID6 and AID7 submissions"""
    return SubmissionContainer(
        assertion_criteria=amp_asco_cap_assertion_criteria,
        clinical_impact_submission=[
            SubmissionClinicalImpactSubmission(
                record_status=RecordStatus.NOVEL,
                local_id="civic.mpid:33",
                submitted_assembly=Assembly.GRCH37,
                local_key="civic.aid:6",
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
                            name="Lung Non-small Cell Carcinoma",
                        )
                    ]
                ),
                variant_set=SubmissionVariantSet(
                    variant=[
                        SubmissionVariant(
                            hgvs="NM_005228.4:c.2573T>G",
                            gene=[SubmissionVariantGene(symbol="EGFR")],
                            alternate_designations=["LEU858ARG", "L813R", "LEU813ARG"],
                        )
                    ]
                ),
                clinical_impact_classification=SomaticClinicalImpactClassification(
                    clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
                    assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
                    comment="L858R is among the most common sensitizing EGFR mutations in NSCLC, and is assessed via DNA mutational analysis, including Sanger sequencing and next generation sequencing methods. Tyrosine kinase inhibitor afatinib is FDA approved as a first line systemic therapy in NSCLC with sensitizing EGFR mutation (civic.EID:2997).",
                    citation=[
                        # SubmissionCitation(url="https://identifiers.org/civic.mpid:33"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/2997"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/23982599"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/2629"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/18408761"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/982"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/24439929"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/968"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/26515464"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/883"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/22452895"),
                        SubmissionCitation(url="https://civicdb.org/links/evidence/879"),
                        SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/23816960"),
                    ],
                    drug_for_therapeutic_assertion="Afatinib",
                ),
            ),
            civic_aid7_submission,
        ],
    )


@pytest.fixture(scope="module")
def civic_aid9_submission():
    """Create test fixture for CIViC AID9 submission"""
    return SubmissionClinicalImpactSubmission(
        record_status=RecordStatus.NOVEL,
        local_id="civic.mpid:1594",
        submitted_assembly=Assembly.GRCH37,
        local_key="civic.aid:9",
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
                    name="Diffuse Midline Glioma, H3 K27M-mutant",
                )
            ]
        ),
        variant_set=SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs="NM_001105.4:c.983G>T",
                    gene=[SubmissionVariantGene(symbol="ACVR1")],
                    alternate_designations=["GLY328VAL"],
                )
            ]
        ),
        clinical_impact_classification=SomaticClinicalImpactClassification(
            clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.POTENTIAL,
            assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.DIAGNOSTIC_SUPPORTS_DIAGNOSIS,
            comment="ACVR1 G328V mutations occur within the kinase domain, leading to activation of downstream signaling. Exclusively seen in high-grade pediatric gliomas, supporting diagnosis of diffuse intrinsic pontine glioma.",
            citation=[
                # SubmissionCitation(url="https://identifiers.org/civic.mpid:1594"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/4846"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/24705250"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6955"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/24705254"),
            ],
        ),
    )


@pytest.fixture(scope="module")
def civic_diagnostic_submissions(civic_aid9_submission, amp_asco_cap_assertion_criteria):
    """Create test fixture for CIViC AID9 submission"""
    return SubmissionContainer(
        assertion_criteria=amp_asco_cap_assertion_criteria,
        clinical_impact_submission=[civic_aid9_submission],
    )


@pytest.fixture(scope="module")
def civic_aid20_submission():
    """Create test fixture for CIViC AID20 submission"""
    return SubmissionClinicalImpactSubmission(
        record_status=RecordStatus.NOVEL,
        local_id="civic.mpid:12",
        submitted_assembly=Assembly.GRCH37,
        local_key="civic.aid:20",
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
                    name="Colorectal Cancer",
                )
            ]
        ),
        variant_set=SubmissionVariantSet(
            variant=[
                SubmissionVariant(
                    hgvs="NM_004333.4:c.1799T>A",
                    gene=[SubmissionVariantGene(symbol="BRAF")],
                    alternate_designations=["VAL600GLU", "V640E", "VAL640GLU"],
                )
            ]
        ),
        clinical_impact_classification=SomaticClinicalImpactClassification(
            clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
            assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.PROGNOSTIC_POOR_OUTCOME,
            comment="BRAF V600E was associated with worse prognosis in Phase II and III colorectal cancer, with a stronger effect in MSI-Low or MSI-Stable tumors. In metastatic CRC, V600E was associated with worse prognosis, and meta-analysis showed BRAF mutation in CRC associated with multiple negative prognostic markers.",
            citation=[
                # SubmissionCitation(url="https://identifiers.org/civic.mpid:12"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/7159"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/24112392"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/7158"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/21641636"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/7157"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/21502544"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/7156"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/20008640"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/103"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/24594804"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/1552"),
                SubmissionCitation(url="https://pubmed.ncbi.nlm.nih.gov/27404270"),
            ],
        ),
    )


@pytest.fixture(scope="module")
def civic_prognostic_submissions(civic_aid20_submission, amp_asco_cap_assertion_criteria):
    """Create test fixture for CIViC AID20 submission"""
    return SubmissionContainer(
        assertion_criteria=amp_asco_cap_assertion_criteria,
        clinical_impact_submission=[civic_aid20_submission],
    )


def test_read_file(
    clinical_impact_transformer,
    civic_aid6,
    civic_aid7,
    civic_aid9,
    civic_aid20,
    civic_aid200,
):
    """Ensure that read_file method works correctly"""
    path = DATA_DIR / "clinical_impact_civic_assertions.json"
    expected_assertions = [
        civic_aid6,
        civic_aid7,
        civic_aid20,
        civic_aid9,
        civic_aid200,
    ]

    with path.open("rt") as inputf:
        actual = clinical_impact_transformer.read_file(file=inputf)
    assert actual == expected_assertions

    actual = clinical_impact_transformer.read_file(path=path)
    assert actual == expected_assertions

    with pytest.raises(exceptions.InvalidFormat, match="Error decoding GKS JSON"):
        clinical_impact_transformer.read_file(path=DATA_DIR / "invalid_json.json")

    with pytest.raises(
        exceptions.InvalidFormat,
        match=re.escape(
            "Invalid GKS JSON: missing required key `gks_records` (must be a list of statements)"
        ),
    ):
        clinical_impact_transformer.read_file(path=DATA_DIR / "no_gks_records.json")


def test_records_to_submission_container(
    clinical_impact_transformer,
    civic_metadata,
    civic_aid6,
    civic_aid7,
    civic_tr_submissions,
):
    """Ensure that records_to_submission_container works correctly for therapeutic statements"""
    # Test single therapy and CombinationTherapy
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid6, civic_aid7], civic_metadata
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_tr_submissions.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}

    # Test TherapeuticSubstituteGroup
    civic_aid7_cpy = civic_aid7.model_copy()
    civic_aid7_cpy.hasEvidenceLines[
        0
    ].targetProposition.objectTherapeutic.root.membershipOperator = "OR"
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid7_cpy], civic_metadata
    )
    civic_tr_submissions_cpy = civic_tr_submissions.model_copy()
    civic_tr_submissions_cpy.clinical_impact_submission.pop(0)
    civic_tr_submissions_cpy = civic_tr_submissions_cpy.model_dump(exclude_none=True)
    civic_tr_submissions_cpy["clinical_impact_submission"][0]["clinical_impact_classification"][
        "comment"
    ] = f"{civic_tr_submissions_cpy['clinical_impact_submission'][0]['clinical_impact_classification']['comment']} NOTE: These therapies are in substitution."
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_tr_submissions_cpy,
        ignore_order=True,
    )
    assert diff == {}


def test_records_to_submission_container_diagnostic(
    clinical_impact_transformer,
    civic_metadata,
    civic_aid9,
    civic_diagnostic_submissions,
):
    """Ensure that records_to_submission_container works correctly for diagnostic statements"""
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid9], civic_metadata
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_diagnostic_submissions.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}


def test_records_to_submission_container_prognostic(
    clinical_impact_transformer,
    civic_metadata,
    civic_aid20,
    civic_prognostic_submissions,
):
    """Ensure that records_to_submission_container works correctly for prognostic statements"""
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid20], civic_metadata
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_prognostic_submissions.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}


def test_citations(
    clinical_impact_transformer,
    civic_metadata,
    civic_aid20,
    amp_asco_cap_assertion_criteria,
    civic_aid20_submission,
):
    """Test that citations work correctly when evidence line does not have `hasEvidenceItems` and only has extension for citations"""
    # only care about civic urls in this use case
    civic_aid20_submission_cpy = civic_aid20_submission.model_dump()
    new_citations = []
    for citation in civic_aid20_submission_cpy["clinical_impact_classification"]["citation"]:
        if citation["url"].startswith("https://civicdb.org"):
            new_citations.append(citation)
    civic_aid20_submission_cpy["clinical_impact_classification"]["citation"] = new_citations

    expected = SubmissionContainer(
        assertion_criteria=amp_asco_cap_assertion_criteria,
        clinical_impact_submission=[civic_aid20_submission_cpy],
    )

    # Ensure hasEvidenceItems is none
    civic_aid20_cpy = civic_aid20.model_copy(deep=True)
    civic_aid20_cpy.hasEvidenceLines[0].hasEvidenceItems = None
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid20_cpy], civic_metadata
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        expected.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}


def test_contributions(clinical_impact_transformer, civic_aid200, civic_metadata):
    """Test that contributions work correctly"""
    actual = clinical_impact_transformer.records_to_submission_container(
        [civic_aid200], civic_metadata
    )
    assert len(actual.clinical_impact_submission) == 1

    assert (
        actual.clinical_impact_submission[0].clinical_impact_classification.date_last_evaluated
        == "2026-04-16"
    )


def test_no_evidence_lines(
    clinical_impact_transformer,
    civic_aid20_submission,
    amp_asco_cap_assertion_criteria,
    civic_aid20,
    civic_metadata,
):
    """Test that statement with no evidence lines works correctly

    Do not expect any citation or assertion_type_for_clinical_impact
    """
    assertion_copy = civic_aid20.model_copy(deep=True)
    assertion_copy.hasEvidenceLines = None

    actual = clinical_impact_transformer.records_to_submission_container(
        [assertion_copy], civic_metadata
    )

    submission_copy = civic_aid20_submission.model_dump()
    submission_copy["clinical_impact_classification"]["citation"] = []
    submission_copy["clinical_impact_classification"]["assertion_type_for_clinical_impact"] = None
    expected = SubmissionContainer(
        assertion_criteria=amp_asco_cap_assertion_criteria,
        clinical_impact_submission=[SubmissionClinicalImpactSubmission(**submission_copy)],
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        expected.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}
