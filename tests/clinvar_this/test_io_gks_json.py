"""Module for testing GKS JSON transformer"""

from deepdiff import DeepDiff
import json
import pathlib
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
import pytest
from ga4gh.va_spec.aac_2017 import (
    VariantTherapeuticResponseStudyStatement,
)

from clinvar_this import exceptions
from clinvar_this.io.gks_json import GksJsonTransformer


DATA_DIR = pathlib.Path(__file__).parent / "data/io_gks_json"


@pytest.fixture(scope="module")
def gks_json_transformer():
    """Create test fixture for GksJsonTransformer"""
    return GksJsonTransformer()


@pytest.fixture(scope="module")
def gks_json_data():
    """Create test fixture for GKS JSON data"""
    with (DATA_DIR / "civic_assertions.json").open() as f:
        return json.load(f)


@pytest.fixture(scope="module")
def civic_aid6(gks_json_data):
    """Create test fixture for CIViC AID6"""
    return VariantTherapeuticResponseStudyStatement(**gks_json_data[0])


@pytest.fixture(scope="module")
def civic_aid7(gks_json_data):
    """Create test fixture for CIViC AID7"""
    return VariantTherapeuticResponseStudyStatement(**gks_json_data[1])


@pytest.fixture(scope="module")
def civic_aid7_submission():
    """Create test fixture for CIViC AID7 submission"""
    return SubmissionClinicalImpactSubmission(
        record_status=RecordStatus.NOVEL,
        local_id="civic.mpid:12",
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
                )
            ]
        ),
        clinical_impact_classification=SomaticClinicalImpactClassification(
            clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
            assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
            comment="Combination treatment of BRAF inhibitor dabrafenib and MEK inhibitor trametinib is recommended for adjuvant treatment of stage III or recurrent melanoma with BRAF V600E mutation detected by the approved THxID kit, as well as first line treatment for metastatic melanoma. The treatments are FDA approved based on studies including the Phase III COMBI-V, COMBI-D and COMBI-AD Trials. Combination therapy is now recommended above BRAF inhibitor monotherapy. Cutaneous squamous-cell carcinoma and keratoacanthoma occur at lower rates with combination therapy than with BRAF inhibitor alone.",
            citation=[
                SubmissionCitation(url="https://identifiers.org/civic.mpid:12"),
                SubmissionCitation(url="https://civicdb.org/links/evidence/3758"),
                SubmissionCitation(id="25399551", db=CitationDb.PUBMED),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6178"),
                SubmissionCitation(id="28891408", db=CitationDb.PUBMED),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6940"),
                SubmissionCitation(id="23020132", db=CitationDb.PUBMED),
                SubmissionCitation(url="https://civicdb.org/links/evidence/6938"),
                SubmissionCitation(id="25265492", db=CitationDb.PUBMED),
            ],
            drug_for_therapeutic_assertion="Trametinib;Dabrafenib",
            date_last_evaluated="2018-05-15",
        ),
    )


@pytest.fixture(scope="module")
def civic_tr_submissions(civic_aid7_submission):
    """Create test fixture for CIViC AID6 and AID7 submissions"""
    return SubmissionContainer(
        assertion_criteria=SubmissionAssertionCriteria(
            db=CitationDb.PUBMED,
            id="27993330",  # AMP/ASCO/CAP
        ),
        clinical_impact_submission=[
            SubmissionClinicalImpactSubmission(
                record_status=RecordStatus.NOVEL,
                local_id="civic.mpid:33",
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
                        )
                    ]
                ),
                clinical_impact_classification=SomaticClinicalImpactClassification(
                    clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
                    assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.THERAPEUTIC_SENSITIVITY_RESPONSE,
                    comment="L858R is among the most common sensitizing EGFR mutations in NSCLC, and is assessed via DNA mutational analysis, including Sanger sequencing and next generation sequencing methods. Tyrosine kinase inhibitor afatinib is FDA approved as a first line systemic therapy in NSCLC with sensitizing EGFR mutation (civic.EID:2997).",
                    citation=[
                        SubmissionCitation(url="https://identifiers.org/civic.mpid:33"),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/2997"
                        ),
                        SubmissionCitation(id="23982599", db=CitationDb.PUBMED),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/2629"
                        ),
                        SubmissionCitation(id="18408761", db=CitationDb.PUBMED),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/982"
                        ),
                        SubmissionCitation(id="24439929", db=CitationDb.PUBMED),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/968"
                        ),
                        SubmissionCitation(id="26515464", db=CitationDb.PUBMED),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/883"
                        ),
                        SubmissionCitation(id="22452895", db=CitationDb.PUBMED),
                        SubmissionCitation(
                            url="https://civicdb.org/links/evidence/879"
                        ),
                        SubmissionCitation(id="23816960", db=CitationDb.PUBMED),
                    ],
                    drug_for_therapeutic_assertion="Afatinib",
                    date_last_evaluated="2018-02-23",
                ),
            ),
            civic_aid7_submission,
        ],
    )


def test_read_file(gks_json_transformer, civic_aid6, civic_aid7):
    """Ensure that read_file method works correctly"""
    path = DATA_DIR / "civic_assertions.json"
    expected_assertions = [civic_aid6, civic_aid7]

    with path.open("rt") as inputf:
        actual = gks_json_transformer.read_file(file=inputf)
    assert actual == expected_assertions

    actual = gks_json_transformer.read_file(path=path)
    assert actual == expected_assertions

    with pytest.raises(exceptions.InvalidFormat, match=r"Error decoding JSON"):
        gks_json_transformer.read_file(path=DATA_DIR / "example_bad.json")


def test_records_to_submission_container(
    gks_json_transformer, civic_aid6, civic_aid7, civic_tr_submissions
):
    """Ensure that records_to_submission_container works correctly"""
    # Test single therapy and CombinationTherapy
    actual = gks_json_transformer.records_to_submission_container(
        [civic_aid6, civic_aid7]
    )
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_tr_submissions.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}

    # Test TherapeuticSubstituteGroup
    civic_aid7_cpy = civic_aid7.model_copy()
    civic_aid7_cpy.proposition.objectTherapeutic.root.groupType.name = (
        "TherapeuticSubstituteGroup"
    )
    actual = gks_json_transformer.records_to_submission_container([civic_aid7_cpy])
    civic_tr_submissions_cpy = civic_tr_submissions.model_copy()
    civic_tr_submissions_cpy.clinical_impact_submission.pop(0)
    civic_tr_submissions_cpy = civic_tr_submissions_cpy.model_dump(exclude_none=True)
    civic_tr_submissions_cpy["clinical_impact_submission"][0][
        "clinical_impact_classification"
    ][
        "comment"
    ] = f"{civic_tr_submissions_cpy['clinical_impact_submission'][0]['clinical_impact_classification']['comment']} NOTE: These therapies are in substitution."
    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        civic_tr_submissions_cpy,
        ignore_order=True,
    )
    assert diff == {}
