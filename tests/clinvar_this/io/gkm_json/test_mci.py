"""Test clinical-impact conversion of an MCI statement."""

from conftest import NATIVE_DIR
from deepdiff import DeepDiff

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
from clinvar_this.io.gkm_json import (
    BatchMetadata,
    ClinicalImpactTransformer,
    read_gkm_json_file,
)


def test_mci_clinical_impact_submission():
    """Convert an Allele with a gene symbol and evidence-line citations."""
    path = NATIVE_DIR / "clinical_impact_mci_statements.json"
    _, records = read_gkm_json_file(str(path))
    metadata = BatchMetadata(
        affected_status=AffectedStatus.UNKNOWN,
        collection_method=CollectionMethod.CURATION,
        submitted_assembly=Assembly.NOT_APPLICABLE,
    )
    expected = SubmissionContainer(
        assertion_criteria=SubmissionAssertionCriteria(
            db=CitationDb.PUBMED, id="27993330"
        ),
        clinical_impact_submission=[
            SubmissionClinicalImpactSubmission(
                record_status=RecordStatus.NOVEL,
                local_id="ga4gh:VA.5ajmZMdX9tr2XjgGNlw-wNl4_yG_FaDe",
                submitted_assembly=Assembly.NOT_APPLICABLE,
                observed_in=[
                    SubmissionObservedInSomatic(
                        affected_status=AffectedStatus.UNKNOWN,
                        allele_origin=AlleleOrigin.SOMATIC,
                        collection_method=CollectionMethod.CURATION,
                    )
                ],
                condition_set=SubmissionConditionSetSomatic(
                    condition=[SubmissionCondition(db="MONDO", id="MONDO:0006412")]
                ),
                variant_set=SubmissionVariantSet(
                    variant=[
                        SubmissionVariant(
                            hgvs="NM_002524.5:c.175G>A",
                            gene=[SubmissionVariantGene(symbol="NRAS")],
                        )
                    ]
                ),
                clinical_impact_classification=SomaticClinicalImpactClassification(
                    clinical_impact_classification_description=SomaticClinicalImpactClassificationDescription.STRONG,
                    assertion_type_for_clinical_impact=SomaticClinicalImpactAssertionType.DIAGNOSTIC_SUPPORTS_DIAGNOSIS,
                    date_last_evaluated="2024-11-19",
                    comment=records[0].description,
                    citation=[
                        SubmissionCitation(db=CitationDb.PUBMED, id="26566875"),
                        SubmissionCitation(db=CitationDb.PUBMED, id="31768065"),
                        SubmissionCitation(db=CitationDb.PUBMED, id="29720485"),
                        SubmissionCitation(db=CitationDb.PUBMED, id="32601132"),
                    ],
                ),
            )
        ],
    )

    actual = ClinicalImpactTransformer().records_to_submission_container(
        records, metadata
    )

    diff = DeepDiff(
        actual.model_dump(exclude_none=True),
        expected.model_dump(exclude_none=True),
        ignore_order=True,
    )
    assert diff == {}
