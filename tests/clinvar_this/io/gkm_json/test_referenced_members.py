"""Regression coverage for referenced categorical-variant members."""

from pathlib import Path

import pytest
from ga4gh.cat_vrs.models import CategoricalVariant
from ga4gh.core.models import iriReference

from clinvar_api.models import CollectionMethod
from clinvar_this.io.gkm_json import read_gkm_json_file
from clinvar_this.io.gkm_json.base import BatchMetadata


@pytest.mark.parametrize(
    "filename",
    ["clinical_impact_civic_assertions.json", "oncogenicity_civic_assertions.json"],
)
def test_referenced_member_preserves_submission(filename: str) -> None:
    """An unresolved member does not prevent conversion of an embedded variant."""
    path = Path(__file__).parents[2] / "data/io_gkm_json/native" / filename
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    metadata = BatchMetadata(collection_method=CollectionMethod.CURATION)
    expected = transformer.records_to_submission_container([statement], metadata)
    submissions = getattr(expected, transformer.submission_container_attribute)
    assert len(submissions) == 1
    assert submissions[0].variant_set.variant[0].hgvs

    variant = statement.proposition.subject
    assert isinstance(variant, CategoricalVariant)
    variant.members = [
        iriReference("https://example.org/variants/unresolved"),
        *(variant.members or []),
    ]

    actual = transformer.records_to_submission_container([statement], metadata)

    assert actual.model_dump(exclude_none=True) == expected.model_dump(exclude_none=True)
