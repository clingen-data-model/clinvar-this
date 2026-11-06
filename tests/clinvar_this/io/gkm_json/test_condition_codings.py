"""Regression coverage for condition codings without identifiers."""

from pathlib import Path

import pytest
from ga4gh.core.models import Coding, MappableConcept

from clinvar_api.models import CollectionMethod, ConditionDb
from clinvar_this import exceptions
from clinvar_this.io.gkm_json import read_gkm_json_file
from clinvar_this.io.gkm_json.base import BatchMetadata


def test_hp_condition_coding_without_id() -> None:
    """Preserve HP code recognition when the optional coding ID is absent."""
    path = (
        Path(__file__).parents[2]
        / "data/io_gkm_json/native/oncogenicity_civic_assertions.json"
    )
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    statement.proposition.object = MappableConcept(
        primaryCoding=Coding(system="HPO", code="HP:0002664")
    )

    container = transformer.records_to_submission_container(
        [statement], BatchMetadata(collection_method=CollectionMethod.CURATION)
    )

    condition = container.oncogenicity_submission[0].condition_set.condition[0]
    assert condition.db == ConditionDb.HP
    assert condition.id == "HP:0002664"
    assert condition.name is None


@pytest.mark.parametrize("concept_name", [None, "Preferred condition name"])
def test_unknown_condition_system_without_id(concept_name: str | None) -> None:
    """Retain a condition name when its coding system is unsupported."""
    path = (
        Path(__file__).parents[2]
        / "data/io_gkm_json/native/oncogenicity_civic_assertions.json"
    )
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    statement.proposition.object = MappableConcept(
        name=concept_name,
        primaryCoding=Coding(
            system="https://example.org/conditions",
            code="123",
            name="Coding condition name",
        ),
    )

    container = transformer.records_to_submission_container(
        [statement], BatchMetadata(collection_method=CollectionMethod.CURATION)
    )

    condition = container.oncogenicity_submission[0].condition_set.condition[0]
    assert condition.name == (concept_name or "Coding condition name")
    assert condition.db is None
    assert condition.id is None


def test_condition_without_identifier_or_name() -> None:
    """Reject a condition that cannot be represented without guessing its system."""
    path = (
        Path(__file__).parents[2]
        / "data/io_gkm_json/native/oncogenicity_civic_assertions.json"
    )
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    statement.proposition.object = MappableConcept(
        primaryCoding=Coding(system="OMIM", code="155240")
    )

    with pytest.raises(exceptions.InvalidFormat, match="supported identifier"):
        transformer.records_to_submission_container(
            [statement], BatchMetadata(collection_method=CollectionMethod.CURATION)
        )
