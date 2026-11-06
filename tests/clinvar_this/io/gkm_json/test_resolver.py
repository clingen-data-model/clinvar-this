"""Shared native and bundle reader coverage."""

import json
from pathlib import Path

import pytest
from conftest import BUNDLE_DIR, NATIVE_DIR

from clinvar_api.models import AffectedStatus, CollectionMethod
from clinvar_api.msg.sub_payload import Assembly
from clinvar_this.io.gkm_json import read_gkm_json_file
from clinvar_this.io.gkm_json.base import BatchMetadata
from clinvar_this import exceptions
from clinvar_this.io.gkm_json.resolver import UnsupportedGkmStatementType


@pytest.mark.parametrize(
    "filename",
    ["clinical_impact_civic_assertions.json", "oncogenicity_civic_assertions.json"],
)
def test_native_bundle_equivalence(filename: str) -> None:
    """Both representations produce identical statements and ClinVar payloads."""
    native = NATIVE_DIR / filename
    bundle = BUNDLE_DIR / filename
    schema = BUNDLE_DIR / "bundle.schema.json"

    native_transformer, native_records = read_gkm_json_file(str(native))
    bundle_transformer, bundle_records = read_gkm_json_file(str(bundle), str(schema))

    assert type(native_transformer) is type(bundle_transformer)
    assert native_records == bundle_records

    metadata = BatchMetadata(
        affected_status=AffectedStatus.UNKNOWN,
        collection_method=CollectionMethod.CURATION,
        submitted_assembly=Assembly.NOT_APPLICABLE,
    )
    native_payload = native_transformer.records_to_submission_container(
        native_records, metadata
    )
    bundle_payload = bundle_transformer.records_to_submission_container(
        bundle_records, metadata
    )

    assert native_payload.model_dump(exclude_none=True) == bundle_payload.model_dump(
        exclude_none=True
    )


@pytest.mark.parametrize(
    "filename",
    ["clinical_impact_civic_assertions.json", "oncogenicity_civic_assertions.json"],
)
def test_bundle_statement_collections_are_discovered(filename: str, tmp_path: Path) -> None:
    """Find supported records regardless of collection names."""
    document = json.loads((BUNDLE_DIR / filename).read_text())
    schema = json.loads((BUNDLE_DIR / "bundle.schema.json").read_text())
    records = document.pop("assertion")
    record_schema = schema["properties"].pop("assertion")
    schema["required"].remove("assertion")
    # Split records to catch a loader that stops at the first matching collection.
    for index, (key, record) in enumerate(records.items()):
        name = f"interpretations_{index}"
        document[name] = {key: record}
        schema["properties"][name] = record_schema
    bundle_path = tmp_path / "bundle.json"
    schema_path = tmp_path / "schema.json"
    bundle_path.write_text(json.dumps(document))
    schema_path.write_text(json.dumps(schema))

    _, actual = read_gkm_json_file(str(bundle_path), str(schema_path))
    _, expected = read_gkm_json_file(str(NATIVE_DIR / filename))

    assert actual == expected


def test_bundle_without_supported_statements(tmp_path: Path) -> None:
    """Reject bundles containing only generic evidence statements."""
    document = json.loads((BUNDLE_DIR / "clinical_impact_civic_assertions.json").read_text())
    document["assertion"] = {}
    bundle_path = tmp_path / "bundle.json"
    bundle_path.write_text(json.dumps(document))

    with pytest.raises(exceptions.InvalidFormat, match="nonempty"):
        read_gkm_json_file(str(bundle_path), str(BUNDLE_DIR / "bundle.schema.json"))


def test_bundle_with_mixed_supported_profiles(tmp_path: Path) -> None:
    """Reject mixed profiles even when they occur in separate collections."""
    schema = json.loads((BUNDLE_DIR / "bundle.schema.json").read_text())
    record_schema = schema["properties"]["assertion"]
    document = {}
    for name, filename in (
        ("clinical", "clinical_impact_civic_assertions.json"),
        ("oncogenic", "oncogenicity_civic_assertions.json"),
    ):
        record = json.loads((NATIVE_DIR / filename).read_text())["gkm_records"][0]
        document[name] = {record["id"]: record}
    schema["properties"] = {name: record_schema for name in document}
    schema["required"] = list(document)
    bundle_path = tmp_path / "bundle.json"
    schema_path = tmp_path / "schema.json"
    bundle_path.write_text(json.dumps(document))
    schema_path.write_text(json.dumps(schema))

    with pytest.raises(UnsupportedGkmStatementType):
        read_gkm_json_file(str(bundle_path), str(schema_path))
