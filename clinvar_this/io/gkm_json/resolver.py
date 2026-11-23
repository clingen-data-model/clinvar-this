"""Utilities for resolving GKM JSON files.

This module provides helper functions for determining which GKM JSON
transformer should be used for a given input file. Each input file must
contain records corresponding to exactly one supported GKM statement type.

Currently supported statement types:
    * VariantClinicalSignificanceStatement
    * VariantOncogenicityStatement
"""

import io
import json

from ga4gh.gkm.bundles import BundleError, load_bundle

from clinvar_this import exceptions
from clinvar_this.io.gkm_json.clinical_impact_transformer import (
    ClinicalImpactTransformer,
)
from clinvar_this.io.gkm_json.oncogenicity_transformer import (
    OncogenicityTransformer,
)

GkmTransformer = ClinicalImpactTransformer | OncogenicityTransformer


class UnsupportedGkmStatementType(Exception):
    """Raised when a GKM statement type is not supported by clinvar-this."""


class AmbiguousGkmStatementType(Exception):
    """Raised when GKM records match multiple supported statement types."""


def read_gkm_json_file(
    path: str,
    schema_path: str | None = None,
) -> tuple[GkmTransformer, list]:
    """Read GKM native JSON or bundle statements with the matching transformer.

    A schema path selects bundle input. The loader expands references and
    identifies supported statements by model type in every collection.
    All selected records must match one supported statement profile.

    :param path: Path to GKM native JSON or bundle input file
    :param schema_path: Producer schema path for bundle input, if supplied
    :return: Matching transformer and validated statements.
    :raises exceptions.InvalidFormat: If input cannot be read or decoded, records
        are missing or empty, or bundle schema or reference validation fails.
    :raises UnsupportedGkmStatementType: If no profile validates all records.
    :raises AmbiguousGkmStatementType: If multiple profiles validate all records.
    """
    transformers: list[GkmTransformer] = [
        ClinicalImpactTransformer(),
        OncogenicityTransformer(),
    ]
    supported_types = tuple(transformer.gkm_statement_cls for transformer in transformers)

    try:
        if schema_path is not None:
            bundle = load_bundle(path, schema=schema_path)
            records = [
                bundle.export(record, deep=True)
                for name in bundle.collection_names()
                for record in bundle.collection(name).values()
                if isinstance(record, supported_types)
            ]
        else:
            with open(path) as inputf:
                data = json.load(inputf)

            if not isinstance(data, dict) or "gkm_records" not in data:
                raise exceptions.InvalidFormat("Expected GKM native JSON with gkm_records; GKM bundles require --schema PATH")

            records = data["gkm_records"]
        if not isinstance(records, list) or not records:
            raise exceptions.InvalidFormat("GKM input must contain a nonempty list of statements")
    except (BundleError, OSError, json.JSONDecodeError) as error:
        raise exceptions.InvalidFormat(str(error)) from error

    matches = []

    for transformer in transformers:
        try:
            statements = transformer._read_file(io.StringIO(json.dumps({"gkm_records": records})))
        except (
            exceptions.InvalidFormat,
            ValueError,
        ):
            continue

        matches.append((transformer, statements))

    if not matches:
        msg = f"No supported GKM statement type found in GKM JSON file: {path}"
        raise UnsupportedGkmStatementType(msg)

    if len(matches) > 1:
        msg = f"GKM JSON file ambiguously matched multiple supported statement types: {path}"
        raise AmbiguousGkmStatementType(msg)

    return matches[0]
