"""Utilities for resolving GKS JSON files.

This module provides helper functions for determining which GKS JSON
transformer should be used for a given input file. Each input file must
contain records corresponding to exactly one supported GKS statement type.

Currently supported statement types:
    * VariantClinicalSignificanceStatement
    * VariantOncogenicityStatement
"""

import io
import json

from ga4gh.gkm.bundles import BundleError, load_bundle

from clinvar_this import exceptions
from clinvar_this.io.gks_json.clinical_impact_transformer import (
    ClinicalImpactTransformer,
)
from clinvar_this.io.gks_json.oncogenicity_transformer import (
    OncogenicityTransformer,
)

GksTransformer = ClinicalImpactTransformer | OncogenicityTransformer


class UnsupportedGksStatementType(Exception):
    """Raised when a GKS statement type is not supported by clinvar-this."""


class AmbiguousGksStatementType(Exception):
    """Raised when GKS records match multiple supported statement types."""


def read_gks_json_file(
    path: str,
    schema_path: str | None = None,
) -> tuple[GksTransformer, list]:
    """Read GKS native JSON or bundle statements with the matching transformer.

    A schema path selects bundle input. The loader expands references and
    identifies supported statements by model type in every collection.
    All selected records must match one supported statement profile.

    :param path: Path to GKS native JSON or bundle input file
    :param schema_path: Producer schema path for bundle input, if supplied
    :return: Matching transformer and validated statements.
    :raises exceptions.InvalidFormat: If input cannot be read or decoded, records
        are missing or empty, or bundle schema or reference validation fails.
    :raises UnsupportedGksStatementType: If no profile validates all records.
    :raises AmbiguousGksStatementType: If multiple profiles validate all records.
    """
    transformers: list[GksTransformer] = [
        ClinicalImpactTransformer(),
        OncogenicityTransformer(),
    ]
    supported_types = tuple(transformer.gks_statement_cls for transformer in transformers)

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

            if not isinstance(data, dict) or "gks_records" not in data:
                raise exceptions.InvalidFormat("Expected GKS native JSON with gks_records; GKS bundles require --schema PATH")

            records = data["gks_records"]
        if not isinstance(records, list) or not records:
            raise exceptions.InvalidFormat("GKS input must contain a nonempty list of statements")
    except (BundleError, OSError, json.JSONDecodeError) as error:
        raise exceptions.InvalidFormat(str(error)) from error

    matches = []

    for transformer in transformers:
        try:
            statements = transformer._read_file(io.StringIO(json.dumps({"gks_records": records})))
        except (
            exceptions.InvalidFormat,
            ValueError,
        ):
            continue

        matches.append((transformer, statements))

    if not matches:
        msg = f"No supported GKS statement type found in GKS JSON file: {path}"
        raise UnsupportedGksStatementType(msg)

    if len(matches) > 1:
        msg = f"GKS JSON file ambiguously matched multiple supported statement types: {path}"
        raise AmbiguousGksStatementType(msg)

    return matches[0]
