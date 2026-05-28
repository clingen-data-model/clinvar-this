"""Utilities for resolving GKS JSON files.

This module provides helper functions for determining which GKS JSON
transformer should be used for a given input file. Each input file must
contain records corresponding to exactly one supported GKS statement type.

Currently supported statement types:
    * VariantClinicalSignificanceStatement
    * VariantOncogenicityStatement
"""

from clinvar_this import exceptions
from clinvar_this.io.gks_json.clinical_impact_transformer import (
    ClinicalImpactTransformer,
)
from clinvar_this.io.gks_json.oncogenicity_transformer import OncogenicityTransformer

GksTransformer = ClinicalImpactTransformer | OncogenicityTransformer


class UnsupportedGksStatementType(Exception):
    """Raised when a GKS statement type is not supported by clinvar-this."""


class AmbiguousGksStatementType(Exception):
    """Raised when GKS records match multiple supported statement types."""


def read_gks_json_file(
    path: str,
) -> tuple[GksTransformer, list]:
    """Read a GKS JSON file using the matching transformer.

    Attempts to validate the input file using each supported GKS transformer's statement
    class. All records in the file MUST correspond to only one GKS statement class.

    :param path: Path to GKS JSON input file
    :raise UnsupportedGksStatementType: If no supported GKS statement type is detected
    :raise AmbiguousGksStatementType: If multiple records validate against more than
        one GKS statement type
    :return:
        Tuple containing:
            * Matching GKS transformer instance
            * List of GKS statements
    """
    transformers: list[GksTransformer] = [
        ClinicalImpactTransformer(),
        OncogenicityTransformer(),
    ]

    matches = []

    for transformer in transformers:
        try:
            with open(path) as inputf:
                statements = transformer._read_file(inputf)
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
