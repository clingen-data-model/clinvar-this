"""Init namespace module"""

from .base import BatchMetadata, GksJsonTransformer, batch_metadata_from_mapping
from .clinical_impact_transformer import ClinicalImpactTransformer
from .oncogenicity_transformer import OncogenicityTransformer
from .resolver import (
    AmbiguousGksStatementType,
    UnsupportedGksStatementType,
    read_gks_json_file,
)

__all__ = [
    "GksJsonTransformer",
    "BatchMetadata",
    "batch_metadata_from_mapping",
    "ClinicalImpactTransformer",
    "OncogenicityTransformer",
    "AmbiguousGksStatementType",
    "UnsupportedGksStatementType",
    "read_gks_json_file",
]
