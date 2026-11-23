"""Init namespace module"""

from .base import BatchMetadata, GkmJsonTransformer, batch_metadata_from_mapping
from .clinical_impact_transformer import ClinicalImpactTransformer
from .oncogenicity_transformer import OncogenicityTransformer
from .resolver import (
    AmbiguousGkmStatementType,
    UnsupportedGkmStatementType,
    read_gkm_json_file,
)

__all__ = [
    "GkmJsonTransformer",
    "BatchMetadata",
    "batch_metadata_from_mapping",
    "ClinicalImpactTransformer",
    "OncogenicityTransformer",
    "AmbiguousGkmStatementType",
    "UnsupportedGkmStatementType",
    "read_gkm_json_file",
]
