"""Versioned, offline data source adapters."""

from .base import Adapter, AdapterRegistry, NullAdapter
from .builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from .clinvar import (
    ClinVarRecord,
    ClinVarSnapshot,
    ClinVarSnapshotAdapter,
)
from .clinvar_import import ParseStats, build_snapshot, parse_variant_summary
from .tabular import (
    FrequencyAdapter,
    FunctionalAdapter,
    PredictorAdapter,
    SpliceAdapter,
    TableAdapter,
)

__all__ = [
    "Adapter",
    "AdapterRegistry",
    "ClinVarRecord",
    "ClinVarSnapshot",
    "ClinVarSnapshotAdapter",
    "ParseStats",
    "build_snapshot",
    "FrequencyAdapter",
    "FunctionalAdapter",
    "NullAdapter",
    "PredictorAdapter",
    "RegionAdapter",
    "SpliceAdapter",
    "TableAdapter",
    "TranscriptAdapter",
    "VariantAdapter",
    "parse_variant_summary",
]
