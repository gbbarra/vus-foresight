"""Applying a change to the reference and describing what it does."""

from .annotator import (
    CONSEQUENCE_RANK,
    CodingEdit,
    annotate_coding_edits,
    annotate_intronic_substitution,
    genomic_deletion,
    genomic_vcf,
    most_severe,
)
from .indel import (
    ProteinChange,
    StopIndex,
    annotate_inframe_deletion,
    build_stop_index,
    coding_and_utr3,
    normalize_deletion_3prime,
    protein_change,
    ptc_after_indel,
)

__all__ = [
    "CONSEQUENCE_RANK",
    "CodingEdit",
    "ProteinChange",
    "StopIndex",
    "annotate_coding_edits",
    "annotate_inframe_deletion",
    "annotate_intronic_substitution",
    "build_stop_index",
    "coding_and_utr3",
    "genomic_deletion",
    "genomic_vcf",
    "most_severe",
    "normalize_deletion_3prime",
    "protein_change",
    "ptc_after_indel",
]
