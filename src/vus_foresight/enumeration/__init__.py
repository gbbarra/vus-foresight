"""Variant enumeration (spec section 3).

Tier 1 is mandatory: every coding SNV and every SNV in the flanking window of
every junction. Tier 2 adds intra-codon MNVs, frameshift PTC classes, in-frame
deletions of one and two codons, and exon-level CNVs.

Explicitly out of scope, and not silently omitted:

* arbitrary insertions -- they grow as ``4**k`` per position;
* deep intronic variants outside the declared window;
* complex rearrangements;
* combinations of variants at distinct loci, which depend on phase. Phase is a
  property of a patient, not of a variant, and this system never sees patients.
"""

from __future__ import annotations

from collections.abc import Iterator
from enum import Enum

from ..genome.transcript import FlankSequences, Transcript
from ..variant import Variant
from .cnv import enumerate_exon_cnvs, exon_interval_count
from .frameshift import (
    FrameshiftClass,
    FrameshiftEnumeration,
    MinimalIndel,
    build_frameshift_classes,
    enumerate_frameshift_classes,
    indel_hgvs,
    verify_class_membership,
)
from .inframe import INFRAME_DELETION_LENGTHS, enumerate_inframe_deletions
from .mnv import MNV_PER_CODON, alternative_codons, enumerate_intracodon_mnvs, mnv_count
from .snv import (
    coding_snv_count,
    enumerate_coding_snvs,
    enumerate_intronic_snvs,
    intronic_snv_positions,
)

__all__ = [
    "INFRAME_DELETION_LENGTHS",
    "MNV_PER_CODON",
    "EnumerationClass",
    "FrameshiftClass",
    "FrameshiftEnumeration",
    "MinimalIndel",
    "alternative_codons",
    "build_frameshift_classes",
    "coding_snv_count",
    "enumerate_all",
    "enumerate_coding_snvs",
    "enumerate_exon_cnvs",
    "enumerate_frameshift_classes",
    "enumerate_inframe_deletions",
    "enumerate_intracodon_mnvs",
    "enumerate_intronic_snvs",
    "exon_interval_count",
    "indel_hgvs",
    "intronic_snv_positions",
    "mnv_count",
    "verify_class_membership",
]


class EnumerationClass(str, Enum):
    """Selectable enumeration units, so a run can be scoped from the CLI."""

    CODING_SNV = "coding_snv"
    INTRONIC_SNV = "intronic_snv"
    INTRACODON_MNV = "intracodon_mnv"
    FRAMESHIFT_CLASS = "frameshift_class"
    INFRAME_DELETION = "inframe_deletion"
    EXON_CNV = "exon_cnv"


TIER1: tuple[EnumerationClass, ...] = (
    EnumerationClass.CODING_SNV,
    EnumerationClass.INTRONIC_SNV,
)

TIER2: tuple[EnumerationClass, ...] = (
    EnumerationClass.INTRACODON_MNV,
    EnumerationClass.FRAMESHIFT_CLASS,
    EnumerationClass.INFRAME_DELETION,
    EnumerationClass.EXON_CNV,
)


def enumerate_all(
    transcript: Transcript,
    *,
    classes: tuple[EnumerationClass, ...] = TIER1 + TIER2,
    flanks: FlankSequences | None = None,
    flank_bp: int = 50,
) -> Iterator[Variant]:
    """Yield the selected enumeration classes in a fixed order.

    Order is ``classes`` order, and within each class the enumerator's own
    deterministic order. Nothing here iterates a set.
    """
    for klass in classes:
        if klass is EnumerationClass.CODING_SNV:
            yield from enumerate_coding_snvs(transcript)
        elif klass is EnumerationClass.INTRONIC_SNV:
            yield from enumerate_intronic_snvs(transcript, flanks, flank_bp=flank_bp)
        elif klass is EnumerationClass.INTRACODON_MNV:
            yield from enumerate_intracodon_mnvs(transcript)
        elif klass is EnumerationClass.FRAMESHIFT_CLASS:
            yield from enumerate_frameshift_classes(transcript)
        elif klass is EnumerationClass.INFRAME_DELETION:
            yield from enumerate_inframe_deletions(transcript)
        elif klass is EnumerationClass.EXON_CNV:
            yield from enumerate_exon_cnvs(transcript)
        else:  # pragma: no cover - exhaustive over the enum
            raise ValueError(f"unknown enumeration class {klass}")
