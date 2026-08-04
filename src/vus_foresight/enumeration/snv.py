"""Tier 1: every possible single-nucleotide variant.

"Every possible" is meant literally -- this is not a set of observed variants.
The cardinality is therefore a closed-form invariant (``3 x CDS length`` for the
coding set), which is what makes the layer-1 tests exact rather than
approximate.
"""

from __future__ import annotations

from typing import Iterator

from ..annotate.annotator import annotate_coding_edits, annotate_intronic_substitution
from ..annotate.annotator import CodingEdit
from ..genome.sequence import alternatives
from ..genome.transcript import FlankSequences, Transcript
from ..variant import Variant

__all__ = [
    "enumerate_coding_snvs",
    "enumerate_intronic_snvs",
    "coding_snv_count",
    "intronic_snv_positions",
]


def coding_snv_count(transcript: Transcript) -> int:
    """Closed form: three alternatives at every CDS base."""
    return 3 * transcript.cds_length


def enumerate_coding_snvs(transcript: Transcript) -> Iterator[Variant]:
    """Yield every coding SNV, in CDS position then canonical base order.

    Ordering is part of the contract: the output must be byte-reproducible
    across runs (spec section 11), and that starts with the enumerator being a
    deterministic sequence rather than a set.
    """
    for cds_position, ref in transcript.iter_cds_positions():
        for alt in alternatives(ref):
            yield annotate_coding_edits(
                transcript,
                (CodingEdit(cds_position=cds_position, ref=ref, alt=alt),),
            )


def intronic_snv_positions(transcript: Transcript, *, flank_bp: int = 50) -> list[tuple[int, int]]:
    """``(intron_index, offset_from_donor)`` for every enumerated intronic base.

    Both sides of every junction are covered: offsets ``1..flank_bp`` from the
    donor and the mirror positions counted back from the acceptor. Introns
    shorter than ``2 x flank_bp`` would otherwise have their two windows
    overlap, so positions are deduplicated -- an intron of length 30 contributes
    30 positions, not 60.
    """
    seen: set[tuple[int, int]] = set()
    ordered: list[tuple[int, int]] = []
    for index, intron in enumerate(transcript.introns):
        donor_side = range(1, min(flank_bp, intron.length) + 1)
        acceptor_side = range(max(1, intron.length - flank_bp + 1), intron.length + 1)
        for n in list(donor_side) + list(acceptor_side):
            key = (index, n)
            if key in seen:
                continue
            seen.add(key)
            ordered.append(key)
    ordered.sort()
    return ordered


def enumerate_intronic_snvs(
    transcript: Transcript,
    flanks: FlankSequences | None = None,
    *,
    flank_bp: int = 50,
) -> Iterator[Variant]:
    """Yield every SNV in the flanking window of every junction.

    Without flanking sequence the reference base is unknown; all four bases are
    then emitted as alternatives minus none, and each row is marked
    ``reference_base_unknown``. This keeps the *positional* coverage complete --
    which is what the gap map needs, since splice criteria are driven by
    position and by SpliceAI, not by the reference base -- while making the
    missing information explicit instead of guessed.
    """
    flanks = flanks or FlankSequences()
    introns = transcript.introns
    for index, n in intronic_snv_positions(transcript, flank_bp=flank_bp):
        intron = introns[index]
        ref = flanks.base_at(transcript, intron, n)
        alts = alternatives(ref) if ref != "N" else ("A", "C", "G", "T")
        for alt in alts:
            yield annotate_intronic_substitution(transcript, intron, n, ref, alt)
