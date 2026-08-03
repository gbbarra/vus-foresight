"""Tier 2: intra-codon multi-nucleotide variants.

Each codon has 63 alternatives; 9 differ at a single position and are already
covered by Tier 1, leaving exactly **54** multi-nucleotide alternatives per
codon. They are enumerated because their protein consequence is not a function
of their component substitutions -- ``GCC``(Ala) can become ``ATC``(Ile) while
``ACC``(Thr) and ``GTC``(Val) are what the two single substitutions give -- and
no standard annotator composes that correctly.
"""

from __future__ import annotations

from itertools import product
from typing import Iterator

from ..annotate.annotator import CodingEdit, annotate_coding_edits
from ..genome.sequence import BASES
from ..genome.transcript import Transcript
from ..variant import Variant, VariantKind

__all__ = ["MNV_PER_CODON", "enumerate_intracodon_mnvs", "mnv_count", "alternative_codons"]

#: 4**3 - 1 alternatives, minus the 3 x 3 single-nucleotide ones.
MNV_PER_CODON = 54


def mnv_count(transcript: Transcript) -> int:
    """Closed form: 54 per codon, terminator codon included."""
    return MNV_PER_CODON * transcript.n_codons


def alternative_codons(ref_codon: str) -> list[tuple[str, tuple[int, ...]]]:
    """Every multi-nucleotide alternative to ``ref_codon``.

    Returns ``(alt_codon, differing_offsets)`` in lexicographic codon order, so
    the enumeration is stable.
    """
    out: list[tuple[str, tuple[int, ...]]] = []
    for combo in product(BASES, repeat=3):
        alt = "".join(combo)
        differing = tuple(i for i in range(3) if alt[i] != ref_codon[i])
        if len(differing) >= 2:
            out.append((alt, differing))
    return out


def enumerate_intracodon_mnvs(transcript: Transcript) -> Iterator[Variant]:
    """Yield every intra-codon MNV, codon by codon."""
    for codon_index in range(1, transcript.n_codons + 1):
        ref_codon = transcript.codon_sequence(codon_index)
        first_cds, _, _ = transcript.codon_bounds_cds(codon_index)
        for alt_codon, differing in alternative_codons(ref_codon):
            edits = tuple(
                CodingEdit(
                    cds_position=first_cds + offset,
                    ref=ref_codon[offset],
                    alt=alt_codon[offset],
                )
                for offset in differing
            )
            yield annotate_coding_edits(transcript, edits, kind=VariantKind.MNV)
