"""Tier 2: in-frame deletions of one and two codons.

Directly enumerable and relevant to PM4. Every 3-nt and 6-nt deletion is
generated at every CDS offset -- not only codon-aligned ones, because a deletion
straddling two codons fuses them into a new residue and is a different protein
change from either aligned neighbour.

After 3'-most normalisation many offsets collapse onto the same description, so
the emitted set is deduplicated by normalised HGVS.
"""

from __future__ import annotations

from typing import Iterator

from ..annotate.indel import annotate_inframe_deletion
from ..genome.transcript import Transcript
from ..variant import Variant

__all__ = ["INFRAME_DELETION_LENGTHS", "enumerate_inframe_deletions"]

#: One and two codons.
INFRAME_DELETION_LENGTHS: tuple[int, ...] = (3, 6)


def enumerate_inframe_deletions(
    transcript: Transcript, *, lengths: tuple[int, ...] = INFRAME_DELETION_LENGTHS
) -> Iterator[Variant]:
    """Yield every distinct in-frame deletion of the given lengths.

    The terminator codon is excluded from the span: removing it is a stop-loss,
    which is a different consequence answering to different criteria.
    """
    coding_limit = transcript.cds_length - 3
    seen: set[str] = set()
    for length in lengths:
        for start in range(1, coding_limit - length + 2):
            variant = annotate_inframe_deletion(transcript, start, length)
            if variant.hgvs_c in seen:
                continue
            seen.add(variant.hgvs_c)
            yield variant
