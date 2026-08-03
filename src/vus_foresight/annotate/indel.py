"""Frameshift and in-frame deletion machinery.

Two things live here that the substitution annotator does not need:

* :func:`first_stop_index` -- an O(n) dynamic program giving, for every offset
  into the coding+3'UTR sequence, the first in-frame terminator reachable from
  it. Every frameshift query then costs O(1) instead of re-translating the
  downstream sequence, which is what makes enumerating a hundred thousand
  minimal indels per gene tractable.
* 3'-most normalisation. HGVS requires a deletion in a repeated stretch to be
  shifted as far 3' as possible; skipping this produces descriptions that are
  arithmetically defensible and do not match any curated record.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..genome.sequence import translate, translate_codon
from ..genome.transcript import CPosition, Transcript, format_c_position
from ..variant import Consequence, Variant, VariantKind
from . import hgvs
from .annotator import CONSEQUENCE_RANK, genomic_deletion

__all__ = [
    "StopIndex",
    "build_stop_index",
    "first_stop_index",
    "ptc_after_indel",
    "normalize_deletion_3prime",
    "protein_change",
    "ProteinChange",
    "annotate_inframe_deletion",
    "coding_and_utr3",
]


def coding_and_utr3(transcript: Transcript) -> str:
    """CDS plus 3' UTR, i.e. everything a shifted reading frame can run into.

    A frameshift near the 3' end frequently reads past the normal terminator, so
    stopping the search at the end of the CDS would report "no stop" for
    variants that plainly have one.
    """
    return transcript.sequence[transcript.cds_start_tx - 1 :]


@dataclass(frozen=True, slots=True)
class StopIndex:
    """Precomputed first-terminator lookup over one sequence."""

    sequence: str
    table: tuple[int | None, ...]

    def __getitem__(self, offset: int) -> int | None:
        if offset < 0 or offset >= len(self.table):
            return None
        return self.table[offset]


def build_stop_index(sequence: str) -> StopIndex:
    """``table[j]`` is the smallest ``s >= j`` with ``s % 3 == j % 3`` and a stop.

    Because the recurrence steps by three, the table is frame-agnostic: it
    answers the question for whichever reading frame the query offset implies.
    """
    n = len(sequence)
    table: list[int | None] = [None] * n
    for j in range(n - 1, -1, -1):
        if j + 3 > n:
            table[j] = None
            continue
        if translate_codon(sequence[j : j + 3]) == "*":
            table[j] = j
        else:
            table[j] = table[j + 3] if j + 3 < n else None
    return StopIndex(sequence=sequence, table=tuple(table))


def first_stop_index(sequence: str, offset: int) -> int | None:
    """Convenience wrapper for a one-off query."""
    return build_stop_index(sequence)[offset]


def ptc_after_indel(
    stop_index: StopIndex,
    cds_position: int,
    deleted: int,
    inserted: str = "",
) -> int | None:
    """1-based codon number of the first terminator in the mutant protein.

    ``cds_position`` is the 1-based CDS coordinate of the first deleted base
    (for a pure insertion, the base *after* which the insertion occurs, plus
    one). ``None`` means the shifted frame runs off the end of the available
    sequence without meeting a terminator.

    The walk is short: at most a couple of hybrid codons straddling the edit
    have to be assembled by hand, after which the remaining mutant sequence is
    literally a suffix of the reference and :class:`StopIndex` answers in O(1).
    """
    sequence = stop_index.sequence
    if cds_position < 1:
        raise ValueError("cds_position is 1-based and must be >= 1")

    codon_index0 = (cds_position - 1) // 3
    head = sequence[3 * codon_index0 : cds_position - 1] + inserted
    tail_start = cds_position - 1 + deleted

    k = codon_index0  # 0-based index of the next mutant codon to read
    i_head = 0
    i_tail = 0
    while True:
        if i_head >= len(head):
            base = tail_start + i_tail
            stop = stop_index[base]
            if stop is None:
                return None
            return k + 1 + (stop - base) // 3

        chars: list[str] = []
        for _ in range(3):
            if i_head < len(head):
                chars.append(head[i_head])
                i_head += 1
            else:
                pos = tail_start + i_tail
                if pos >= len(sequence):
                    return None
                chars.append(sequence[pos])
                i_tail += 1
        if translate_codon("".join(chars)) == "*":
            return k + 1
        k += 1


def normalize_deletion_3prime(sequence: str, start: int, length: int, *, limit: int) -> int:
    """Shift a deletion as far 3' as the sequence allows.

    ``start`` and the return value are 1-based offsets into ``sequence``;
    ``limit`` is the last 1-based position the deletion is permitted to occupy,
    which keeps a deletion from sliding out of the region being enumerated.
    """
    if length <= 0:
        raise ValueError("deletion length must be positive")
    while start + length - 1 < limit and sequence[start - 1] == sequence[start + length - 1]:
        start += 1
    return start


@dataclass(frozen=True, slots=True)
class ProteinChange:
    """The trimmed difference between a reference and an alternate peptide."""

    #: 1-based residue number of the first changed reference residue.
    first_residue: int
    #: The reference residues removed or replaced (may be empty for an insertion).
    ref_segment: str
    #: The residues that take their place (empty for a pure deletion).
    alt_segment: str

    @property
    def last_residue(self) -> int:
        return self.first_residue + max(len(self.ref_segment), 1) - 1

    @property
    def is_pure_deletion(self) -> bool:
        return self.alt_segment == "" and self.ref_segment != ""

    @property
    def is_silent(self) -> bool:
        return self.ref_segment == "" and self.alt_segment == ""


def protein_change(ref_protein: str, alt_protein: str) -> ProteinChange:
    """Trim the common prefix, then the common suffix.

    Prefix first is what makes the result 3'-most, which is the HGVS
    requirement: ``MAAA*`` -> ``MAA*`` must be reported as a deletion of the
    third alanine, not the first.
    """
    prefix = 0
    limit = min(len(ref_protein), len(alt_protein))
    while prefix < limit and ref_protein[prefix] == alt_protein[prefix]:
        prefix += 1

    suffix = 0
    while (
        suffix < limit - prefix
        and ref_protein[len(ref_protein) - 1 - suffix] == alt_protein[len(alt_protein) - 1 - suffix]
    ):
        suffix += 1

    return ProteinChange(
        first_residue=prefix + 1,
        ref_segment=ref_protein[prefix : len(ref_protein) - suffix],
        alt_segment=alt_protein[prefix : len(alt_protein) - suffix],
    )


def annotate_inframe_deletion(
    transcript: Transcript,
    cds_start: int,
    length: int,
) -> Variant:
    """Annotate an in-frame deletion of ``length`` (3 or 6) nucleotides.

    The deletion is normalised 3'-most within the region that excludes the
    terminator codon. Deletions overlapping the initiation codon are reported as
    ``start_lost``: their protein effect is not predictable and inventing one
    would be worse than saying so.
    """
    if length % 3 != 0:
        raise ValueError(f"in-frame deletion length must be a multiple of 3, got {length}")
    cds = transcript.cds
    # The terminator codon is excluded: deleting it is a stop-loss event, a
    # different class of variant scored by different criteria.
    coding_limit = transcript.cds_length - 3
    if cds_start < 1 or cds_start + length - 1 > coding_limit:
        raise ValueError(
            f"in-frame deletion {cds_start}..{cds_start + length - 1} falls outside "
            f"CDS positions 1..{coding_limit}"
        )

    start = normalize_deletion_3prime(cds, cds_start, length, limit=coding_limit)
    end = start + length - 1

    alt_cds = cds[: start - 1] + cds[start - 1 + length :]
    ref_protein = translate(cds, stop_at_terminator=True)
    alt_protein = translate(alt_cds, stop_at_terminator=True)
    change = protein_change(ref_protein, alt_protein)

    c_start = format_c_position(CPosition(base=start))
    c_end = format_c_position(CPosition(base=end))
    hgvs_c = hgvs.format_deletion(c_start, c_end)

    if start <= 3:
        consequence = Consequence.START_LOST
        hgvs_p: str | None = hgvs.format_protein_start_lost()
    elif change.is_silent:
        # Deleting a repeat unit can leave the peptide unchanged in principle;
        # for a genuine in-frame deletion this should not happen, so it is worth
        # surfacing rather than swallowing.
        consequence = Consequence.INFRAME_DELETION
        hgvs_p = None
    elif change.is_pure_deletion:
        consequence = Consequence.INFRAME_DELETION
        hgvs_p = hgvs.format_protein_deletion(
            change.ref_segment[0],
            change.first_residue,
            change.ref_segment[-1],
            change.first_residue + len(change.ref_segment) - 1,
        )
    else:
        consequence = Consequence.INFRAME_DELETION
        hgvs_p = hgvs.format_protein_delins(
            change.ref_segment[0],
            change.first_residue,
            change.ref_segment[-1],
            change.first_residue + len(change.ref_segment) - 1,
            change.alt_segment,
        )

    terms = {consequence}
    tx_positions = [transcript.cds_position_to_tx(p) for p in range(start, end + 1)]
    if any(t in transcript.exonic_splice_region_tx for t in tx_positions):
        terms.add(Consequence.SPLICE_REGION)
    ordered = tuple(sorted(terms, key=lambda t: CONSEQUENCE_RANK[t]))

    genomic = genomic_deletion(transcript, start, end)
    attributes = {"deleted_length": str(length)}
    if genomic is None:
        attributes["genomic_discontiguous"] = "true"

    return Variant(
        gene=transcript.gene,
        transcript=transcript.transcript_id,
        kind=VariantKind.INFRAME_DELETION,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        grch38_pos=genomic,
        consequence=ordered[0],
        consequence_terms=ordered,
        mutational_distance=length,
        codon_index=(start - 1) // 3 + 1,
        cds_position=start,
        attributes=attributes,
    )
