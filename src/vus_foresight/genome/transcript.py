"""Transcript structure and coordinate arithmetic.

This is the highest-risk module in the project. Spec section 11 says it plainly:
a bug in enumeration produces variants that do not exist and is caught by
counting, but a bug *here* produces HGVS that looks right and is wrong.

Two rules keep it honest:

* MANE Select is the only source of truth. Legacy BIC numbering never enters
  here, not even as a fallback.
* Strand is handled in exactly one place -- :meth:`Transcript.genomic_at` -- and
  every other coordinate function is expressed in transcript space. That is what
  makes the BRCA1(-)/BRCA2(+) symmetry test a real control rather than a
  coincidence.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from functools import cached_property
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .sequence import reverse_complement

__all__ = [
    "CPosition",
    "Exon",
    "IntronSpan",
    "Transcript",
    "format_c_position",
]


class Exon(BaseModel):
    """One exon in genomic coordinates.

    ``start``/``end`` are 1-based inclusive and always ``start <= end``, i.e.
    genomic order, independent of strand. ``label`` is the *clinical* exon name
    and is a string on purpose: BRCA1 has no exon 4, so its labels run
    1,2,3,5,...,24 and arithmetic on them is meaningless.
    """

    model_config = ConfigDict(frozen=True)

    label: str
    start: int = Field(gt=0)
    end: int = Field(gt=0)

    @model_validator(mode="after")
    def _ordered(self) -> Exon:
        if self.end < self.start:
            raise ValueError(f"exon {self.label}: end {self.end} < start {self.start}")
        return self

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True, slots=True, order=True)
class CPosition:
    """An HGVS ``c.`` coordinate.

    ``base`` is the anchor: positive inside the CDS, negative in the 5' UTR, and
    positive-with-``utr3`` in the 3' UTR (rendered ``*n``). ``offset`` is 0 for
    exonic positions and non-zero for intronic ones (``c.100+3``, ``c.101-2``).
    """

    base: int
    offset: int = 0
    utr3: bool = False

    @property
    def is_intronic(self) -> bool:
        return self.offset != 0

    @property
    def is_coding(self) -> bool:
        return self.offset == 0 and not self.utr3 and self.base > 0

    def __str__(self) -> str:  # pragma: no cover - trivial delegation
        return format_c_position(self)


def format_c_position(pos: CPosition) -> str:
    """Render a :class:`CPosition` in HGVS form."""
    anchor = f"*{pos.base}" if pos.utr3 else str(pos.base)
    if pos.offset == 0:
        return anchor
    return f"{anchor}{pos.offset:+d}"


@dataclass(frozen=True, slots=True)
class IntronSpan:
    """One intron, described in transcript order.

    ``upstream_index``/``downstream_index`` are 0-based positions in the
    transcript-ordered exon list.
    """

    upstream_index: int
    downstream_index: int
    #: Genomic coordinate of the last exonic base *before* the intron, in
    #: transcript direction.
    upstream_last_genomic: int
    #: Genomic coordinate of the first exonic base *after* the intron.
    downstream_first_genomic: int
    length: int
    #: ``c.`` anchor for ``+n`` numbering (last base of the upstream exon).
    upstream_anchor: CPosition
    #: ``c.`` anchor for ``-n`` numbering (first base of the downstream exon).
    downstream_anchor: CPosition


class Transcript(BaseModel):
    """A MANE Select transcript with its spliced sequence.

    ``exons`` are given in *transcript* order (5' to 3' of the mRNA), so for a
    minus-strand gene their genomic coordinates descend. ``sequence`` is the
    spliced transcript on the coding strand, so it can be indexed directly by
    transcript position without any strand reasoning.
    """

    model_config = ConfigDict(frozen=True)

    transcript_id: str
    gene: str
    chrom: str
    strand: Literal["+", "-"]
    exons: tuple[Exon, ...]
    #: 1-based transcript position of the A of the initiation codon.
    cds_start_tx: int = Field(gt=0)
    #: 1-based transcript position of the last base of the termination codon.
    cds_end_tx: int = Field(gt=0)
    sequence: str

    # ------------------------------------------------------------------
    # validation
    # ------------------------------------------------------------------
    @model_validator(mode="after")
    def _check(self) -> Transcript:
        if len(self.exons) == 0:
            raise ValueError("transcript must have at least one exon")

        total = sum(e.length for e in self.exons)
        if len(self.sequence) != total:
            raise ValueError(
                f"{self.transcript_id}: sequence length {len(self.sequence)} does not "
                f"match summed exon length {total}"
            )

        # Exons must be non-overlapping and ordered along the transcript.
        for prev, nxt in zip(self.exons, self.exons[1:]):
            if self.strand == "+":
                if nxt.start <= prev.end:
                    raise ValueError(
                        f"{self.transcript_id}: exons {prev.label}/{nxt.label} overlap "
                        "or are out of transcript order on the plus strand"
                    )
            else:
                if nxt.end >= prev.start:
                    raise ValueError(
                        f"{self.transcript_id}: exons {prev.label}/{nxt.label} overlap "
                        "or are out of transcript order on the minus strand"
                    )

        if self.cds_end_tx <= self.cds_start_tx:
            raise ValueError(f"{self.transcript_id}: CDS end is not after CDS start")
        if self.cds_end_tx > total:
            raise ValueError(f"{self.transcript_id}: CDS end runs past the transcript")
        cds_len = self.cds_end_tx - self.cds_start_tx + 1
        if cds_len % 3 != 0:
            raise ValueError(f"{self.transcript_id}: CDS length {cds_len} is not a multiple of 3")

        labels = [e.label for e in self.exons]
        if len(set(labels)) != len(labels):
            raise ValueError(f"{self.transcript_id}: duplicate exon labels {labels}")
        return self

    # ------------------------------------------------------------------
    # basic geometry
    # ------------------------------------------------------------------
    @property
    def length(self) -> int:
        """Spliced transcript length in nucleotides."""
        return len(self.sequence)

    @property
    def cds_length(self) -> int:
        """CDS length including the termination codon."""
        return self.cds_end_tx - self.cds_start_tx + 1

    @property
    def n_codons(self) -> int:
        """Number of codons including the terminator."""
        return self.cds_length // 3

    @property
    def protein_length(self) -> int:
        """Number of translated residues, excluding the terminator."""
        return self.n_codons - 1

    @cached_property
    def cds(self) -> str:
        """The spliced coding sequence, initiation codon through terminator."""
        return self.sequence[self.cds_start_tx - 1 : self.cds_end_tx]

    @property
    def step(self) -> int:
        """Genomic increment corresponding to +1 in transcript space."""
        return 1 if self.strand == "+" else -1

    def _exon_first_genomic(self, exon: Exon) -> int:
        return exon.start if self.strand == "+" else exon.end

    def _exon_last_genomic(self, exon: Exon) -> int:
        return exon.end if self.strand == "+" else exon.start

    @cached_property
    def _exon_tx_bounds(self) -> tuple[tuple[int, int], ...]:
        """``(tx_start, tx_end)`` 1-based inclusive for each exon, in tx order."""
        bounds: list[tuple[int, int]] = []
        cursor = 1
        for exon in self.exons:
            bounds.append((cursor, cursor + exon.length - 1))
            cursor += exon.length
        return tuple(bounds)

    # ------------------------------------------------------------------
    # transcript <-> genomic
    # ------------------------------------------------------------------
    def genomic_at(self, tx_pos: int) -> int:
        """Genomic coordinate of a 1-based transcript position.

        The single point in the codebase where strand becomes arithmetic.
        """
        if not 1 <= tx_pos <= self.length:
            raise ValueError(
                f"transcript position {tx_pos} outside 1..{self.length} for {self.transcript_id}"
            )
        for exon, (tx_start, tx_end) in zip(self.exons, self._exon_tx_bounds):
            if tx_start <= tx_pos <= tx_end:
                return self._exon_first_genomic(exon) + self.step * (tx_pos - tx_start)
        raise AssertionError("unreachable: exon bounds do not cover the transcript")

    @cached_property
    def _genomic_to_tx(self) -> dict[int, int]:
        return {self.genomic_at(t): t for t in range(1, self.length + 1)}

    def tx_at_genomic(self, genomic_pos: int) -> int | None:
        """Transcript position for a genomic coordinate, or ``None`` if intronic."""
        return self._genomic_to_tx.get(genomic_pos)

    def exon_index_at_tx(self, tx_pos: int) -> int:
        """0-based index, in transcript order, of the exon containing ``tx_pos``."""
        for idx, (tx_start, tx_end) in enumerate(self._exon_tx_bounds):
            if tx_start <= tx_pos <= tx_end:
                return idx
        raise ValueError(f"transcript position {tx_pos} is out of range")

    # ------------------------------------------------------------------
    # transcript <-> c.
    # ------------------------------------------------------------------
    def c_at_tx(self, tx_pos: int) -> CPosition:
        """``c.`` coordinate of an exonic transcript position."""
        if tx_pos < self.cds_start_tx:
            return CPosition(base=tx_pos - self.cds_start_tx)
        if tx_pos > self.cds_end_tx:
            return CPosition(base=tx_pos - self.cds_end_tx, utr3=True)
        return CPosition(base=tx_pos - self.cds_start_tx + 1)

    def tx_at_c(self, pos: CPosition) -> int:
        """Transcript position of an *exonic* ``c.`` coordinate."""
        if pos.offset != 0:
            raise ValueError(f"c.{format_c_position(pos)} is intronic, not exonic")
        if pos.utr3:
            return self.cds_end_tx + pos.base
        if pos.base < 0:
            return self.cds_start_tx + pos.base
        if pos.base == 0:
            raise ValueError("c.0 does not exist")
        return self.cds_start_tx + pos.base - 1

    def cds_position_to_tx(self, cds_pos: int) -> int:
        """Transcript position of a 1-based CDS offset."""
        if not 1 <= cds_pos <= self.cds_length:
            raise ValueError(f"CDS position {cds_pos} outside 1..{self.cds_length}")
        return self.cds_start_tx + cds_pos - 1

    def genomic_at_cds(self, cds_pos: int) -> int:
        """Genomic coordinate of a 1-based CDS offset."""
        return self.genomic_at(self.cds_position_to_tx(cds_pos))

    # ------------------------------------------------------------------
    # introns
    # ------------------------------------------------------------------
    @cached_property
    def introns(self) -> tuple[IntronSpan, ...]:
        spans: list[IntronSpan] = []
        bounds = self._exon_tx_bounds
        for i in range(len(self.exons) - 1):
            up, down = self.exons[i], self.exons[i + 1]
            up_last = self._exon_last_genomic(up)
            down_first = self._exon_first_genomic(down)
            length = abs(down_first - up_last) - 1
            if length <= 0:
                raise ValueError(
                    f"{self.transcript_id}: exons {up.label}/{down.label} are "
                    "contiguous; an intron of length 0 is not representable"
                )
            spans.append(
                IntronSpan(
                    upstream_index=i,
                    downstream_index=i + 1,
                    upstream_last_genomic=up_last,
                    downstream_first_genomic=down_first,
                    length=length,
                    upstream_anchor=self.c_at_tx(bounds[i][1]),
                    downstream_anchor=self.c_at_tx(bounds[i + 1][0]),
                )
            )
        return tuple(spans)

    def intron_c_position(self, intron: IntronSpan, n: int) -> CPosition:
        """``c.`` coordinate of the ``n``-th intronic base, counting from the donor.

        Applies the HGVS numbering split: bases in the 5' half are numbered
        ``+n`` from the upstream exon, bases in the 3' half ``-m`` from the
        downstream exon. For an odd-length intron the central base belongs to
        the 5' half, per HGVS.
        """
        if not 1 <= n <= intron.length:
            raise ValueError(f"intronic offset {n} outside 1..{intron.length}")
        if n <= (intron.length + 1) // 2:
            anchor = intron.upstream_anchor
            return CPosition(base=anchor.base, offset=n, utr3=anchor.utr3)
        anchor = intron.downstream_anchor
        return CPosition(base=anchor.base, offset=-(intron.length - n + 1), utr3=anchor.utr3)

    def intron_genomic(self, intron: IntronSpan, n: int) -> int:
        """Genomic coordinate of the ``n``-th intronic base, counting from the donor."""
        if not 1 <= n <= intron.length:
            raise ValueError(f"intronic offset {n} outside 1..{intron.length}")
        return intron.upstream_last_genomic + self.step * n

    # ------------------------------------------------------------------
    # zoning used by consequence assignment
    # ------------------------------------------------------------------
    @cached_property
    def exonic_splice_region_tx(self) -> frozenset[int]:
        """Transcript positions within 3 nt of an exon edge that abuts an intron."""
        positions: set[int] = set()
        bounds = self._exon_tx_bounds
        for idx, (tx_start, tx_end) in enumerate(bounds):
            if idx > 0:  # acceptor side
                positions.update(range(tx_start, min(tx_start + 3, tx_end + 1)))
            if idx < len(bounds) - 1:  # donor side
                positions.update(range(max(tx_end - 2, tx_start), tx_end + 1))
        return frozenset(positions)

    @cached_property
    def nmd_escape_cds_start(self) -> int:
        """First CDS position from which a PTC escapes nonsense-mediated decay.

        The standard 50-nucleotide rule: a termination codon in the last exon,
        or within the last 50 nt of the penultimate exon, is not recognised by
        the EJC-dependent NMD machinery.

        Returned as a 1-based CDS offset. For a single-exon transcript there is
        no EJC at all, so nothing is degraded and the boundary is CDS position 1.
        """
        bounds = self._exon_tx_bounds
        if len(bounds) < 2:
            return 1
        penultimate_end_tx = bounds[-2][1]
        boundary_tx = penultimate_end_tx - 50 + 1
        boundary_cds = boundary_tx - self.cds_start_tx + 1
        return max(1, boundary_cds)

    def ptc_escapes_nmd(self, ptc_codon: int) -> bool:
        """Whether a PTC at the given 1-based codon escapes NMD."""
        first_cds_pos, _, _ = self.codon_bounds_cds(ptc_codon)
        return first_cds_pos >= self.nmd_escape_cds_start

    def codon_bounds_cds(self, codon_index: int) -> tuple[int, int, int]:
        """The three 1-based CDS positions of a 1-based codon index."""
        if not 1 <= codon_index <= self.n_codons:
            raise ValueError(f"codon {codon_index} outside 1..{self.n_codons}")
        first = (codon_index - 1) * 3 + 1
        return (first, first + 1, first + 2)

    def codon_sequence(self, codon_index: int) -> str:
        first, _, last = self.codon_bounds_cds(codon_index)
        return self.cds[first - 1 : last]

    def iter_cds_positions(self) -> Iterator[tuple[int, str]]:
        """Yield ``(cds_position, ref_base)`` for the whole CDS, in order."""
        for i, base in enumerate(self.cds, start=1):
            yield i, base


class FlankSequences(BaseModel):
    """Intronic flanking sequence, keyed by genomic coordinate.

    Kept separate from :class:`Transcript` because the transcript sequence is a
    spliced mRNA and has no intronic bases at all. Supplying flanks is optional:
    without them the intronic enumeration still produces every position (its
    cardinality does not depend on the reference base) but cannot name a
    reference allele, and such rows are emitted with ``ref = "N"`` and marked so
    downstream.
    """

    model_config = ConfigDict(frozen=True)

    #: genomic position -> plus-strand reference base
    bases: dict[int, str] = Field(default_factory=dict)

    def base_at(self, transcript: Transcript, intron: IntronSpan, n: int) -> str:
        genomic = transcript.intron_genomic(intron, n)
        plus_base = self.bases.get(genomic)
        if plus_base is None:
            return "N"
        return plus_base if transcript.strand == "+" else reverse_complement(plus_base)

    @property
    def is_empty(self) -> bool:
        return not self.bases
