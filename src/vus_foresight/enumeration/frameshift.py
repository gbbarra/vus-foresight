"""Tier 2: frameshifts, enumerated as equivalence classes of resulting PTC.

Indels are not enumerated -- there are unboundedly many and they are not what
PVS1 reads. What PVS1 reads is *where the premature termination codon lands*,
so the enumeration unit is the PTC position, and each class carries the minimal
indels that produce it.

Volume collapses from "all indels" to one class per residue.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator

from ..annotate import hgvs
from ..annotate.indel import (
    StopIndex,
    build_stop_index,
    coding_and_utr3,
    normalize_deletion_3prime,
    ptc_after_indel,
)
from ..genome.sequence import BASES, translate, translate_codon
from ..genome.transcript import CPosition, Transcript, format_c_position
from ..variant import Consequence, Variant, VariantKind

__all__ = [
    "MinimalIndel",
    "FrameshiftClass",
    "FrameshiftEnumeration",
    "normalize_insertion_3prime",
    "enumerate_minimal_indels",
    "build_frameshift_classes",
    "enumerate_frameshift_classes",
    "MAX_REPRESENTATIVES",
]

#: How many member indels each class stores. The full membership is a set of
#: known size; listing all of it would bloat every row without adding
#: information, but a class must be able to prove itself, so a handful of
#: verifiable representatives are kept alongside the exact count.
MAX_REPRESENTATIVES = 5


@dataclass(frozen=True, slots=True, order=True)
class MinimalIndel:
    """A smallest edit that shifts the reading frame.

    Ordering is by edit size, then position, then deleted-before-inserted, then
    inserted sequence -- so "the representative" of a class is well defined and
    stable across runs.
    """

    size: int
    cds_position: int
    deleted: int
    inserted: str = ""

    @property
    def net_shift(self) -> int:
        return len(self.inserted) - self.deleted


def normalize_insertion_3prime(
    sequence: str, position: int, inserted: str, *, limit: int
) -> tuple[int, str]:
    """Shift an insertion as far 3' as the sequence allows.

    ``position`` is the 1-based offset *before which* the sequence is inserted.
    Sliding rotates the inserted string, which is what keeps
    ``ins`` of ``AG`` before an ``A`` equal to ``ins`` of ``GA`` one base later.
    """
    if not inserted:
        raise ValueError("insertion must be non-empty")
    while position <= limit and sequence[position - 1] == inserted[0]:
        inserted = inserted[1:] + inserted[0]
        position += 1
    return position, inserted


def _deletion_hgvs(transcript: Transcript, indel: MinimalIndel) -> str:
    start = normalize_deletion_3prime(
        transcript.cds, indel.cds_position, indel.deleted, limit=transcript.cds_length
    )
    end = start + indel.deleted - 1
    return hgvs.format_deletion(
        format_c_position(CPosition(base=start)), format_c_position(CPosition(base=end))
    )


def _insertion_hgvs(transcript: Transcript, indel: MinimalIndel) -> str:
    cds = transcript.cds
    position, inserted = normalize_insertion_3prime(
        cds, indel.cds_position, indel.inserted, limit=transcript.cds_length
    )
    length = len(inserted)
    # A duplication is an insertion whose content equals the bases immediately
    # 5' of it; HGVS mandates the ``dup`` form in that case.
    preceding = cds[position - 1 - length : position - 1]
    if preceding == inserted:
        first = position - length
        last = position - 1
        span = (
            format_c_position(CPosition(base=first))
            if first == last
            else f"{format_c_position(CPosition(base=first))}_"
            f"{format_c_position(CPosition(base=last))}"
        )
        return f"c.{span}dup"
    return (
        f"c.{format_c_position(CPosition(base=position - 1))}_"
        f"{format_c_position(CPosition(base=position))}ins{inserted}"
    )


def indel_hgvs(transcript: Transcript, indel: MinimalIndel) -> str:
    """HGVS ``c.`` description of a minimal indel, 3'-normalised."""
    if indel.deleted and not indel.inserted:
        return _deletion_hgvs(transcript, indel)
    if indel.inserted and not indel.deleted:
        return _insertion_hgvs(transcript, indel)
    raise ValueError("minimal indels are pure deletions or pure insertions")


def enumerate_minimal_indels(transcript: Transcript) -> Iterator[MinimalIndel]:
    """Every 1-2 nt deletion and 1-2 nt insertion inside the CDS.

    These are the frame-shifting edits of smallest size. Larger indels reach the
    same PTC classes and add nothing to the map, which is the entire argument
    for classing by PTC in the first place.
    """
    n = transcript.cds_length
    for position in range(1, n + 1):
        yield MinimalIndel(size=1, cds_position=position, deleted=1)
    for position in range(1, n):
        yield MinimalIndel(size=2, cds_position=position, deleted=2)
    # Insertions before position 1 would precede the initiation codon.
    for position in range(2, n + 1):
        for base in BASES:
            yield MinimalIndel(size=1, cds_position=position, deleted=0, inserted=base)
    for position in range(2, n + 1):
        for first in BASES:
            for second in BASES:
                yield MinimalIndel(
                    size=2, cds_position=position, deleted=0, inserted=first + second
                )


@dataclass(slots=True)
class FrameshiftClass:
    """All frameshifts whose premature terminator lands on the same codon."""

    ptc_codon: int
    member_count: int = 0
    representatives: list[MinimalIndel] = field(default_factory=list)

    def add(self, indel: MinimalIndel) -> None:
        self.member_count += 1
        self.representatives.append(indel)
        self.representatives.sort()
        del self.representatives[MAX_REPRESENTATIVES:]


@dataclass(frozen=True, slots=True)
class FrameshiftEnumeration:
    """The classes, plus the codons no minimal indel can reach.

    The two lists together always account for every residue of the protein. A
    non-empty ``unreachable`` is a real property of the sequence, not a gap in
    the enumeration: some codons simply cannot be the *first* stop of any 1-2 nt
    frameshift. Reporting it is what keeps the count honest.
    """

    classes: tuple[FrameshiftClass, ...]
    unreachable: tuple[int, ...]

    @property
    def total_codon_positions(self) -> int:
        return len(self.classes) + len(self.unreachable)


def build_frameshift_classes(transcript: Transcript) -> FrameshiftEnumeration:
    """Group every minimal frame-shifting indel by the PTC it produces."""
    search_space = coding_and_utr3(transcript)
    stop_index = build_stop_index(search_space)

    by_ptc: dict[int, FrameshiftClass] = {}
    for indel in enumerate_minimal_indels(transcript):
        if indel.net_shift % 3 == 0:
            continue
        ptc = ptc_after_indel(
            stop_index, indel.cds_position, indel.deleted, indel.inserted
        )
        if ptc is None or ptc > transcript.protein_length:
            # A terminator at or beyond the normal stop is not premature.
            continue
        by_ptc.setdefault(ptc, FrameshiftClass(ptc_codon=ptc)).add(indel)

    classes = tuple(by_ptc[k] for k in sorted(by_ptc))
    unreachable = tuple(
        codon for codon in range(1, transcript.protein_length + 1) if codon not in by_ptc
    )
    return FrameshiftEnumeration(classes=classes, unreachable=unreachable)


def _frameshift_protein_hgvs(
    transcript: Transcript, indel: MinimalIndel, ptc_codon: int
) -> str:
    """Derive ``p.Arg100SerfsTer12`` by translating the mutant, not by rule."""
    cds = coding_and_utr3(transcript)
    window = 3 * ptc_codon
    mutant = (
        cds[: indel.cds_position - 1]
        + indel.inserted
        + cds[indel.cds_position - 1 + indel.deleted : indel.cds_position - 1 + indel.deleted + window]
    )
    ref_protein = translate(cds[: window + 3])
    alt_protein = translate(mutant[:window])

    for i in range(min(len(ref_protein), len(alt_protein))):
        if ref_protein[i] != alt_protein[i]:
            residue = i + 1
            alt_aa = alt_protein[i]
            if alt_aa == "*":
                return hgvs.format_protein_substitution(ref_protein[i], residue, "*")
            return hgvs.format_protein_frameshift(
                ref_protein[i], residue, alt_aa, ptc_codon - residue + 1
            )
    # Identical up to the PTC: the shift only manifests at the terminator.
    residue = ptc_codon
    return hgvs.format_protein_substitution(ref_protein[residue - 1], residue, "*")


def enumerate_frameshift_classes(transcript: Transcript) -> Iterator[Variant]:
    """Yield one row per reachable PTC position."""
    enumeration = build_frameshift_classes(transcript)
    for cls in enumeration.classes:
        representative = cls.representatives[0]
        yield Variant(
            gene=transcript.gene,
            transcript=transcript.transcript_id,
            kind=VariantKind.FRAMESHIFT_CLASS,
            hgvs_c=indel_hgvs(transcript, representative),
            hgvs_p=_frameshift_protein_hgvs(transcript, representative, cls.ptc_codon),
            grch38_pos=None,
            consequence=Consequence.FRAMESHIFT,
            consequence_terms=(Consequence.FRAMESHIFT,),
            mutational_distance=0,
            codon_index=(representative.cds_position - 1) // 3 + 1,
            cds_position=representative.cds_position,
            ptc_codon=cls.ptc_codon,
            attributes={
                "member_count": str(cls.member_count),
                "representative_indels": ";".join(
                    indel_hgvs(transcript, i) for i in cls.representatives
                ),
                "nmd_escape": str(transcript.ptc_escapes_nmd(cls.ptc_codon)).lower(),
            },
        )


def verify_class_membership(
    transcript: Transcript, indel: MinimalIndel, ptc_codon: int
) -> bool:
    """Independently confirm an indel produces a PTC at ``ptc_codon``.

    Rebuilds the mutant sequence and translates it from scratch -- no
    :class:`StopIndex`, no metadata. This is the oracle the layer-1 test uses,
    and it must not share an implementation with the thing under test.
    """
    cds = coding_and_utr3(transcript)
    mutant = (
        cds[: indel.cds_position - 1]
        + indel.inserted
        + cds[indel.cds_position - 1 + indel.deleted :]
    )
    for codon_number, offset in enumerate(range(0, len(mutant) - 2, 3), start=1):
        if translate_codon(mutant[offset : offset + 3]) == "*":
            return codon_number == ptc_codon
    return False
