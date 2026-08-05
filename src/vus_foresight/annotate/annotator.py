"""Turn a coordinate edit into an annotated :class:`~vus_foresight.variant.Variant`.

The one non-negotiable rule of this module: **consequences are derived by
translating the whole codon after applying every edit**, never by composing
per-position effects. A multi-nucleotide variant inside one codon routinely
produces a third amino acid that neither component substitution produces, and
composing position-wise silently gets it wrong. Spec section 11 calls this "the
most likely error in the entire layer" and it is guarded by a dedicated test.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..genome.sequence import reverse_complement, translate_codon
from ..genome.transcript import CPosition, IntronSpan, Transcript, format_c_position
from ..variant import Consequence, Variant, VariantKind
from . import hgvs

__all__ = [
    "CONSEQUENCE_RANK",
    "CodingEdit",
    "annotate_coding_edits",
    "annotate_intronic_substitution",
    "genomic_deletion",
    "genomic_vcf",
    "most_severe",
]


@dataclass(frozen=True, slots=True, order=True)
class CodingEdit:
    """A single-base substitution at a 1-based CDS position."""

    cds_position: int
    ref: str
    alt: str


#: Severity ranking, most severe first. Follows the ordering used by Ensembl VEP
#: so that "primary consequence" means the same thing here as in every external
#: annotation this project is compared against.
CONSEQUENCE_RANK: dict[Consequence, int] = {
    Consequence.EXON_DELETION: 0,
    Consequence.EXON_DUPLICATION: 1,
    Consequence.SPLICE_ACCEPTOR: 2,
    Consequence.SPLICE_DONOR: 3,
    Consequence.NONSENSE: 4,
    Consequence.FRAMESHIFT: 5,
    Consequence.STOP_LOST: 6,
    Consequence.START_LOST: 7,
    Consequence.INFRAME_INSERTION: 8,
    Consequence.INFRAME_DELETION: 9,
    Consequence.MISSENSE: 10,
    Consequence.SPLICE_REGION: 11,
    Consequence.STOP_RETAINED: 12,
    Consequence.SYNONYMOUS: 13,
    Consequence.INTRONIC: 14,
    Consequence.UTR5: 15,
    Consequence.UTR3: 16,
}


def most_severe(terms: list[Consequence]) -> Consequence:
    """The highest-ranking term of a set."""
    return min(terms, key=lambda t: CONSEQUENCE_RANK[t])


def _order_terms(terms: set[Consequence]) -> tuple[Consequence, ...]:
    return tuple(sorted(terms, key=lambda t: CONSEQUENCE_RANK[t]))


def genomic_vcf(
    transcript: Transcript,
    cds_first: int,
    cds_last: int,
    alt_coding: str,
) -> str | None:
    """VCF-style ``chrom-pos-ref-alt`` on the genomic plus strand.

    Returns ``None`` when the CDS span is not genomically contiguous, which
    happens for a multi-nucleotide change inside a codon that straddles an exon
    junction. Such a change is one variant in ``c.`` space and two records in
    VCF space; emitting a single fabricated span would be wrong, so the field is
    left empty and the row is flagged via ``attributes``.
    """
    tx_first = transcript.cds_position_to_tx(cds_first)
    tx_last = transcript.cds_position_to_tx(cds_last)
    genomic = [transcript.genomic_at(t) for t in range(tx_first, tx_last + 1)]
    expected_step = transcript.step
    for a, b in zip(genomic, genomic[1:]):
        if b - a != expected_step:
            return None

    ref_coding = transcript.cds[cds_first - 1 : cds_last]
    if transcript.strand == "+":
        pos, ref, alt = genomic[0], ref_coding, alt_coding
    else:
        pos = genomic[-1]
        ref = reverse_complement(ref_coding)
        alt = reverse_complement(alt_coding)
    return f"{transcript.chrom}-{pos}-{ref}-{alt}"


def genomic_deletion(transcript: Transcript, cds_first: int, cds_last: int) -> str | None:
    """``chrom-pos-REF-del`` for a deletion of a contiguous CDS span.

    A VCF record for a deletion needs an anchor base one position 5' of the
    deleted span, which for a deletion starting at an exon boundary is intronic
    and therefore not in the transcript sequence. Rather than fabricate it, the
    deleted span is reported directly with an explicit ``del`` alt, and the
    convention is documented on :class:`~vus_foresight.variant.Variant`.
    """
    tx_first = transcript.cds_position_to_tx(cds_first)
    tx_last = transcript.cds_position_to_tx(cds_last)
    genomic = [transcript.genomic_at(t) for t in range(tx_first, tx_last + 1)]
    for a, b in zip(genomic, genomic[1:]):
        if b - a != transcript.step:
            return None
    ref_coding = transcript.cds[cds_first - 1 : cds_last]
    if transcript.strand == "+":
        return f"{transcript.chrom}-{genomic[0]}-{ref_coding}-del"
    return f"{transcript.chrom}-{genomic[-1]}-{reverse_complement(ref_coding)}-del"


def _coding_consequence(
    transcript: Transcript, codon_index: int, ref_aa: str, alt_aa: str
) -> Consequence:
    if codon_index == transcript.n_codons:
        # The terminator codon itself.
        return Consequence.STOP_RETAINED if alt_aa == "*" else Consequence.STOP_LOST
    if codon_index == 1 and ref_aa == "M" and alt_aa != "M":
        return Consequence.START_LOST
    if alt_aa == "*":
        return Consequence.NONSENSE
    if alt_aa == ref_aa:
        return Consequence.SYNONYMOUS
    return Consequence.MISSENSE


def _protein_hgvs(
    transcript: Transcript,
    codon_index: int,
    ref_aa: str,
    alt_aa: str,
    consequence: Consequence,
) -> str:
    if consequence is Consequence.START_LOST:
        return hgvs.format_protein_start_lost()
    if consequence is Consequence.STOP_LOST:
        return hgvs.format_protein_stop_lost(codon_index, alt_aa)
    if consequence in (Consequence.SYNONYMOUS, Consequence.STOP_RETAINED):
        return hgvs.format_protein_synonymous(ref_aa, codon_index)
    return hgvs.format_protein_substitution(ref_aa, codon_index, alt_aa)


def annotate_coding_edits(
    transcript: Transcript,
    edits: tuple[CodingEdit, ...],
    *,
    kind: VariantKind | None = None,
) -> Variant:
    """Annotate one or more substitutions confined to a single codon.

    ``edits`` must be non-empty, must not repeat a CDS position, must lie in one
    codon, and each ``alt`` must differ from its ``ref``. Those are enforced
    rather than assumed: a no-op "variant" or a cross-codon edit reaching this
    function means the enumerator is broken, and silently annotating it would
    hide that.
    """
    if not edits:
        raise ValueError("at least one edit is required")
    ordered = tuple(sorted(edits))
    positions = [e.cds_position for e in ordered]
    if len(set(positions)) != len(positions):
        raise ValueError(f"duplicate CDS positions in edits: {positions}")

    codon_indices = {(p - 1) // 3 + 1 for p in positions}
    if len(codon_indices) != 1:
        raise ValueError(
            f"edits span {len(codon_indices)} codons ({sorted(codon_indices)}); "
            "annotate_coding_edits handles intra-codon changes only"
        )
    codon_index = codon_indices.pop()
    codon_first_cds, _, _ = transcript.codon_bounds_cds(codon_index)

    ref_codon = transcript.codon_sequence(codon_index)
    codon_chars = list(ref_codon)
    for edit in ordered:
        offset = edit.cds_position - codon_first_cds
        observed = ref_codon[offset]
        if observed != edit.ref:
            raise ValueError(
                f"reference mismatch at CDS {edit.cds_position}: config says "
                f"{observed!r}, edit says {edit.ref!r}"
            )
        if edit.alt == edit.ref:
            raise ValueError(f"no-op edit at CDS {edit.cds_position}: {edit.ref}>{edit.alt}")
        codon_chars[offset] = edit.alt
    alt_codon = "".join(codon_chars)

    # The whole point: translate the assembled codon, do not compose effects.
    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)

    consequence = _coding_consequence(transcript, codon_index, ref_aa, alt_aa)
    terms = {consequence}

    tx_positions = [transcript.cds_position_to_tx(p) for p in positions]
    if any(t in transcript.exonic_splice_region_tx for t in tx_positions):
        terms.add(Consequence.SPLICE_REGION)

    ordered_terms = _order_terms(terms)
    primary = ordered_terms[0]

    # Minimal differing span, per HGVS: trim identical flanking bases.
    differing = [i for i in range(3) if ref_codon[i] != alt_codon[i]]
    span_first_offset, span_last_offset = differing[0], differing[-1]
    cds_first = codon_first_cds + span_first_offset
    cds_last = codon_first_cds + span_last_offset
    c_first = format_c_position(CPosition(base=cds_first))
    c_last = format_c_position(CPosition(base=cds_last))
    alt_span = alt_codon[span_first_offset : span_last_offset + 1]
    ref_span = ref_codon[span_first_offset : span_last_offset + 1]

    if cds_first == cds_last:
        hgvs_c = hgvs.format_substitution(c_first, ref_span, alt_span)
    else:
        hgvs_c = hgvs.format_delins(c_first, c_last, alt_span)

    genomic = genomic_vcf(transcript, cds_first, cds_last, alt_span)
    attributes: dict[str, str] = {}
    if genomic is None:
        attributes["genomic_discontiguous"] = "true"

    distance = len(differing)
    resolved_kind = kind or (VariantKind.SNV if distance == 1 else VariantKind.MNV)

    return Variant(
        gene=transcript.gene,
        transcript=transcript.transcript_id,
        kind=resolved_kind,
        hgvs_c=hgvs_c,
        hgvs_p=_protein_hgvs(transcript, codon_index, ref_aa, alt_aa, consequence),
        grch38_pos=genomic,
        consequence=primary,
        consequence_terms=ordered_terms,
        mutational_distance=distance,
        codon_index=codon_index,
        cds_position=cds_first,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        ptc_codon=codon_index if consequence is Consequence.NONSENSE else None,
        attributes=attributes,
    )


def _intronic_consequence(offset: int) -> Consequence:
    magnitude = abs(offset)
    if magnitude <= 2:
        return Consequence.SPLICE_DONOR if offset > 0 else Consequence.SPLICE_ACCEPTOR
    if magnitude <= 8:
        return Consequence.SPLICE_REGION
    return Consequence.INTRONIC


def annotate_intronic_substitution(
    transcript: Transcript,
    intron: IntronSpan,
    n: int,
    ref: str,
    alt: str,
) -> Variant:
    """Annotate a substitution at the ``n``-th base of an intron (donor-side count).

    ``ref`` is on the coding strand. Callers that have no flanking sequence pass
    ``"N"``; the resulting row carries ``reference_base_unknown`` so that no
    downstream consumer mistakes it for a called reference allele.
    """
    c_pos = transcript.intron_c_position(intron, n)
    consequence = _intronic_consequence(c_pos.offset)
    genomic_pos = transcript.intron_genomic(intron, n)

    if transcript.strand == "+":
        g_ref, g_alt = ref, alt
    else:
        g_ref, g_alt = reverse_complement(ref), reverse_complement(alt)

    attributes: dict[str, str] = {
        "intron_index": str(intron.upstream_index),
        "intron_offset": str(c_pos.offset),
    }
    if ref == "N":
        attributes["reference_base_unknown"] = "true"

    return Variant(
        gene=transcript.gene,
        transcript=transcript.transcript_id,
        kind=VariantKind.SNV,
        hgvs_c=hgvs.format_substitution(format_c_position(c_pos), ref, alt),
        hgvs_p=None,
        grch38_pos=f"{transcript.chrom}-{genomic_pos}-{g_ref}-{g_alt}",
        consequence=consequence,
        consequence_terms=(consequence,),
        mutational_distance=1,
        codon_index=None,
        cds_position=None,
        attributes=attributes,
    )
