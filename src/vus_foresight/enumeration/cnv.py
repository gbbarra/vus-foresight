"""Tier 2: exon-level copy number variants.

Every contiguous run of exons, deleted and duplicated: ``n(n+1)/2`` intervals
per gene per direction. The volume is trivial, but the *scoring* is not the
point system used everywhere else -- copy number changes are scored by the
ClinGen CNV framework (Riggs et al.), so they are produced by this module and
evaluated by :mod:`vus_foresight.engine.cnv_scoring`, then emitted into the same
output schema.
"""

from __future__ import annotations

from collections.abc import Iterator

from ..annotate import hgvs
from ..genome.transcript import Transcript
from ..variant import Consequence, Variant, VariantKind

__all__ = ["enumerate_exon_cnvs", "exon_interval_count"]


def exon_interval_count(transcript: Transcript) -> int:
    """``n(n+1)/2`` contiguous intervals over ``n`` exons, per direction."""
    n = len(transcript.exons)
    return n * (n + 1) // 2


def _spanned_cds(transcript: Transcript, first: int, last: int) -> tuple[int, int] | None:
    """CDS positions covered by exons ``first..last`` (0-based, transcript order)."""
    bounds = transcript._exon_tx_bounds
    tx_start, tx_end = bounds[first][0], bounds[last][1]
    cds_lo = max(tx_start, transcript.cds_start_tx) - transcript.cds_start_tx + 1
    cds_hi = min(tx_end, transcript.cds_end_tx) - transcript.cds_start_tx + 1
    if cds_hi < cds_lo:
        return None
    return cds_lo, cds_hi


def enumerate_exon_cnvs(transcript: Transcript) -> Iterator[Variant]:
    """Yield every contiguous exon interval as a deletion and as a duplication."""
    exons = transcript.exons
    for first in range(len(exons)):
        for last in range(first, len(exons)):
            spanned = _spanned_cds(transcript, first, last)
            coding_nt = 0 if spanned is None else spanned[1] - spanned[0] + 1
            genomic_lo = min(exons[first].start, exons[last].start)
            genomic_hi = max(exons[first].end, exons[last].end)
            for kind, consequence in (
                ("deletion", Consequence.EXON_DELETION),
                ("duplication", Consequence.EXON_DUPLICATION),
            ):
                yield Variant(
                    gene=transcript.gene,
                    transcript=transcript.transcript_id,
                    kind=VariantKind.CNV,
                    hgvs_c=hgvs.format_exon_cnv(
                        transcript.gene, exons[first].label, exons[last].label, kind
                    ),
                    hgvs_p=None,
                    grch38_pos=f"{transcript.chrom}-{genomic_lo}-{genomic_hi}-{kind}",
                    consequence=consequence,
                    consequence_terms=(consequence,),
                    mutational_distance=0,
                    codon_index=None,
                    cds_position=None if spanned is None else spanned[0],
                    attributes={
                        "exon_first": exons[first].label,
                        "exon_last": exons[last].label,
                        "exon_count": str(last - first + 1),
                        "coding_nt": str(coding_nt),
                        "in_frame": str(coding_nt % 3 == 0).lower(),
                        "spans_whole_gene": str(first == 0 and last == len(exons) - 1).lower(),
                    },
                )
