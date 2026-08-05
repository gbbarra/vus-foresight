"""Adapters that compute their evidence rather than looking it up.

These three need no external dataset: everything they publish is a function of
the variant, the transcript geometry, or the gene config. They are still
adapters, so that a specification addresses ``variant.consequence`` exactly the
way it addresses ``frequency.gnomad.popmax_af`` and the engine keeps no
privileged channel for its own facts.
"""

from __future__ import annotations

from typing import Any

from ..genome.reference import GeneConfig
from ..genome.transcript import Transcript
from ..variant import Variant
from .base import Adapter

__all__ = ["RegionAdapter", "TranscriptAdapter", "VariantAdapter"]


class VariantAdapter(Adapter):
    """The variant's own annotation, addressable from a rule."""

    namespace = "variant"
    name = "annotation"

    def __init__(self, version: str = "1") -> None:
        super().__init__(version)

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        return {
            "gene": variant.gene,
            "kind": variant.kind.value,
            "consequence": variant.consequence.value,
            "consequence_terms": [t.value for t in variant.consequence_terms],
            "mutational_distance": variant.mutational_distance,
            "codon_index": variant.codon_index,
            "cds_position": variant.cds_position,
            "ref_aa": variant.ref_aa,
            "alt_aa": variant.alt_aa,
            "ptc_codon": variant.ptc_codon,
            "hgvs_c": variant.hgvs_c,
            "hgvs_p": variant.hgvs_p,
        }


class TranscriptAdapter(Adapter):
    """Transcript geometry facts a specification may want to gate on."""

    namespace = "transcript"
    name = "mane"

    def __init__(self, version: str = "1") -> None:
        super().__init__(version)

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "id": transcript.transcript_id,
            "strand": transcript.strand,
            "exon_count": len(transcript.exons),
            "protein_length": transcript.protein_length,
            "cds_length": transcript.cds_length,
        }
        if variant.cds_position is not None:
            tx_pos = transcript.cds_position_to_tx(variant.cds_position)
            exon_index = transcript.exon_index_at_tx(tx_pos)
            payload["exon_label"] = transcript.exons[exon_index].label
            payload["is_last_exon"] = exon_index == len(transcript.exons) - 1
            payload["is_penultimate_exon"] = exon_index == len(transcript.exons) - 2
        if variant.ptc_codon is not None:
            payload["ptc_escapes_nmd"] = transcript.ptc_escapes_nmd(variant.ptc_codon)
        return payload


class RegionAdapter(Adapter):
    """Functional regions of the protein that contain the affected residue.

    Feeds PM1 and the PVS1 critical-region branch. The regions themselves come
    from the gene config, so adding a gene adds a YAML block, not code.
    """

    namespace = "region"
    name = "gene_config"

    def __init__(self, gene: GeneConfig, version: str = "1") -> None:
        super().__init__(version)
        self._gene = gene

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        residue = variant.codon_index
        if residue is None:
            return {"names": [], "tags": []}
        hits = [r for r in self._gene.functional_regions if r.contains(residue)]
        tags = sorted({tag for r in hits for tag in r.tags})
        return {
            "names": sorted(r.name for r in hits),
            "tags": tags,
            "count": len(hits),
        }
