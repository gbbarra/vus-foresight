"""Exon-level copy number scoring -- a separate framework, the same schema.

Spec section 3 is explicit: CNVs are scored by the ClinGen CNV framework (Riggs
et al. 2020), not by the Tavtigian point model, so they get their own module.
They are emitted into :class:`~vus_foresight.gapmap.GapMapRow` all the same, but
into their own output file, because two scoring scales sharing one
``points_current`` column would silently corrupt every aggregation over it.

The section codes below (2A, 2C-1, 2E, 3A) are the framework's own; the values
attached to them are configuration.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime

from pydantic import BaseModel, ConfigDict

from ..acmg import (
    ACMGClass,
    BlockingReason,
    Direction,
    EvidenceClass,
    Strength,
)
from ..gapmap import AppliedCriterion, GapMapRow
from ..genome.reference import GeneConfig
from ..genome.transcript import Transcript
from ..variant import Consequence, Variant, VariantKind

__all__ = ["CNVScore", "CNVScoringConfig", "cnv_row", "score_cnv"]

#: The framework works in hundredths; the schema's integer point column stores
#: score x 100 so that no floating point enters a reproducibility comparison.
SCALE = 100


class CNVScoringConfig(BaseModel):
    """ClinGen CNV section values, as configuration."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    framework_id: str = "clingen_cnv"
    framework_version: str = "1.0"

    threshold_pathogenic: int = 99
    threshold_likely_pathogenic: int = 90
    threshold_likely_benign: int = -90
    threshold_benign: int = -99

    #: Loss, section 2A: complete overlap of a gene with established HI.
    loss_full_gene_established: int = 100
    #: Loss, section 2A applied to a gene whose LoF mechanism is not established.
    loss_full_gene_unestablished: int = 15
    #: Loss, section 2E: both breakpoints inside the gene, reading frame lost.
    loss_intragenic_frame_disrupting: int = 90
    #: Loss of an in-frame stretch of exons.
    loss_intragenic_in_frame: int = 15
    #: Gain, section 3A: a whole-gene duplication is not evidence of LoF.
    gain_full_gene: int = 0
    #: Gain that disrupts the reading frame within the gene.
    gain_intragenic_frame_disrupting: int = 90
    gain_intragenic_in_frame: int = 0

    def classify(self, score: int) -> ACMGClass:
        if score >= self.threshold_pathogenic:
            return ACMGClass.PATHOGENIC
        if score >= self.threshold_likely_pathogenic:
            return ACMGClass.LIKELY_PATHOGENIC
        if score <= self.threshold_benign:
            return ACMGClass.BENIGN
        if score <= self.threshold_likely_benign:
            return ACMGClass.LIKELY_BENIGN
        return ACMGClass.UNCERTAIN

    @property
    def version_string(self) -> str:
        return f"{self.framework_id}@{self.framework_version}"


@dataclass(frozen=True, slots=True)
class CNVScore:
    score: int
    section: str
    rationale: str

    @property
    def as_fraction(self) -> float:
        return self.score / SCALE


def score_cnv(variant: Variant, gene: GeneConfig, config: CNVScoringConfig) -> CNVScore:
    """Score one exon-level deletion or duplication."""
    if variant.kind is not VariantKind.CNV:
        raise ValueError(f"{variant.variant_id} is not a copy number variant")

    whole_gene = variant.attributes.get("spans_whole_gene") == "true"
    in_frame = variant.attributes.get("in_frame") == "true"
    established = gene.lof_mechanism == "established"
    exon_count = variant.attributes.get("exon_count", "?")

    if variant.consequence is Consequence.EXON_DELETION:
        if whole_gene:
            if established:
                return CNVScore(
                    config.loss_full_gene_established,
                    "2A",
                    f"complete deletion of {gene.gene}, an established haploinsufficient gene",
                )
            return CNVScore(
                config.loss_full_gene_unestablished,
                "2A",
                f"complete deletion of {gene.gene}, whose loss-of-function "
                f"mechanism is {gene.lof_mechanism}",
            )
        if in_frame:
            return CNVScore(
                config.loss_intragenic_in_frame,
                "2E",
                f"intragenic deletion of {exon_count} exon(s) preserving the reading frame",
            )
        return CNVScore(
            config.loss_intragenic_frame_disrupting,
            "2E",
            f"intragenic deletion of {exon_count} exon(s) disrupting the reading frame",
        )

    if whole_gene:
        return CNVScore(
            config.gain_full_gene,
            "3A",
            f"whole-gene duplication of {gene.gene}; not evidence of loss of function",
        )
    if in_frame:
        return CNVScore(
            config.gain_intragenic_in_frame,
            "3A",
            f"intragenic duplication of {exon_count} exon(s) preserving the reading frame",
        )
    return CNVScore(
        config.gain_intragenic_frame_disrupting,
        "2E",
        f"intragenic duplication of {exon_count} exon(s) disrupting the reading frame",
    )


def cnv_row(
    variant: Variant,
    transcript: Transcript,
    gene: GeneConfig,
    config: CNVScoringConfig,
    *,
    computed_at: datetime,
    clinvar_snapshot: date | None = None,
    gnomad_version: str | None = None,
) -> GapMapRow:
    """Emit a scored CNV into the shared output schema.

    ``points_current`` is the ClinGen score times 100, and ``spec_version``
    names the CNV framework rather than the VCEP specification -- so a query
    that mixes the two files can always tell which scale it is looking at.
    """
    scored = score_cnv(variant, gene, config)
    acmg_class = config.classify(scored.score)
    applied = (
        AppliedCriterion(
            code=f"CNV_{scored.section}",
            direction=Direction.PATHOGENIC if scored.score >= 0 else Direction.BENIGN,
            strength=Strength.SUPPORTING,
            points=scored.score,
            evidence_class=EvidenceClass.INTRINSIC,
            evidence=scored.rationale,
            source=config.version_string,
        ),
    )
    gap_lp = (
        None
        if scored.score >= config.threshold_likely_pathogenic
        else (config.threshold_likely_pathogenic - scored.score)
    )
    gap_lb = (
        None
        if scored.score <= config.threshold_likely_benign
        else (scored.score - config.threshold_likely_benign)
    )
    return GapMapRow(
        gene=variant.gene,
        transcript=variant.transcript,
        hgvs_c=variant.hgvs_c,
        hgvs_p=None,
        grch38_pos=variant.grch38_pos,
        consequence=variant.consequence,
        variant_kind=variant.kind,
        equivalence_class_id=f"cnv_{variant.attributes.get('exon_first')}_"
        f"{variant.attributes.get('exon_last')}_{variant.consequence.value}",
        mutational_distance=0,
        criteria_applied=applied,
        criteria_evaluated_not_applied=(),
        points_current=scored.score,
        class_current=acmg_class,
        points_ceiling_intrinsic=scored.score,
        class_ceiling_intrinsic=acmg_class,
        gap_to_LP=gap_lp,
        gap_to_LB=gap_lb,
        minimum_sufficient_sets=(),
        blocking_reason=(
            BlockingReason.RESOLVED_NOT_BLOCKED
            if acmg_class is not ACMGClass.UNCERTAIN
            else BlockingReason.MISSING_CASE_CONTROL
        ),
        spec_version=config.version_string,
        clinvar_snapshot=clinvar_snapshot,
        gnomad_version=gnomad_version,
        source_versions={"cnv_framework": config.version_string},
        computed_at=computed_at,
    )
