"""The PVS1 decision tree, parameterised by the specification.

Follows the ClinGen SVI framework (Abou Tayoun et al. 2018) as implemented by
autoPVS1: where the premature termination codon lands relative to NMD, how much
of the protein is removed, and whether what is removed is critical.

Nothing here is BRCA-specific. The three gene-level facts it needs -- whether
loss of function is an established mechanism, where the functional regions are,
and which of them are critical -- all arrive from the gene config, and the
strengths each branch grants arrive from the spec.

One simplification is stated rather than hidden: for a canonical splice-site
variant the tree does not model which exon is skipped or whether skipping
preserves frame. Without transcript-level splicing evidence that determination
is not available, so the criterion asks for a splicing prediction and reports
itself unevaluable without one. That turns an un-ingested SpliceAI table into a
visible, actionable gap instead of an invisible assumption.
"""

from __future__ import annotations

from typing import Any

from ..genome.reference import GeneConfig
from ..genome.transcript import Transcript
from ..variant import Consequence, Variant
from .context import MISSING, EvidenceContext
from .spec import PVS1Config

__all__ = ["PVS1_NAMESPACE", "compute_pvs1"]

PVS1_NAMESPACE = "pvs1"

_TRUNCATING_WITH_PTC = (Consequence.NONSENSE, Consequence.FRAMESHIFT)
_CANONICAL_SPLICE = (Consequence.SPLICE_DONOR, Consequence.SPLICE_ACCEPTOR)


def _ptc_codon(variant: Variant) -> int | None:
    if variant.ptc_codon is not None:
        return variant.ptc_codon
    if variant.consequence is Consequence.NONSENSE:
        return variant.codon_index
    return None


def _removes_critical_region(
    gene: GeneConfig, config: PVS1Config, first_removed_residue: int, protein_length: int
) -> tuple[bool, tuple[str, ...]]:
    critical = set(config.critical_region_tags)
    hits = [
        region.name
        for region in gene.functional_regions
        if critical.intersection(region.tags)
        and region.end_aa >= first_removed_residue
        and region.start_aa <= protein_length
    ]
    return bool(hits), tuple(sorted(hits))


def compute_pvs1(
    variant: Variant,
    transcript: Transcript,
    gene: GeneConfig,
    config: PVS1Config,
    context: EvidenceContext,
) -> dict[str, Any]:
    """Evaluate the tree and return the ``pvs1`` namespace for the context.

    ``strength`` is absent from the returned mapping whenever the tree cannot
    reach a verdict. A criterion declaring ``requires: [pvs1.strength]`` then
    reports NOT_EVALUABLE with the reason attached, rather than quietly not
    firing.
    """
    result: dict[str, Any] = {"applicable": False}
    if not config.enabled:
        result["undetermined_reason"] = "pvs1_disabled_in_spec"
        return result

    terms = set(variant.consequence_terms) or {variant.consequence}

    if variant.consequence is Consequence.START_LOST:
        result.update(
            applicable=True,
            strength=config.start_lost_strength.value,
            nmd_escape=False,
            fraction_removed=1.0,
            rationale="initiation codon lost; downstream start use is not predictable",
        )
        return result

    if terms & set(_CANONICAL_SPLICE):
        prediction = context.get(config.splice_prediction_field)
        if config.splice_requires_prediction and prediction is MISSING:
            result.update(
                applicable=False,
                undetermined_reason="splice_prediction_missing",
                missing_field=config.splice_prediction_field,
                rationale=(
                    "canonical splice site; a splicing prediction is required before "
                    "PVS1 can be assigned"
                ),
            )
            return result
        disrupted = (
            prediction is not MISSING and prediction >= config.splice_prediction_threshold
        )
        if not disrupted:
            result.update(
                applicable=False,
                undetermined_reason="splice_not_predicted_disrupted",
                rationale=(
                    f"{config.splice_prediction_field}={prediction} below "
                    f"{config.splice_prediction_threshold}"
                ),
            )
            return result
        result.update(
            applicable=True,
            strength=config.splice_disrupted_strength.value,
            nmd_escape=False,
            rationale=(
                f"canonical splice site predicted disrupted "
                f"({config.splice_prediction_field}={prediction})"
            ),
        )
        return result

    if variant.consequence not in _TRUNCATING_WITH_PTC:
        result["undetermined_reason"] = "consequence_not_truncating"
        return result

    ptc = _ptc_codon(variant)
    if ptc is None:
        result["undetermined_reason"] = "ptc_position_unknown"
        return result

    protein_length = transcript.protein_length
    escapes = transcript.ptc_escapes_nmd(ptc)
    fraction_removed = (protein_length - ptc + 1) / protein_length
    critical, region_names = _removes_critical_region(
        gene, config, ptc, protein_length
    )

    result.update(
        applicable=True,
        nmd_escape=escapes,
        fraction_removed=round(fraction_removed, 6),
        removes_critical_region=critical,
        critical_regions_removed=list(region_names),
        ptc_codon=ptc,
    )

    if not escapes:
        established = gene.lof_mechanism == "established"
        strength = (
            config.nmd_triggered_strength
            if established
            else config.nmd_triggered_unestablished_strength
        )
        result["strength"] = strength.value
        result["rationale"] = (
            f"PTC at codon {ptc} is NMD-competent; loss of function is "
            f"{gene.lof_mechanism} for {gene.gene}"
        )
        return result

    if critical:
        result["strength"] = config.nmd_escape_critical_strength.value
        result["rationale"] = (
            f"PTC at codon {ptc} escapes NMD but removes critical region(s) "
            f"{', '.join(region_names)}"
        )
    elif fraction_removed >= config.min_fraction_removed:
        result["strength"] = config.nmd_escape_large_strength.value
        result["rationale"] = (
            f"PTC at codon {ptc} escapes NMD and removes "
            f"{fraction_removed:.1%} of the protein"
        )
    else:
        result["strength"] = config.nmd_escape_small_strength.value
        result["rationale"] = (
            f"PTC at codon {ptc} escapes NMD and removes only "
            f"{fraction_removed:.1%} of the protein"
        )
    return result
