"""Assemble a gap map row from a variant, the adapters, and a specification.

The pipeline is the only place the pieces meet, and it is deliberately thin:
build the context, evaluate, analyse the gap, stamp provenance. Anything that
needed to know about a gene was decided before this point, in configuration.

Determinism (spec section 11) is a property of this function. ``computed_at`` is
injected rather than read from the clock, adapters are ordered, criteria are
sorted, and nothing iterates a set.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date, datetime
from typing import Iterable, Iterator

from ..adapters.base import Adapter, AdapterRegistry
from ..gapmap import GapMapRow
from ..genome.reference import GeneConfig
from ..genome.transcript import Transcript
from ..variant import Variant, VariantKind
from .context import EvidenceContext
from .equivalence import equivalence_class_id
from .evaluator import Evaluation, evaluate_variant
from .gap import GapAnalysis, analyse_gap
from .pvs1 import PVS1_NAMESPACE, compute_pvs1
from .spec import VCEPSpec

__all__ = ["MapRunner", "VariantResult"]


@dataclass(frozen=True, slots=True)
class VariantResult:
    """A row plus the intermediate objects, so callers can report on the run."""

    variant: Variant
    context: EvidenceContext
    evaluation: Evaluation
    gap: GapAnalysis
    row: GapMapRow


@dataclass(slots=True)
class MapRunner:
    """Runs one gene against one specification with one set of adapters."""

    transcript: Transcript
    gene: GeneConfig
    spec: VCEPSpec
    adapters: AdapterRegistry
    computed_at: datetime
    clinvar_snapshot: date | None = None
    gnomad_version: str | None = None
    _extra_sources: dict[str, str] = field(default_factory=dict)

    def build_context(self, variant: Variant) -> EvidenceContext:
        """Merge every adapter's view, then derive PVS1 from the result.

        PVS1 is computed last because its splice branch reads a prediction that
        an adapter supplies; deriving it first would make the tree blind to the
        very data whose absence it is supposed to report.
        """
        context = EvidenceContext()
        for adapter in self.adapters:
            context.merge(
                adapter.namespace, adapter.lookup(variant, self.transcript), adapter.source_id
            )
        context.merge(
            PVS1_NAMESPACE,
            compute_pvs1(variant, self.transcript, self.gene, self.spec.pvs1, context),
            f"pvs1@{self.spec.spec_version}",
        )
        return context

    def evaluate(self, variant: Variant) -> VariantResult:
        context = self.build_context(variant)
        evaluation = evaluate_variant(variant, context, self.spec)
        gap = analyse_gap(evaluation, self.spec, context)
        row = GapMapRow(
            gene=variant.gene,
            transcript=variant.transcript,
            hgvs_c=variant.hgvs_c,
            hgvs_p=variant.hgvs_p,
            grch38_pos=variant.grch38_pos,
            consequence=variant.consequence,
            variant_kind=variant.kind,
            equivalence_class_id=equivalence_class_id(variant, evaluation, self.spec),
            mutational_distance=variant.mutational_distance,
            criteria_applied=tuple(evaluation.applied),
            criteria_evaluated_not_applied=tuple(evaluation.skipped),
            points_current=evaluation.points,
            class_current=evaluation.acmg_class,
            points_ceiling_intrinsic=gap.points_ceiling_intrinsic,
            class_ceiling_intrinsic=gap.class_ceiling_intrinsic,
            gap_to_LP=gap.gap_to_lp,
            gap_to_LB=gap.gap_to_lb,
            minimum_sufficient_sets=gap.minimum_sufficient_sets,
            blocking_reason=gap.blocking_reason,
            spec_version=self.spec.spec_version,
            clinvar_snapshot=self.clinvar_snapshot,
            gnomad_version=self.gnomad_version,
            source_versions=dict(sorted({**self.adapters.versions, **self._extra_sources}.items())),
            computed_at=self.computed_at,
        )
        return VariantResult(
            variant=variant, context=context, evaluation=evaluation, gap=gap, row=row
        )

    def run(self, variants: Iterable[Variant]) -> Iterator[VariantResult]:
        """Evaluate a stream of variants, preserving the enumerator's order."""
        for variant in variants:
            if variant.kind is VariantKind.CNV:
                # Copy number changes are scored by a different framework and
                # emitted by a different module; letting them through here would
                # mix two point scales in one column.
                continue
            yield self.evaluate(variant)

    def declare_source(self, namespace: str, source_id: str) -> None:
        """Record an extra provenance entry that is not a context namespace."""
        self._extra_sources[namespace] = source_id


def default_registry(
    gene: GeneConfig, *, extra: Iterable[Adapter] = ()
) -> AdapterRegistry:
    """The three computed adapters, plus whatever file-backed ones are supplied."""
    from ..adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter

    registry = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(gene)])
    for adapter in extra:
        registry.add(adapter)
    return registry
