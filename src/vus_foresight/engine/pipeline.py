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
from .signature import SignatureCache, evidence_signature
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
    #: Evaluate once per distinct evidence profile instead of once per variant
    #: (spec section 5). Turning it off is not a fallback -- it is the reference
    #: implementation the cache is tested against.
    reuse_by_signature: bool = True
    #: Record every context path an evaluation reads, so a test can prove the
    #: specification's declared footprint is complete. Off by default: it costs
    #: a set insertion per lookup and is only meaningful under test.
    audit_context_reads: bool = False
    cache: SignatureCache = field(default_factory=SignatureCache)
    _extra_sources: dict[str, str] = field(default_factory=dict)
    #: Union of every context path read, when auditing is on.
    audited_paths: set[str] = field(default_factory=set)

    def build_context(self, variant: Variant) -> EvidenceContext:
        """Merge every adapter's view, then derive PVS1 from the result.

        PVS1 is computed last because its splice branch reads a prediction that
        an adapter supplies; deriving it first would make the tree blind to the
        very data whose absence it is supposed to report.
        """
        context = EvidenceContext(audit=self.audited_paths if self.audit_context_reads else None)
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

    def _assess(self, variant: Variant, context: EvidenceContext) -> tuple[Evaluation, GapAnalysis]:
        """Evaluate and analyse, reusing the result across identical profiles.

        The reuse is exact rather than approximate: the signature covers every
        context path the specification is able to read, so a hit means the
        engine had no input on which the two variants differ.
        """
        if not self.reuse_by_signature:
            evaluation = evaluate_variant(variant, context, self.spec)
            return evaluation, analyse_gap(evaluation, self.spec, context)

        signature = evidence_signature(variant, context, self.spec.field_footprint)
        cached = self.cache.get(signature)
        if cached is not None:
            return cached
        evaluation = evaluate_variant(variant, context, self.spec)
        entry = (evaluation, analyse_gap(evaluation, self.spec, context))
        self.cache.put(signature, entry)
        return entry

    def evaluate(self, variant: Variant) -> VariantResult:
        context = self.build_context(variant)
        evaluation, gap = self._assess(variant, context)
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


def default_registry(gene: GeneConfig, *, extra: Iterable[Adapter] = ()) -> AdapterRegistry:
    """The three computed adapters, plus whatever file-backed ones are supplied."""
    from ..adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter

    registry = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(gene)])
    for adapter in extra:
        registry.add(adapter)
    return registry
