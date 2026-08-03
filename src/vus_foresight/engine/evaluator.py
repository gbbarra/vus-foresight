"""Criterion evaluation, producing the complete trace.

Principle 1 of the spec: every evaluation returns every criterion that was
tested, applied or not, with the reason. The map is derived from the trace --
without the trace there is no map. So this module's return value is not a
classification with some diagnostics attached; it is a trace from which a
classification happens to be computable.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..acmg import (
    ACMGClass,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    SkipReason,
    Strength,
)
from ..gapmap import AppliedCriterion, SkippedCriterion
from ..variant import Variant
from .context import MISSING, EvidenceContext, LookupLog
from .predicates import evaluate_rule
from .spec import CriterionSpec, VCEPSpec

__all__ = ["Evaluation", "evaluate_variant", "render_evidence"]

_TEMPLATE_TOKEN = re.compile(r"\{([A-Za-z0-9_.]+)\}")


def render_evidence(template: str | None, context: EvidenceContext, fallback: str) -> str:
    """Fill ``{dotted.path}`` tokens from the context.

    An unresolvable token renders as ``?`` rather than raising: the evidence
    string is documentation, and a template typo must not be able to abort a
    map-wide run.
    """
    if template is None:
        return fallback

    def substitute(match: re.Match[str]) -> str:
        value = context.get(match.group(1))
        return "?" if value is MISSING else str(value)

    return _TEMPLATE_TOKEN.sub(substitute, template)


@dataclass(slots=True)
class Evaluation:
    """The full trace for one variant under one specification."""

    applied: list[AppliedCriterion] = field(default_factory=list)
    skipped: list[SkippedCriterion] = field(default_factory=list)
    points: int = 0
    acmg_class: ACMGClass = ACMGClass.UNCERTAIN
    #: code -> the paths that were absent, for criteria stopped by missing data.
    missing_by_code: dict[str, tuple[str, ...]] = field(default_factory=dict)

    @property
    def applied_codes(self) -> tuple[str, ...]:
        return tuple(c.code for c in self.applied)

    def outcome_of(self, code: str) -> CriterionOutcome | None:
        for applied in self.applied:
            if applied.code == code:
                return CriterionOutcome.APPLIED
        for skipped in self.skipped:
            if skipped.code == code:
                return skipped.outcome
        return None

    def has_directional_conflict(self, min_points_each: int = 1) -> bool:
        """Whether both directions contribute at least ``min_points_each``.

        The threshold matters: a lone supporting criterion on each side is
        ordinary curation, not a variant blocked by contradictory evidence.
        """
        pathogenic = sum(
            c.points for c in self.applied if c.direction is Direction.PATHOGENIC
        )
        benign = -sum(c.points for c in self.applied if c.direction is Direction.BENIGN)
        return pathogenic >= min_points_each and benign >= min_points_each


def _consequence_gate(criterion: CriterionSpec, variant: Variant) -> bool:
    if not criterion.applies_to:
        return True
    terms = set(variant.consequence_terms) or {variant.consequence}
    return bool(terms.intersection(criterion.applies_to))


def _resolve_strength(
    criterion: CriterionSpec, context: EvidenceContext, log: LookupLog
) -> Strength:
    for rung in criterion.strength_ladder:
        if rung.when.evaluate(context, log):
            return rung.strength
    return criterion.strength


def _evaluate_one(
    criterion: CriterionSpec,
    variant: Variant,
    context: EvidenceContext,
    spec: VCEPSpec,
) -> AppliedCriterion | SkippedCriterion:
    if not _consequence_gate(criterion, variant):
        return SkippedCriterion(
            code=criterion.code,
            direction=criterion.direction,
            evidence_class=criterion.evidence_class,
            outcome=CriterionOutcome.NOT_APPLICABLE,
            reason=SkipReason.CONSEQUENCE_OUT_OF_SCOPE,
            detail=(
                f"applies to {', '.join(c.value for c in criterion.applies_to)}; "
                f"variant is {variant.consequence.value}"
            ),
        )

    if criterion.evidence_class is EvidenceClass.EXTRINSIC:
        return SkippedCriterion(
            code=criterion.code,
            direction=criterion.direction,
            evidence_class=criterion.evidence_class,
            outcome=CriterionOutcome.NOT_APPLICABLE,
            reason=SkipReason.EXTRINSIC_EVIDENCE_REQUIRED,
            detail=criterion.description,
        )

    prerequisite_log = LookupLog()
    for path in criterion.requires:
        context.get(path, prerequisite_log)
    if prerequisite_log.has_missing:
        return SkippedCriterion(
            code=criterion.code,
            direction=criterion.direction,
            evidence_class=criterion.evidence_class,
            outcome=CriterionOutcome.NOT_EVALUABLE,
            reason=SkipReason.DATA_MISSING,
            missing_fields=tuple(prerequisite_log.missing),
            detail=_undetermined_detail(context, criterion),
        )

    fired, log = evaluate_rule(criterion.rule, context)
    if not fired and log.has_missing:
        # A rule that could not see its inputs did not "fail" -- it was never
        # answerable. Conflating the two would understate the gap.
        return SkippedCriterion(
            code=criterion.code,
            direction=criterion.direction,
            evidence_class=criterion.evidence_class,
            outcome=CriterionOutcome.NOT_EVALUABLE,
            reason=SkipReason.DATA_MISSING,
            missing_fields=tuple(log.missing),
            detail=_undetermined_detail(context, criterion),
        )
    if not fired:
        return SkippedCriterion(
            code=criterion.code,
            direction=criterion.direction,
            evidence_class=criterion.evidence_class,
            outcome=CriterionOutcome.NOT_MET,
            reason=SkipReason.RULE_NOT_MET,
            detail=None,
        )

    strength = _resolve_strength(criterion, context, log)
    points = spec.points_for(criterion, strength)
    namespaces = {path.split(".", 1)[0] for path in log.touched} or {"engine"}
    source = ", ".join(
        sorted({context.sources.get(ns, ns) for ns in namespaces})
    )
    return AppliedCriterion(
        code=criterion.code,
        direction=criterion.direction,
        strength=strength,
        points=points,
        evidence_class=criterion.evidence_class,
        evidence=render_evidence(
            criterion.evidence_template,
            context,
            fallback=criterion.description or f"{criterion.code} rule satisfied",
        ),
        source=source,
    )


def _undetermined_detail(context: EvidenceContext, criterion: CriterionSpec) -> str | None:
    """Surface a structured reason when one is available (currently PVS1's)."""
    if criterion.code.startswith("PVS1"):
        reason = context.get("pvs1.undetermined_reason")
        if reason is not MISSING:
            return str(reason)
    return None


def _resolve_mutex(
    applied: list[AppliedCriterion],
    skipped: list[SkippedCriterion],
    spec: VCEPSpec,
) -> None:
    """Keep only the strongest applied member of each mutually exclusive group."""
    groups: dict[str, list[AppliedCriterion]] = {}
    for entry in applied:
        group = spec.by_code(entry.code).mutex_group
        if group:
            groups.setdefault(group, []).append(entry)

    for group, members in groups.items():
        if len(members) < 2:
            continue
        winner = max(members, key=lambda m: (abs(m.points), m.code))
        for member in members:
            if member is winner:
                continue
            applied.remove(member)
            skipped.append(
                SkippedCriterion(
                    code=member.code,
                    direction=member.direction,
                    evidence_class=member.evidence_class,
                    outcome=CriterionOutcome.NOT_MET,
                    reason=SkipReason.SUPERSEDED_BY_MUTEX,
                    detail=f"superseded by {winner.code} in mutex group {group!r}",
                )
            )


def evaluate_variant(
    variant: Variant, context: EvidenceContext, spec: VCEPSpec
) -> Evaluation:
    """Test every criterion in the specification against one variant."""
    applied: list[AppliedCriterion] = []
    skipped: list[SkippedCriterion] = []
    missing_by_code: dict[str, tuple[str, ...]] = {}

    for criterion in spec.criteria:
        outcome = _evaluate_one(criterion, variant, context, spec)
        if isinstance(outcome, AppliedCriterion):
            applied.append(outcome)
        else:
            skipped.append(outcome)
            if outcome.missing_fields:
                missing_by_code[criterion.code] = outcome.missing_fields

    _resolve_mutex(applied, skipped, spec)

    applied.sort(key=lambda c: c.code)
    skipped.sort(key=lambda c: c.code)
    points = sum(c.points for c in applied)
    return Evaluation(
        applied=applied,
        skipped=skipped,
        points=points,
        acmg_class=spec.point_system.classify(points),
        missing_by_code=missing_by_code,
    )
