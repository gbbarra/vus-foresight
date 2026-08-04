"""The gap: ceiling, distance to a verdict, and what would close it.

This is the product (spec section 7). The point model makes it tractable:
combining criteria is addition, so "what is still missing" is subtraction and
"what would be enough" is a subset-sum over a few dozen items. Exhaustive search
with pruning is correct and fast at this size; a solver would add a dependency
and take away the ability to explain the answer.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from ..acmg import (
    ACMGClass,
    BlockingReason,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FEASIBILITY_COST,
    FeasibilityTag,
    SkipReason,
)
from ..gapmap import EvidenceRequirement, EvidenceSet
from .context import EvidenceContext, LookupLog
from .evaluator import Evaluation
from .spec import CriterionSpec, VCEPSpec

__all__ = [
    "Candidate",
    "GapAnalysis",
    "candidate_requirements",
    "intrinsic_ceiling",
    "minimum_sufficient_sets",
    "derive_blocking_reason",
    "analyse_gap",
]


@dataclass(frozen=True, slots=True)
class Candidate:
    """A criterion that is not applied but could be, and at what cost."""

    criterion: CriterionSpec
    requirement: EvidenceRequirement

    @property
    def code(self) -> str:
        return self.criterion.code

    @property
    def points(self) -> int:
        return self.requirement.points

    @property
    def mutex_group(self) -> str | None:
        return self.criterion.mutex_group


def _resolve_feasibility(criterion: CriterionSpec, context: EvidenceContext) -> FeasibilityTag:
    log = LookupLog()
    for rung in criterion.feasibility_ladder:
        if rung.when.evaluate(context, log):
            return rung.feasibility
    if criterion.feasibility is not None:
        return criterion.feasibility
    # Intrinsic and semi-intrinsic evidence needs no patient by definition, so
    # the honest default is "somebody already has this, we have not loaded it".
    return FeasibilityTag.AVAILABLE_UNINGESTED


def candidate_requirements(
    evaluation: Evaluation, spec: VCEPSpec, context: EvidenceContext
) -> list[Candidate]:
    """Criteria that could still contribute points, with feasibility resolved.

    A criterion whose rule was evaluated and *not met* is excluded: its data was
    present and the answer was no. Only criteria blocked by missing data, and
    extrinsic criteria this system structurally cannot evaluate, are candidates.
    """
    applied_codes = set(evaluation.applied_codes)
    occupied_groups = {
        spec.by_code(code).mutex_group for code in applied_codes if spec.by_code(code).mutex_group
    }

    candidates: list[Candidate] = []
    for skipped in evaluation.skipped:
        if skipped.reason not in (
            SkipReason.DATA_MISSING,
            SkipReason.EXTRINSIC_EVIDENCE_REQUIRED,
        ):
            continue
        criterion = spec.by_code(skipped.code)
        if criterion.mutex_group and criterion.mutex_group in occupied_groups:
            continue
        strength = criterion.ceiling_strength
        candidates.append(
            Candidate(
                criterion=criterion,
                requirement=EvidenceRequirement(
                    code=criterion.code,
                    direction=criterion.direction,
                    strength=strength,
                    points=spec.points_for(criterion, strength),
                    evidence_class=criterion.evidence_class,
                    feasibility=_resolve_feasibility(criterion, context),
                    required_observations=criterion.observations_for(strength),
                    description=criterion.description,
                ),
            )
        )
    candidates.sort(key=lambda c: c.code)
    return candidates


def intrinsic_ceiling(
    evaluation: Evaluation, candidates: list[Candidate], spec: VCEPSpec
) -> tuple[int, ACMGClass]:
    """Highest point total reachable without ever observing a patient.

    "Maximum" is meant in the pathogenic direction, matching the field's name in
    the output schema. Only intrinsic and semi-intrinsic candidates count, one
    per mutually exclusive group, and every one of them is obtainable by
    definition -- that is what makes them intrinsic.

    The benign side of the same question is answered by ``gap_to_LB`` together
    with the minimum sufficient sets that target LB.
    """
    best_by_group: dict[str, int] = {}
    ungrouped = 0
    for candidate in candidates:
        if candidate.criterion.evidence_class is EvidenceClass.EXTRINSIC:
            continue
        if candidate.requirement.direction is not Direction.PATHOGENIC:
            continue
        group = candidate.mutex_group
        if group:
            best_by_group[group] = max(best_by_group.get(group, 0), candidate.points)
        else:
            ungrouped += candidate.points
    ceiling = evaluation.points + ungrouped + sum(best_by_group.values())
    return ceiling, spec.point_system.classify(ceiling)


def _acquisition_cost(requirements: tuple[EvidenceRequirement, ...], excess: int) -> int:
    """Deterministic ranking key: feasibility dominates, then size, then slack."""
    feasibility = sum(FEASIBILITY_COST[r.feasibility] for r in requirements)
    return feasibility * 1000 + len(requirements) * 10 + max(0, excess)


def minimum_sufficient_sets(
    candidates: list[Candidate],
    gap: int,
    target: ACMGClass,
    direction: Direction,
    spec: VCEPSpec,
) -> tuple[EvidenceSet, ...]:
    """Minimal subsets whose combined points close ``gap``, cheapest first.

    ``gap`` is a positive magnitude. Minimality is by set inclusion, not by
    cardinality: a two-criterion set is still reported when no single criterion
    suffices, but a set containing an already-sufficient subset never is.
    """
    if gap <= 0:
        return ()

    usable = [
        c
        for c in candidates
        if c.requirement.direction is direction
        and (
            spec.gap.include_intractable
            or c.requirement.feasibility is not FeasibilityTag.INTRACTABLE
        )
    ]
    # Points are signed; work in magnitudes so both directions share the search.
    usable = [c for c in usable if abs(c.points) > 0]
    usable.sort(key=lambda c: (-abs(c.points), c.code))

    accepted: list[tuple[frozenset[str], EvidenceSet]] = []
    for size in range(1, min(spec.gap.max_cardinality, len(usable)) + 1):
        for combo in combinations(usable, size):
            groups = [c.mutex_group for c in combo if c.mutex_group]
            if len(groups) != len(set(groups)):
                continue
            total = sum(abs(c.points) for c in combo)
            if total < gap:
                continue
            codes = frozenset(c.code for c in combo)
            if any(previous <= codes for previous, _ in accepted):
                continue
            requirements = tuple(sorted((c.requirement for c in combo), key=lambda r: r.code))
            worst = max(requirements, key=lambda r: FEASIBILITY_COST[r.feasibility])
            accepted.append(
                (
                    codes,
                    EvidenceSet(
                        target=target,
                        requirements=requirements,
                        total_points=sum(r.points for r in requirements),
                        feasibility=worst.feasibility,
                        acquisition_cost=_acquisition_cost(requirements, total - gap),
                    ),
                )
            )

    sets = [entry for _, entry in accepted]
    sets.sort(key=lambda s: (s.acquisition_cost, s.codes))
    return tuple(sets[: spec.gap.max_sets_per_target])


def derive_blocking_reason(
    evaluation: Evaluation,
    spec: VCEPSpec,
    candidates: list[Candidate],
    sets: tuple[EvidenceSet, ...],
) -> BlockingReason:
    """Categorise *why* this variant is stuck, using only the spec's own mapping.

    The enum exists so the map aggregates: "how many BRCA2 VUS are blocked for
    want of functional data" has to be one ``GROUP BY``. That only works if the
    assignment is deterministic and does not depend on criterion codes being
    special-cased in Python -- so every criterion declares its own
    ``blocks_as``, and this function only ranks.
    """
    if evaluation.acmg_class is not ACMGClass.UNCERTAIN:
        return BlockingReason.RESOLVED_NOT_BLOCKED
    if evaluation.has_directional_conflict(spec.blocking.conflict_min_points_each):
        return spec.blocking.conflict_reason

    signals: set[BlockingReason] = set()

    for skipped in evaluation.skipped:
        if skipped.outcome is CriterionOutcome.NOT_EVALUABLE:
            reason = spec.by_code(skipped.code).blocks_as
            if reason is not None:
                signals.add(reason)

    cheapest = sets[0] if sets else None
    if cheapest is not None:
        for requirement in cheapest.requirements:
            reason = spec.by_code(requirement.code).blocks_as
            if reason is not None:
                signals.add(reason)

    if not signals:
        for skipped in evaluation.skipped:
            if skipped.reason is SkipReason.RULE_NOT_MET:
                reason = spec.by_code(skipped.code).blocks_as
                if reason is not None:
                    signals.add(reason)

    for reason in spec.blocking.priority:
        if reason in signals:
            return reason
    return spec.blocking.default


@dataclass(frozen=True, slots=True)
class GapAnalysis:
    """Everything section 7 asks for, for one variant."""

    gap_to_lp: int | None
    gap_to_lb: int | None
    points_ceiling_intrinsic: int
    class_ceiling_intrinsic: ACMGClass
    minimum_sufficient_sets: tuple[EvidenceSet, ...]
    blocking_reason: BlockingReason
    candidates: tuple[Candidate, ...]

    @property
    def available_uningested(self) -> tuple[str, ...]:
        """Codes whose evidence exists publicly and simply is not loaded.

        The most actionable output the system produces: engineering work, not
        bench work.
        """
        return tuple(
            sorted(
                c.code
                for c in self.candidates
                if c.requirement.feasibility is FeasibilityTag.AVAILABLE_UNINGESTED
            )
        )


def analyse_gap(evaluation: Evaluation, spec: VCEPSpec, context: EvidenceContext) -> GapAnalysis:
    """Run the whole of section 7 for one evaluated variant."""
    candidates = candidate_requirements(evaluation, spec, context)
    points = evaluation.points
    gap_lp = spec.point_system.gap_to_likely_pathogenic(points)
    gap_lb = spec.point_system.gap_to_likely_benign(points)

    sets: list[EvidenceSet] = []
    if gap_lp is not None:
        sets.extend(
            minimum_sufficient_sets(
                candidates, gap_lp, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
            )
        )
    if gap_lb is not None:
        sets.extend(
            minimum_sufficient_sets(
                candidates, gap_lb, ACMGClass.LIKELY_BENIGN, Direction.BENIGN, spec
            )
        )
    ordered = tuple(sorted(sets, key=lambda s: (s.acquisition_cost, s.target.value, s.codes)))

    ceiling_points, ceiling_class = intrinsic_ceiling(evaluation, candidates, spec)
    return GapAnalysis(
        gap_to_lp=gap_lp,
        gap_to_lb=gap_lb,
        points_ceiling_intrinsic=ceiling_points,
        class_ceiling_intrinsic=ceiling_class,
        minimum_sufficient_sets=ordered,
        blocking_reason=derive_blocking_reason(evaluation, spec, candidates, ordered),
        candidates=tuple(candidates),
    )
