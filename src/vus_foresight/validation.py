"""The falsifiable claim (spec section 10).

The map asserts something that can be wrong, and checking it is what turns a
taxonomy into science: compute the map with the evidence state of date T, take
the variants that were VUS at T and were reclassified in ClinVar by T+n, and ask
four questions.

1. **Recall of resolvability.** Among those that were resolved, what fraction did
   the map mark resolvable at all?
2. **Accuracy of cause.** Did the predicted ``blocking_reason`` match the kind of
   evidence that actually turned up in the submission record? This is the
   interesting one, and the publishable one.
3. **Direction.** Did the map point towards pathogenic when the variant went
   pathogenic, and benign when it went benign?
4. **Temporal calibration.** Were variants with a smaller gap resolved sooner?

This is a validation *study*, not a unit test. It runs on demand against
``vus-hindsight`` outcomes; nothing in CI depends on it.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence

from .acmg import ACMGClass, BlockingReason, FeasibilityTag
from .gapmap import GapMapRow

__all__ = [
    "Outcome",
    "ValidationResult",
    "read_outcomes",
    "is_resolvable",
    "predicted_direction",
    "spearman",
    "validate",
]

#: Which blocking reason each kind of newly-submitted evidence corresponds to.
#: A study can override it; the default reflects how ClinVar submissions
#: describe their evidence.
EVIDENCE_TO_BLOCKING: dict[str, BlockingReason] = {
    "functional": BlockingReason.MISSING_FUNCTIONAL,
    "segregation": BlockingReason.MISSING_SEGREGATION,
    "case_control": BlockingReason.MISSING_CASE_CONTROL,
    "population_frequency": BlockingReason.FREQUENCY_UNINFORMATIVE,
    "computational": BlockingReason.PREDICTOR_INDETERMINATE,
}


@dataclass(frozen=True, slots=True)
class Outcome:
    """What actually happened to one variant between T and T+n."""

    variant_id: str
    class_at_t: ACMGClass
    class_at_t_plus_n: ACMGClass
    #: Free-form key into :data:`EVIDENCE_TO_BLOCKING`; empty when unrecorded.
    evidence_type: str = ""
    resolved_on: date | None = None

    @property
    def was_resolved(self) -> bool:
        return (
            self.class_at_t is ACMGClass.UNCERTAIN
            and self.class_at_t_plus_n is not ACMGClass.UNCERTAIN
        )

    @property
    def direction(self) -> str | None:
        if self.class_at_t_plus_n in (ACMGClass.PATHOGENIC, ACMGClass.LIKELY_PATHOGENIC):
            return "pathogenic"
        if self.class_at_t_plus_n in (ACMGClass.BENIGN, ACMGClass.LIKELY_BENIGN):
            return "benign"
        return None


@dataclass(slots=True)
class ValidationResult:
    """The four metrics, plus the material to inspect the failures."""

    considered: int = 0
    resolved: int = 0
    resolvable_recall_numerator: int = 0
    cause_correct: int = 0
    cause_scored: int = 0
    direction_correct: int = 0
    direction_scored: int = 0
    #: Spearman rank correlation between gap size and days to resolution.
    temporal_correlation: float | None = None
    #: Cases worth a human look, never summarised away.
    misses: list[str] = field(default_factory=list)
    cause_confusion: dict[tuple[str, str], int] = field(default_factory=dict)

    @property
    def resolvability_recall(self) -> float | None:
        if self.resolved == 0:
            return None
        return self.resolvable_recall_numerator / self.resolved

    @property
    def cause_accuracy(self) -> float | None:
        if self.cause_scored == 0:
            return None
        return self.cause_correct / self.cause_scored

    @property
    def direction_accuracy(self) -> float | None:
        if self.direction_scored == 0:
            return None
        return self.direction_correct / self.direction_scored

    def summary(self) -> str:
        def pct(value: float | None) -> str:
            return "n/a" if value is None else f"{value:.1%}"

        return "\n".join(
            [
                f"variants considered      {self.considered}",
                f"resolved by T+n          {self.resolved}",
                f"1. resolvability recall  {pct(self.resolvability_recall)}",
                f"2. cause accuracy        {pct(self.cause_accuracy)} "
                f"(n={self.cause_scored})",
                f"3. direction accuracy    {pct(self.direction_accuracy)} "
                f"(n={self.direction_scored})",
                "4. temporal calibration  "
                + (
                    "n/a"
                    if self.temporal_correlation is None
                    else f"rho={self.temporal_correlation:+.3f} "
                    "(negative means a smaller gap resolved sooner)"
                ),
            ]
        )


def read_outcomes(path: str | Path) -> list[Outcome]:
    """Read a hindsight outcome table.

    Columns: ``variant_id``, ``class_at_t``, ``class_at_t_plus_n``, and
    optionally ``evidence_type`` and ``resolved_on`` (ISO date).
    """
    outcomes: list[Outcome] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            resolved_on = row.get("resolved_on") or ""
            outcomes.append(
                Outcome(
                    variant_id=row["variant_id"],
                    class_at_t=ACMGClass(row["class_at_t"]),
                    class_at_t_plus_n=ACMGClass(row["class_at_t_plus_n"]),
                    evidence_type=(row.get("evidence_type") or "").strip(),
                    resolved_on=date.fromisoformat(resolved_on) if resolved_on else None,
                )
            )
    return outcomes


def is_resolvable(row: GapMapRow) -> bool:
    """Whether the map said this variant could be resolved at all.

    A finite gap is not enough: the gap has to be closable by evidence somebody
    could actually obtain, so a variant whose only sufficient sets are
    intractable counts as *not* resolvable.
    """
    if row.gap_to_LP is None or row.gap_to_LB is None:
        return True  # already at or past a verdict in one direction
    return any(
        evidence_set.feasibility is not FeasibilityTag.INTRACTABLE
        for evidence_set in row.minimum_sufficient_sets
    )


def predicted_direction(row: GapMapRow) -> str | None:
    """Which way the map leaned: the target of its cheapest sufficient set."""
    if not row.minimum_sufficient_sets:
        return None
    cheapest = min(row.minimum_sufficient_sets, key=lambda s: s.acquisition_cost)
    return (
        "pathogenic" if cheapest.target is ACMGClass.LIKELY_PATHOGENIC else "benign"
    )


def _rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def spearman(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Rank correlation, implemented here to keep scipy out of the dependency set."""
    if len(xs) != len(ys) or len(xs) < 3:
        return None
    rx, ry = _rank(xs), _rank(ys)
    n = len(xs)
    mean_x, mean_y = sum(rx) / n, sum(ry) / n
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(rx, ry))
    denominator = (
        sum((a - mean_x) ** 2 for a in rx) * sum((b - mean_y) ** 2 for b in ry)
    ) ** 0.5
    if denominator == 0:
        return None
    return numerator / denominator


def validate(
    rows: Iterable[GapMapRow],
    outcomes: Iterable[Outcome],
    *,
    reference_date: date | None = None,
    evidence_mapping: dict[str, BlockingReason] | None = None,
) -> ValidationResult:
    """Run the section 10 protocol.

    ``rows`` must be the map computed with the evidence state of T, and
    ``outcomes`` what ClinVar recorded by T+n.
    """
    mapping = evidence_mapping or EVIDENCE_TO_BLOCKING
    by_id = {row.variant_id: row for row in rows}
    result = ValidationResult()

    gaps: list[float] = []
    delays: list[float] = []

    for outcome in outcomes:
        row = by_id.get(outcome.variant_id)
        if row is None:
            continue
        result.considered += 1
        if not outcome.was_resolved:
            continue
        result.resolved += 1

        if is_resolvable(row):
            result.resolvable_recall_numerator += 1
        else:
            result.misses.append(
                f"{outcome.variant_id}: resolved to "
                f"{outcome.class_at_t_plus_n.value} but the map offered no feasible "
                f"evidence set (blocking_reason={row.blocking_reason.value})"
            )

        if outcome.evidence_type:
            expected = mapping.get(outcome.evidence_type)
            if expected is not None:
                result.cause_scored += 1
                key = (row.blocking_reason.value, expected.value)
                result.cause_confusion[key] = result.cause_confusion.get(key, 0) + 1
                if row.blocking_reason is expected:
                    result.cause_correct += 1

        actual = outcome.direction
        predicted = predicted_direction(row)
        if actual and predicted:
            result.direction_scored += 1
            if actual == predicted:
                result.direction_correct += 1

        if reference_date and outcome.resolved_on:
            gap = row.gap_to_LP if actual == "pathogenic" else row.gap_to_LB
            if gap is not None:
                gaps.append(float(gap))
                delays.append(float((outcome.resolved_on - reference_date).days))

    result.temporal_correlation = spearman(gaps, delays)
    return result
