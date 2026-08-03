"""Recomputing the map over time, and attributing what moved (spec section 4).

The semi-intrinsic criteria make the map temporally dynamic: PS1 and PM5 depend
on the state of a public database on a date, so a variant can leave VUS
**without any new evidence about itself appearing** -- purely because a
different variant in the same codon got classified. That observation is the
reason the intrinsic/semi-intrinsic distinction exists, and it is only visible
if the map is recomputed against dated snapshots and the results are diffed.

So the diff does not just report that a variant moved. It reports *what kind of
evidence moved it*, which is the difference between a changelog and a finding.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Sequence

from ..acmg import ACMGClass, BlockingReason, EvidenceClass
from ..gapmap import GapMapRow

__all__ = [
    "TransitionCause",
    "Transition",
    "TimelineDiff",
    "compare_maps",
    "compare_series",
]


class TransitionCause(str, Enum):
    """Why a variant's classification moved between two snapshots."""

    #: Only semi-intrinsic criteria changed. Nothing new was learned about this
    #: variant -- a *neighbour* was classified. The section 4 headline.
    NEIGHBOUR_EVIDENCE = "neighbour_evidence"
    #: Intrinsic criteria changed: new frequency, predictor or assay data.
    INTRINSIC_EVIDENCE = "intrinsic_evidence"
    #: Both moved in the same interval; the two cannot be separated.
    MIXED_EVIDENCE = "mixed_evidence"
    #: The same criteria applied, at different strengths.
    STRENGTH_ONLY = "strength_only"
    #: The class moved with no criterion delta at all, which means the
    #: specification or its thresholds changed rather than the evidence.
    SPEC_CHANGE = "spec_change"


@dataclass(frozen=True, slots=True)
class Transition:
    """One variant, between two snapshots."""

    variant_id: str
    gene: str
    hgvs_c: str
    hgvs_p: str | None
    consequence: str
    class_before: ACMGClass
    class_after: ACMGClass
    points_before: int
    points_after: int
    blocking_before: BlockingReason
    blocking_after: BlockingReason
    gap_to_lp_before: int | None
    gap_to_lp_after: int | None
    criteria_gained: tuple[tuple[str, str], ...]
    criteria_lost: tuple[tuple[str, str], ...]
    cause: TransitionCause
    spec_before: str
    spec_after: str

    @property
    def class_changed(self) -> bool:
        return self.class_before is not self.class_after

    @property
    def resolved(self) -> bool:
        return (
            self.class_before is ACMGClass.UNCERTAIN
            and self.class_after is not ACMGClass.UNCERTAIN
        )

    @property
    def resolved_without_self_evidence(self) -> bool:
        """Left VUS on the strength of a neighbour's classification alone.

        The falsifiable, publishable observation: nobody generated any evidence
        about this variant, and its classification changed anyway.
        """
        return self.resolved and self.cause is TransitionCause.NEIGHBOUR_EVIDENCE


def _applied_index(row: GapMapRow) -> dict[str, tuple[str, EvidenceClass]]:
    return {
        criterion.code: (criterion.strength.value, criterion.evidence_class)
        for criterion in row.criteria_applied
    }


def _attribute(
    gained: Sequence[tuple[str, str, EvidenceClass]],
    lost: Sequence[tuple[str, str, EvidenceClass]],
    strength_changed: bool,
    spec_changed: bool,
) -> TransitionCause:
    classes = {entry[2] for entry in list(gained) + list(lost)}
    if not classes:
        if strength_changed:
            return TransitionCause.STRENGTH_ONLY
        # No criterion moved and no strength moved, yet something did. The only
        # remaining explanation is that the rules themselves changed; the
        # transition records both spec versions so a reviewer can confirm it.
        del spec_changed
        return TransitionCause.SPEC_CHANGE
    semi = EvidenceClass.SEMI_INTRINSIC in classes
    other = bool(classes - {EvidenceClass.SEMI_INTRINSIC})
    if semi and other:
        return TransitionCause.MIXED_EVIDENCE
    if semi:
        return TransitionCause.NEIGHBOUR_EVIDENCE
    return TransitionCause.INTRINSIC_EVIDENCE


@dataclass(slots=True)
class TimelineDiff:
    """Everything that changed between two dated maps."""

    label_before: str
    label_after: str
    transitions: list[Transition] = field(default_factory=list)
    appeared: list[str] = field(default_factory=list)
    disappeared: list[str] = field(default_factory=list)
    compared: int = 0

    @property
    def class_changes(self) -> list[Transition]:
        return [t for t in self.transitions if t.class_changed]

    @property
    def resolutions(self) -> list[Transition]:
        return [t for t in self.transitions if t.resolved]

    @property
    def neighbour_driven_resolutions(self) -> list[Transition]:
        return [t for t in self.transitions if t.resolved_without_self_evidence]

    def counts_by_cause(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for transition in self.class_changes:
            counts[transition.cause.value] = counts.get(transition.cause.value, 0) + 1
        return dict(sorted(counts.items()))

    def summary(self) -> str:
        lines = [
            f"{self.label_before} -> {self.label_after}",
            f"  compared              {self.compared:,}",
            f"  criteria moved        {len(self.transitions):,}",
            f"  class changed         {len(self.class_changes):,}",
            f"  left VUS              {len(self.resolutions):,}",
            f"  ... on a neighbour's  {len(self.neighbour_driven_resolutions):,} "
            "(no new evidence about the variant itself)",
        ]
        if self.appeared or self.disappeared:
            lines.append(
                f"  rows appeared {len(self.appeared):,}, "
                f"disappeared {len(self.disappeared):,}"
            )
        for cause, count in self.counts_by_cause().items():
            lines.append(f"    {cause:22s} {count:>8,}")
        return "\n".join(lines)


def compare_maps(
    before: Iterable[GapMapRow],
    after: Iterable[GapMapRow],
    *,
    label_before: str = "before",
    label_after: str = "after",
    include_unchanged: bool = False,
) -> TimelineDiff:
    """Diff two maps of the same gene computed at different dates.

    Rows are matched on ``variant_id``. A variant present in only one snapshot
    is reported separately rather than treated as a change: it means the
    enumeration or the transcript moved, which is a different event from
    evidence moving.
    """
    before_rows = {row.variant_id: row for row in before}
    after_rows = {row.variant_id: row for row in after}
    diff = TimelineDiff(label_before=label_before, label_after=label_after)

    diff.appeared = sorted(set(after_rows) - set(before_rows))
    diff.disappeared = sorted(set(before_rows) - set(after_rows))

    for variant_id in sorted(set(before_rows) & set(after_rows)):
        old, new = before_rows[variant_id], after_rows[variant_id]
        diff.compared += 1

        old_applied, new_applied = _applied_index(old), _applied_index(new)
        gained = [
            (code, strength, evidence_class)
            for code, (strength, evidence_class) in sorted(new_applied.items())
            if code not in old_applied
        ]
        lost = [
            (code, strength, evidence_class)
            for code, (strength, evidence_class) in sorted(old_applied.items())
            if code not in new_applied
        ]
        strength_changed = any(
            old_applied[code][0] != new_applied[code][0]
            for code in set(old_applied) & set(new_applied)
        )
        moved = bool(gained or lost) or strength_changed
        if not moved and old.class_current is new.class_current and not include_unchanged:
            continue

        diff.transitions.append(
            Transition(
                variant_id=variant_id,
                gene=new.gene,
                hgvs_c=new.hgvs_c,
                hgvs_p=new.hgvs_p,
                consequence=new.consequence.value,
                class_before=old.class_current,
                class_after=new.class_current,
                points_before=old.points_current,
                points_after=new.points_current,
                blocking_before=old.blocking_reason,
                blocking_after=new.blocking_reason,
                gap_to_lp_before=old.gap_to_LP,
                gap_to_lp_after=new.gap_to_LP,
                criteria_gained=tuple((code, strength) for code, strength, _ in gained),
                criteria_lost=tuple((code, strength) for code, strength, _ in lost),
                cause=_attribute(
                    gained,
                    lost,
                    strength_changed,
                    old.spec_version != new.spec_version,
                ),
                spec_before=old.spec_version,
                spec_after=new.spec_version,
            )
        )
    return diff


def compare_series(
    snapshots: Sequence[tuple[str, Iterable[GapMapRow]]],
) -> list[TimelineDiff]:
    """Diff consecutive snapshots of a series, oldest first.

    Consecutive rather than all-against-the-first: the question section 4 asks
    is when a variant moved and on what, and collapsing five years into a single
    before/after loses exactly that.
    """
    materialised = [(label, list(rows)) for label, rows in snapshots]
    if len(materialised) < 2:
        raise ValueError("a series needs at least two snapshots")
    return [
        compare_maps(
            older,
            newer,
            label_before=label_older,
            label_after=label_newer,
        )
        for (label_older, older), (label_newer, newer) in zip(
            materialised, materialised[1:]
        )
    ]
