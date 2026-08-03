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
from typing import TYPE_CHECKING, Iterable, Sequence

from .acmg import ACMGClass, BlockingReason, EvidenceClass, FeasibilityTag
from .gapmap import GapMapRow

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .adapters.clinvar import ClinVarSnapshot
    from .engine.spec import VCEPSpec
    from .engine.timeline import TimelineDiff

__all__ = [
    "Outcome",
    "ValidationResult",
    "read_outcomes",
    "outcomes_from_snapshots",
    "observed_causes_from_diff",
    "is_resolvable",
    "predicted_direction",
    "spearman",
    "validate",
    "run_time_series_study",
    "TimeSeriesStudy",
    "CLINVAR_TO_ACMG",
    "NEIGHBOUR_CLASSIFICATION",
]

#: ClinVar's normalised vocabulary mapped onto the five-tier terminology.
#:
#: ``conflicting`` becomes VUS rather than being dropped. A variant with
#: conflicting submissions has not been resolved -- treating it as resolved
#: would inflate every metric in this module.
CLINVAR_TO_ACMG: dict[str, ACMGClass] = {
    "pathogenic": ACMGClass.PATHOGENIC,
    "likely_pathogenic": ACMGClass.LIKELY_PATHOGENIC,
    "uncertain": ACMGClass.UNCERTAIN,
    "conflicting": ACMGClass.UNCERTAIN,
    "likely_benign": ACMGClass.LIKELY_BENIGN,
    "benign": ACMGClass.BENIGN,
}

#: Observed cause for a resolution driven by a neighbour being classified.
#:
#: Deliberately *not* a member of :class:`BlockingReason`. That enum is an output
#: column and its vocabulary is fixed by the schema; widening it would change
#: what every existing ``GROUP BY`` means. This label lives only in the study,
#: where the observed vocabulary can be one item wider than the predicted one.
#:
#: Its existence is itself a finding. No blocking reason in the schema says "no
#: neighbour has been classified yet", so a variant blocked on exactly that --
#: the temporally dynamic case this whole project is built around -- currently
#: has to be filed under something else. Metric 2 will show that as a systematic
#: miss, which is the right way for the gap to surface.
NEIGHBOUR_CLASSIFICATION = "NEIGHBOUR_CLASSIFICATION"

#: Which blocking reason each kind of newly-submitted evidence corresponds to.
#: A study can override it; the default reflects how ClinVar submissions
#: describe their evidence.
EVIDENCE_TO_BLOCKING: dict[str, str] = {
    "functional": BlockingReason.MISSING_FUNCTIONAL.value,
    "segregation": BlockingReason.MISSING_SEGREGATION.value,
    "case_control": BlockingReason.MISSING_CASE_CONTROL.value,
    "population_frequency": BlockingReason.FREQUENCY_UNINFORMATIVE.value,
    "computational": BlockingReason.PREDICTOR_INDETERMINATE.value,
    "neighbour_classification": NEIGHBOUR_CLASSIFICATION,
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
    #: Resolutions for which no evidence reached this pipeline's sources at all.
    #: Not a failure of the map -- a measurement of how much of the field's
    #: resolution happens on evidence that never becomes public data.
    unobserved_cause: int = 0
    direction_correct: int = 0
    direction_scored: int = 0
    #: Resolutions on which the map had accumulated no evidence either way, so
    #: it made no directional claim to score.
    direction_unpredicted: int = 0
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
                f"(n={self.cause_scored}"
                + (
                    f", {self.unobserved_cause} resolved on evidence this "
                    "pipeline never sees)"
                    if self.unobserved_cause
                    else ")"
                ),
                f"3. direction accuracy    {pct(self.direction_accuracy)} "
                f"(n={self.direction_scored}"
                + (
                    f", {self.direction_unpredicted} with no evidence either way)"
                    if self.direction_unpredicted
                    else ")"
                ),
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


def outcomes_from_snapshots(
    before: "ClinVarSnapshot",
    after: "ClinVarSnapshot",
    *,
    transcript_id: str,
    min_stars: int = 1,
) -> list[Outcome]:
    """Derive the ground truth from the two dated ClinVar snapshots themselves.

    This works because of a property of the design rather than a convenience:
    ``clinvar.self.classification`` -- ClinVar's own call on the variant -- is
    published by the adapter and read by **no criterion**, enforced by a test.
    The oracle for the study is therefore sitting in the same snapshots the
    engine consumes, and is firewalled from the thing being measured.

    Only variants *present and uncertain* at T are eligible. A variant absent
    from the earlier snapshot was not a VUS then, it was unexamined, and
    counting it would answer a different question from the one section 10 asks.
    """
    outcomes: list[Outcome] = []
    for hgvs_c, old in sorted(before.by_nucleotide.items()):
        if old.stars < min_stars:
            continue
        class_before = CLINVAR_TO_ACMG.get(old.classification)
        if class_before is not ACMGClass.UNCERTAIN:
            continue
        new = after.by_nucleotide.get(hgvs_c)
        if new is None:
            # Withdrawn or reclassified out of the archive. Not a resolution.
            continue
        class_after = CLINVAR_TO_ACMG.get(new.classification)
        if class_after is None:
            continue
        outcomes.append(
            Outcome(
                variant_id=f"{transcript_id}:{hgvs_c}",
                class_at_t=class_before,
                class_at_t_plus_n=class_after,
                resolved_on=new.last_evaluated,
            )
        )
    return outcomes


def observed_causes_from_diff(
    diff: "TimelineDiff", spec: "VCEPSpec"
) -> dict[str, str]:
    """What kind of evidence actually turned up, per variant.

    Metric 2 asks whether the predicted ``blocking_reason`` matches the evidence
    that actually appeared. The obvious source for that is the free text of the
    submission record, which would need parsing prose -- and this project uses
    no LLM at any step, so that route means a hand-written keyword classifier
    nobody can calibrate.

    There is an exact source instead. Every criterion already declares
    ``blocks_as``: the blocking reason its absence causes. So a criterion that
    *newly applied* between the two snapshots is, by the specification's own
    mapping, the evidence that relieved that block. No new table, no prose.

    The limitation is worth stating: this measures the evidence that reached
    *this pipeline's sources*, not the evidence a submitter cited. The two
    coincide for functional, frequency and computational data, which arrive here
    from the same public datasets. They can diverge for a submitter's unpublished
    segregation or case-control data, which never reaches an intrinsic criterion
    at all -- so a variant the field resolved on family data shows up here as
    having no observed cause, which is the honest answer rather than a guess.
    """
    causes: dict[str, str] = {}
    for transition in diff.transitions:
        labels = [
            label
            for code, _strength in transition.criteria_gained
            if (label := _observed_label(spec, code)) is not None
        ]
        if not labels:
            continue
        # Rank by the spec's own priority so a variant relieved on two fronts is
        # attributed the same way the blocking reason itself would have been.
        ranking = [reason.value for reason in spec.blocking.priority]
        ranking.append(NEIGHBOUR_CLASSIFICATION)
        for candidate in ranking:
            if candidate in labels:
                causes[transition.variant_id] = candidate
                break
        else:
            causes[transition.variant_id] = labels[0]
    return causes


def _observed_label(spec: "VCEPSpec", code: str) -> str | None:
    """What kind of evidence a newly-applied criterion represents."""
    try:
        criterion = spec.by_code(code)
    except KeyError:
        return None
    if criterion.blocks_as is not None:
        return criterion.blocks_as.value
    if criterion.evidence_class is EvidenceClass.SEMI_INTRINSIC:
        # PS1 and PM5 have no blocking reason to declare, because the schema has
        # no member for "no neighbour classified yet". They are still the most
        # informative thing that can move a variant, so they are labelled rather
        # than dropped.
        return NEIGHBOUR_CLASSIFICATION
    return None


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
    """Which way the map leaned, or ``None`` when it did not lean.

    The sign of the accumulated points, and nothing else.

    Two tempting definitions are wrong and were tried first. The target of the
    *cheapest* sufficient set measures acquisition cost, not evidence: with no
    data loaded, the cheapest route is whichever criterion is easiest to obtain,
    which says nothing about where the variant is heading. And comparing
    ``gap_to_LP`` against ``gap_to_LB`` inherits the asymmetry of the
    thresholds -- at zero points the gaps are 6 and 1, so that rule would call
    every unevidenced variant benign.

    A variant with no applied criteria genuinely points nowhere, so it is not
    scored. That shrinks metric 3's denominator to the variants on which the map
    actually made a directional claim, which is the only population the metric
    means anything for.
    """
    if row.points_current > 0:
        return "pathogenic"
    if row.points_current < 0:
        return "benign"
    return None


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
    evidence_mapping: dict[str, str] | None = None,
    observed_causes: dict[str, str] | None = None,
) -> ValidationResult:
    """Run the section 10 protocol.

    ``rows`` must be the map computed with the evidence state of T, and
    ``outcomes`` what ClinVar recorded by T+n.

    ``observed_causes`` supplies metric 2 from the map's own time series -- see
    :func:`observed_causes_from_diff`. When both it and an outcome's
    ``evidence_type`` are present, the curated ``evidence_type`` wins: a human
    who read the submission record knows something the pipeline's sources do not.
    """
    mapping = evidence_mapping or EVIDENCE_TO_BLOCKING
    # ``None`` and ``{}`` mean different things and must stay apart: no
    # time-series mode at all, versus time-series mode in which nothing moved.
    # Collapsing them would report zero unobserved causes precisely when every
    # cause is unobserved.
    time_series_mode = observed_causes is not None
    observed = observed_causes or {}
    outcomes = list(outcomes)
    # Only the variants ClinVar has an opinion about are ever looked up, and
    # that is a few thousand out of a whole gene's enumeration. Indexing just
    # those lets ``rows`` be a stream rather than a list -- the difference
    # between 9 GB and a few megabytes on real BRCA1.
    wanted = {outcome.variant_id for outcome in outcomes}
    by_id = {row.variant_id: row for row in rows if row.variant_id in wanted}
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

        expected = (
            mapping.get(outcome.evidence_type)
            if outcome.evidence_type
            else observed.get(outcome.variant_id)
        )
        if expected is not None:
            result.cause_scored += 1
            key = (row.blocking_reason.value, expected)
            result.cause_confusion[key] = result.cause_confusion.get(key, 0) + 1
            if row.blocking_reason.value == expected:
                result.cause_correct += 1
        elif time_series_mode:
            # The field resolved it on evidence this pipeline never sees --
            # a submitter's segregation or case-control data. Counted, because
            # a systematically large share here is itself a finding about how
            # much of the resolution happens outside public datasets.
            result.unobserved_cause += 1

        actual = outcome.direction
        predicted = predicted_direction(row)
        if actual and predicted:
            result.direction_scored += 1
            if actual == predicted:
                result.direction_correct += 1
        elif actual:
            result.direction_unpredicted += 1

        if reference_date and outcome.resolved_on:
            gap = row.gap_to_LP if actual == "pathogenic" else row.gap_to_LB
            if gap is not None:
                gaps.append(float(gap))
                delays.append(float((outcome.resolved_on - reference_date).days))

    result.temporal_correlation = spearman(gaps, delays)
    return result


@dataclass(slots=True)
class TimeSeriesStudy:
    """The section 10 protocol run over the map's own recomputation.

    Three things are compared, and keeping them distinct is the whole point:

    * what the map **predicted** at T -- the gap, the blocking reason, the
      direction of the cheapest sufficient set;
    * what the map itself **did** between T and T+n, once new snapshots landed;
    * what ClinVar **concluded** by T+n, which no criterion is allowed to read.
    """

    result: ValidationResult
    #: Resolved in ClinVar, and the map reached a verdict too. The map and the
    #: field agreed, on evidence both could see.
    anticipated: list[str] = field(default_factory=list)
    #: Resolved in ClinVar; the map registered the same new evidence but it was
    #: not enough to cross a threshold. The most diagnostic bucket of the three:
    #: either this specification is stricter than the submitter was, or the
    #: submitter had evidence the public datasets do not carry.
    saw_evidence_only: list[str] = field(default_factory=list)
    #: Resolved in ClinVar with nothing moving here at all. The evidence never
    #: became public data, or a source is missing from the pipeline.
    unmoved: list[str] = field(default_factory=list)
    #: The map moved and ClinVar has not. Not errors -- these are the map's
    #: live, still-unfalsified predictions, and the rows worth reading first.
    ahead_of_clinvar: list[str] = field(default_factory=list)

    @property
    def resolved_total(self) -> int:
        return len(self.anticipated) + len(self.saw_evidence_only) + len(self.unmoved)

    @property
    def anticipation_rate(self) -> float | None:
        """Fraction of the field's resolutions the map also brought to a verdict."""
        return len(self.anticipated) / self.resolved_total if self.resolved_total else None

    @property
    def evidence_visibility_rate(self) -> float | None:
        """Fraction where the map at least saw the same evidence arrive.

        Separating this from the verdict rate matters: a low verdict rate with a
        high visibility rate means the thresholds are the disagreement, while a
        low visibility rate means the data never reached the pipeline. Those
        call for completely different work.
        """
        if not self.resolved_total:
            return None
        return (
            len(self.anticipated) + len(self.saw_evidence_only)
        ) / self.resolved_total

    def summary(self) -> str:
        def pct(value: float | None) -> str:
            return "n/a" if value is None else f"{value:.1%}"

        return "\n".join(
            [
                self.result.summary(),
                f"5. map reached a verdict {pct(self.anticipation_rate)} "
                f"({len(self.anticipated)} of {self.resolved_total})",
                f"   map saw the evidence  {pct(self.evidence_visibility_rate)} "
                f"({len(self.saw_evidence_only)} saw it without crossing a threshold, "
                f"{len(self.unmoved)} never saw it)",
                f"   map ahead of ClinVar  {len(self.ahead_of_clinvar)} "
                "(moved here, not yet in the archive)",
            ]
        )


def run_time_series_study(
    rows_at_t: Iterable[GapMapRow],
    rows_at_t_plus_n: Iterable[GapMapRow],
    clinvar_at_t: "ClinVarSnapshot",
    clinvar_at_t_plus_n: "ClinVarSnapshot",
    spec: "VCEPSpec",
    *,
    transcript_id: str,
    reference_date: date | None = None,
    min_stars: int = 1,
    curated_outcomes: Iterable[Outcome] | None = None,
) -> TimeSeriesStudy:
    """Run section 10 with no external outcome table.

    Everything the protocol needs is already on disk: two maps, two ClinVar
    snapshots, and the specification. ``curated_outcomes`` still takes
    precedence when supplied, because a human who read the submission records
    knows what the pipeline's sources cannot.

    The map at T is walked twice -- once to diff it, once to look up the
    variants ClinVar resolved -- so it must be re-readable. Pass a list, or a
    :class:`~vus_foresight.output.GapMapSource` for a whole gene.
    """
    from .engine.timeline import compare_maps

    for name, rows in (("rows_at_t", rows_at_t), ("rows_at_t_plus_n", rows_at_t_plus_n)):
        if iter(rows) is rows:
            raise TypeError(
                f"{name} is a one-shot iterator; pass a list or a GapMapSource. "
                "Materialising it here would defeat the point of streaming."
            )

    diff = compare_maps(rows_at_t, rows_at_t_plus_n)
    observed = observed_causes_from_diff(diff, spec)

    outcomes = list(curated_outcomes) if curated_outcomes is not None else []
    if not outcomes:
        outcomes = outcomes_from_snapshots(
            clinvar_at_t,
            clinvar_at_t_plus_n,
            transcript_id=transcript_id,
            min_stars=min_stars,
        )

    result = validate(
        rows_at_t,
        outcomes,
        reference_date=reference_date,
        observed_causes=observed,
    )

    reached_verdict = {t.variant_id for t in diff.class_changes}
    saw_evidence = {t.variant_id for t in diff.transitions} - reached_verdict
    resolved_in_archive = {o.variant_id for o in outcomes if o.was_resolved}

    study = TimeSeriesStudy(result=result)
    study.anticipated = sorted(resolved_in_archive & reached_verdict)
    study.saw_evidence_only = sorted(resolved_in_archive & saw_evidence)
    study.unmoved = sorted(resolved_in_archive - reached_verdict - saw_evidence)
    study.ahead_of_clinvar = sorted(reached_verdict - resolved_in_archive)
    return study
