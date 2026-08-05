"""Layer 3: does the classification behave the way a point system must?

Failures here produce classifications that are plausible and wrong, so the
assertions are about structure -- monotonicity, non-overestimation, exclusivity
-- rather than about any particular verdict.
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from vus_foresight.acmg import (
    ACMGClass,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    PointSystem,
    SkipReason,
    Strength,
)
from vus_foresight.enumeration import enumerate_coding_snvs

CLASS_ORDER = {
    ACMGClass.BENIGN: 0,
    ACMGClass.LIKELY_BENIGN: 1,
    ACMGClass.UNCERTAIN: 2,
    ACMGClass.LIKELY_PATHOGENIC: 3,
    ACMGClass.PATHOGENIC: 4,
}


@settings(max_examples=300, deadline=None)
@given(st.integers(min_value=-40, max_value=40), st.integers(min_value=0, max_value=12))
def test_adding_pathogenic_points_never_moves_the_class_towards_benign(base, added):
    system = PointSystem()
    before = system.classify(base)
    after = system.classify(base + added)
    assert CLASS_ORDER[after] >= CLASS_ORDER[before]


@settings(max_examples=300, deadline=None)
@given(st.integers(min_value=-40, max_value=40), st.integers(min_value=0, max_value=12))
def test_adding_benign_points_never_moves_the_class_towards_pathogenic(base, added):
    system = PointSystem()
    assert CLASS_ORDER[system.classify(base - added)] <= CLASS_ORDER[system.classify(base)]


@settings(max_examples=300, deadline=None)
@given(st.integers(min_value=-40, max_value=40))
def test_gap_arithmetic_is_consistent_with_classification(points):
    system = PointSystem()
    gap_lp = system.gap_to_likely_pathogenic(points)
    if gap_lp is None:
        assert CLASS_ORDER[system.classify(points)] >= CLASS_ORDER[ACMGClass.LIKELY_PATHOGENIC]
    else:
        assert gap_lp > 0
        assert system.classify(points + gap_lp) in (
            ACMGClass.LIKELY_PATHOGENIC,
            ACMGClass.PATHOGENIC,
        )

    gap_lb = system.gap_to_likely_benign(points)
    if gap_lb is None:
        assert CLASS_ORDER[system.classify(points)] <= CLASS_ORDER[ACMGClass.LIKELY_BENIGN]
    else:
        assert gap_lb > 0
        assert system.classify(points - gap_lb) in (ACMGClass.LIKELY_BENIGN, ACMGClass.BENIGN)


def test_thresholds_reproduce_the_published_bands():
    system = PointSystem()
    assert system.classify(10) is ACMGClass.PATHOGENIC
    assert system.classify(9) is ACMGClass.LIKELY_PATHOGENIC
    assert system.classify(6) is ACMGClass.LIKELY_PATHOGENIC
    assert system.classify(5) is ACMGClass.UNCERTAIN
    assert system.classify(0) is ACMGClass.UNCERTAIN
    assert system.classify(-1) is ACMGClass.LIKELY_BENIGN
    assert system.classify(-6) is ACMGClass.LIKELY_BENIGN
    assert system.classify(-7) is ACMGClass.BENIGN


def test_point_values_match_the_tavtigian_table():
    system = PointSystem()
    assert system.points_for(Direction.PATHOGENIC, Strength.SUPPORTING) == 1
    assert system.points_for(Direction.PATHOGENIC, Strength.MODERATE) == 2
    assert system.points_for(Direction.PATHOGENIC, Strength.STRONG) == 4
    assert system.points_for(Direction.PATHOGENIC, Strength.VERY_STRONG) == 8
    assert system.points_for(Direction.BENIGN, Strength.SUPPORTING) == -1
    assert system.points_for(Direction.BENIGN, Strength.STRONG) == -4
    assert system.points_for(Direction.BENIGN, Strength.STAND_ALONE) == -8


# --------------------------------------------------------------------------
# Trace completeness -- principle 1: without the trace there is no map.
# --------------------------------------------------------------------------


def _run(gene, config, spec, limit=60):
    from datetime import datetime

    from vus_foresight.engine.pipeline import MapRunner, default_registry

    runner = MapRunner(
        transcript=gene.transcript,
        gene=config,
        spec=spec,
        adapters=default_registry(config),
        computed_at=datetime(1970, 1, 1),
    )
    variants = list(enumerate_coding_snvs(gene.transcript))[:limit]
    return list(runner.run(variants))


def test_every_criterion_appears_in_every_trace(minus_gene, minus_config, toy_spec):
    codes = {c.code for c in toy_spec.criteria}
    for result in _run(minus_gene, minus_config, toy_spec):
        seen = {c.code for c in result.row.criteria_applied} | {
            s.code for s in result.row.criteria_evaluated_not_applied
        }
        assert seen == codes, f"{result.row.hgvs_c} is missing {codes - seen}"


def test_extrinsic_criteria_are_never_applied(minus_gene, minus_config, toy_spec):
    """This system never sees a patient, so it can never satisfy extrinsic evidence."""
    for result in _run(minus_gene, minus_config, toy_spec):
        for applied in result.row.criteria_applied:
            assert applied.evidence_class is not EvidenceClass.EXTRINSIC
        extrinsic = [
            s
            for s in result.row.criteria_evaluated_not_applied
            if s.evidence_class is EvidenceClass.EXTRINSIC
        ]
        assert extrinsic
        assert all(s.reason is SkipReason.EXTRINSIC_EVIDENCE_REQUIRED for s in extrinsic)


def test_missing_data_is_never_reported_as_rule_not_met(minus_gene, minus_config, toy_spec):
    """Absent and false are different evidence states and must stay different."""
    for result in _run(minus_gene, minus_config, toy_spec):
        for skipped in result.row.criteria_evaluated_not_applied:
            if skipped.reason is SkipReason.DATA_MISSING:
                assert skipped.outcome is CriterionOutcome.NOT_EVALUABLE
                assert skipped.missing_fields, skipped.code
            if skipped.reason is SkipReason.RULE_NOT_MET:
                assert skipped.outcome is CriterionOutcome.NOT_MET
                assert not skipped.missing_fields


def test_points_equal_the_sum_of_applied_criteria(minus_gene, minus_config, toy_spec):
    for result in _run(minus_gene, minus_config, toy_spec):
        assert result.row.points_current == sum(c.points for c in result.row.criteria_applied)
        assert result.row.class_current is toy_spec.point_system.classify(result.row.points_current)


def test_mutually_exclusive_criteria_never_apply_together(minus_gene, minus_config, toy_spec):
    groups = {c.code: c.mutex_group for c in toy_spec.criteria if c.mutex_group}
    for result in _run(minus_gene, minus_config, toy_spec):
        applied_groups = [groups[c.code] for c in result.row.criteria_applied if c.code in groups]
        assert len(applied_groups) == len(set(applied_groups))


# --------------------------------------------------------------------------
# The gap itself.
# --------------------------------------------------------------------------


def test_ceiling_is_never_below_the_current_score(minus_gene, minus_config, toy_spec):
    for result in _run(minus_gene, minus_config, toy_spec):
        assert result.row.points_ceiling_intrinsic >= result.row.points_current


def test_minimum_sufficient_sets_actually_close_the_gap(minus_gene, minus_config, toy_spec):
    for result in _run(minus_gene, minus_config, toy_spec):
        row = result.row
        for evidence_set in row.minimum_sufficient_sets:
            total = sum(abs(r.points) for r in evidence_set.requirements)
            gap = (
                row.gap_to_LP
                if evidence_set.target is ACMGClass.LIKELY_PATHOGENIC
                else row.gap_to_LB
            )
            assert gap is not None
            assert total >= gap, f"{evidence_set.codes} sums to {total}, gap is {gap}"


def test_minimum_sufficient_sets_are_minimal_by_inclusion(minus_gene, minus_config, toy_spec):
    for result in _run(minus_gene, minus_config, toy_spec):
        for target in (ACMGClass.LIKELY_PATHOGENIC, ACMGClass.LIKELY_BENIGN):
            sets = [s for s in result.row.minimum_sufficient_sets if s.target is target]
            for i, first in enumerate(sets):
                for j, second in enumerate(sets):
                    if i == j:
                        continue
                    assert not set(first.codes) < set(second.codes), (
                        f"{second.codes} contains the sufficient subset {first.codes}"
                    )


def test_intractable_evidence_is_excluded_from_proposed_sets(minus_gene, minus_config, toy_spec):
    """A set requiring a de novo event in adult-onset cancer is not a research plan."""
    from vus_foresight.acmg import FeasibilityTag

    assert not toy_spec.gap.include_intractable
    for result in _run(minus_gene, minus_config, toy_spec):
        for evidence_set in result.row.minimum_sufficient_sets:
            assert all(
                r.feasibility is not FeasibilityTag.INTRACTABLE for r in evidence_set.requirements
            )


def test_sets_are_ordered_cheapest_first(minus_gene, minus_config, toy_spec):
    for result in _run(minus_gene, minus_config, toy_spec):
        costs = [s.acquisition_cost for s in result.row.minimum_sufficient_sets]
        assert costs == sorted(costs)


def test_resolved_variants_are_not_marked_blocked(minus_gene, minus_config, toy_spec):
    from vus_foresight.acmg import BlockingReason

    for result in _run(minus_gene, minus_config, toy_spec):
        row = result.row
        if row.class_current is not ACMGClass.UNCERTAIN:
            assert row.blocking_reason is BlockingReason.RESOLVED_NOT_BLOCKED
        else:
            assert row.blocking_reason is not BlockingReason.RESOLVED_NOT_BLOCKED
