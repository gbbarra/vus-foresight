"""The edges of the gap engine: the ceiling, the sets, and the blocking label.

Oracle for this file: the module's own stated contract, exercised against
specifications built inside the test rather than against ``config/specs``. Each
assertion names a concrete list of codes, a concrete feasibility tag or a
concrete :class:`BlockingReason`, so a shipped threshold moving cannot make one
of these tests quietly agree with a wrong answer.

How these tests could nevertheless be wrong: by describing a configuration no
YAML file could produce, or by asserting on a hand-built ``Evaluation`` that the
evaluator would never emit. Both are guarded the same way -- every specification
here is built through the very Pydantic models ``load_spec`` validates into,
with their validators running, and every evaluation comes out of
``evaluate_variant`` rather than being assembled by hand. What is deliberately
*not* realistic is the shipped configuration: ``config/specs/toy`` gives every
criterion an explicit ``feasibility`` and a mutex layout in which no group is
ever occupied while a sibling is still missing, so several branches below are
unreachable from the shipped files. That is exactly why they are here: the first
specification that omits a feasibility tag would otherwise discover an untested
default in the one column the whole project exists to report.
"""

from __future__ import annotations

import pytest

from vus_foresight.acmg import (
    ACMGClass,
    BlockingReason,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    PointSystem,
    SkipReason,
    Strength,
)
from vus_foresight.engine.context import EvidenceContext
from vus_foresight.engine.evaluator import evaluate_variant
from vus_foresight.engine.gap import (
    Candidate,
    analyse_gap,
    candidate_requirements,
    derive_blocking_reason,
    minimum_sufficient_sets,
)
from vus_foresight.engine.predicates import Leaf, RuleNode
from vus_foresight.engine.spec import (
    BlockingConfig,
    CriterionSpec,
    FeasibilityRung,
    GapConfig,
    VCEPSpec,
)
from vus_foresight.variant import Consequence, Variant, VariantKind

# --------------------------------------------------------------------------
# Builders. Deliberately thin: they set defaults, never behaviour.
# --------------------------------------------------------------------------


def _leaf(field: str, op: str, value=None) -> RuleNode:
    return RuleNode(leaf=Leaf(field=field, op=op, value=value))


def _criterion(code: str, **overrides) -> CriterionSpec:
    """An intrinsic criterion whose data is absent unless a test supplies it.

    The default rule reads a path no builder here ever populates, so a criterion
    is a *candidate* -- blocked by missing data -- until a test says otherwise.
    """
    fields = {
        "code": code,
        "direction": Direction.PATHOGENIC,
        "evidence_class": EvidenceClass.INTRINSIC,
        "requires": (f"evidence.{code.lower()}",),
        "rule": _leaf(f"evidence.{code.lower()}", "is_true"),
        "feasibility": FeasibilityTag.AVAILABLE_UNINGESTED,
    }
    fields.update(overrides)
    return CriterionSpec(**fields)


def _extrinsic(code: str, **overrides) -> CriterionSpec:
    fields = {
        "code": code,
        "direction": Direction.PATHOGENIC,
        "evidence_class": EvidenceClass.EXTRINSIC,
        "strength": Strength.STRONG,
        "feasibility": FeasibilityTag.NEEDS_CASES,
        "rule": None,
    }
    fields.update(overrides)
    return CriterionSpec(**fields)


def _spec(*criteria: CriterionSpec, **overrides) -> VCEPSpec:
    fields = {
        "spec_id": "gapedges",
        "name": "Specification built for the gap edge tests",
        "version": "0.0.0",
        "verified": True,
        "criteria": criteria,
    }
    fields.update(overrides)
    return VCEPSpec(**fields)


def _variant(consequence: Consequence = Consequence.MISSENSE) -> Variant:
    return Variant(
        gene="TOY1",
        transcript="NM_999003.1",
        kind=VariantKind.SNV,
        hgvs_c="c.10A>G",
        consequence=consequence,
    )


def _context(data: dict | None = None) -> EvidenceContext:
    return EvidenceContext(data if data is not None else {}, {})


def _candidates(spec: VCEPSpec, data: dict | None = None) -> list[Candidate]:
    """Candidates as the engine derives them, from a real evaluation."""
    context = _context(data)
    return candidate_requirements(evaluate_variant(_variant(), context, spec), spec, context)


#: ``count`` interchangeable pathogenic criteria worth one point each, all
#: blocked by missing data. Used where the question is about the *shape* of the
#: search rather than about any one criterion.
def _supporting_spec(count: int, **overrides) -> VCEPSpec:
    return _spec(*(_criterion(f"PP{i}") for i in range(1, count + 1)), **overrides)


# --------------------------------------------------------------------------
# Feasibility resolution: which column of the report a missing criterion
# lands in.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("declared", "with_ladder", "region_assayed", "expected"),
    [
        # The ladder is the reason the mechanism exists: PS3 inside a region an
        # assay already covers is somebody's supplementary table, not bench work.
        (FeasibilityTag.ASSAY_FEASIBLE, True, True, FeasibilityTag.AVAILABLE_UNINGESTED),
        (FeasibilityTag.ASSAY_FEASIBLE, True, False, FeasibilityTag.ASSAY_FEASIBLE),
        # An absent flag is not a false flag, and must not promote the criterion
        # into the "already published" bucket.
        (FeasibilityTag.ASSAY_FEASIBLE, True, None, FeasibilityTag.ASSAY_FEASIBLE),
        # No rung matched and the specification declared nothing: intrinsic
        # evidence needs no patient by definition, so the honest default is
        # "somebody already has this, we have not loaded it".
        (None, True, False, FeasibilityTag.AVAILABLE_UNINGESTED),
        (None, False, None, FeasibilityTag.AVAILABLE_UNINGESTED),
    ],
)
def test_feasibility_falls_from_the_ladder_to_the_declared_tag_to_the_intrinsic_default(
    declared, with_ladder, region_assayed, expected
):
    """This tag decides whether the report asks for engineering or for a bench.

    ``available_uningested`` is the single most actionable thing the system
    emits. Resolving it from a stale default instead of from the ladder would
    either invent published data that does not exist, or bury the finding the
    project exists to surface.
    """
    ladder = (
        (
            FeasibilityRung(
                feasibility=FeasibilityTag.AVAILABLE_UNINGESTED,
                when=_leaf("functional.region_assayed", "is_true"),
            ),
        )
        if with_ladder
        else ()
    )
    spec = _spec(_criterion("PS3", feasibility=declared, feasibility_ladder=ladder))
    data = {} if region_assayed is None else {"functional": {"region_assayed": region_assayed}}

    candidates = _candidates(spec, data)

    assert [c.code for c in candidates] == ["PS3"]
    assert candidates[0].requirement.feasibility is expected


# --------------------------------------------------------------------------
# Which criteria become candidates at all.
# --------------------------------------------------------------------------


def test_a_missing_criterion_is_not_offered_when_its_group_already_scored():
    """Proposing it would promise points the mutex resolution would then discard.

    BA1, BS1 and PM2 share one group: at most one of them can ever contribute.
    Once BA1 has fired, acquiring the data PM2 is waiting on cannot move the
    total by a single point, so listing it as a gap would advertise work with a
    guaranteed return of zero.
    """
    spec = _spec(
        _criterion(
            "BA1",
            direction=Direction.BENIGN,
            strength=Strength.STAND_ALONE,
            mutex_group="frequency",
            requires=("frequency.gnomad.faf95_popmax",),
            rule=_leaf("frequency.gnomad.faf95_popmax", "ge", 0.001),
        ),
        _criterion("PM2_Supporting", mutex_group="frequency"),
        _criterion("PP3", mutex_group="computational"),
    )
    data = {"frequency": {"gnomad": {"faf95_popmax": 0.005}}}

    candidates = _candidates(spec, data)

    assert evaluate_variant(_variant(), _context(data), spec).applied_codes == ("BA1",)
    assert [c.code for c in candidates] == ["PP3"]


def test_a_criterion_answered_no_is_not_a_gap_but_one_never_asked_is():
    """Absent data and a negative answer are different states of the world.

    Only the first can be closed by acquiring something. Conflating them is how
    a gap map starts recommending experiments whose result is already known.
    """
    spec = _spec(
        _criterion("PP3"),
        _criterion("PM1"),
        _extrinsic("PS4"),
    )
    data = {"evidence": {"pm1": False}}

    candidates = _candidates(spec, data)

    evaluation = evaluate_variant(_variant(), _context(data), spec)
    assert evaluation.outcome_of("PM1") is CriterionOutcome.NOT_MET
    assert [c.code for c in candidates] == ["PP3", "PS4"]


# --------------------------------------------------------------------------
# The search for a sufficient set.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("gap", [0, -1, -6])
def test_no_evidence_is_proposed_once_the_target_is_already_reached(gap):
    """A non-positive gap is a variant that already qualifies.

    The candidates below would close a gap of four, so an empty answer here is
    the guard doing its job rather than an empty candidate list. Proposing
    evidence for a class the variant already holds would put work into the
    report that nobody needs to do.
    """
    spec = _spec(_criterion("PS3", strength=Strength.STRONG))
    candidates = _candidates(spec)

    sets = minimum_sufficient_sets(
        candidates, gap, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )

    assert sets == ()
    assert minimum_sufficient_sets(
        candidates, 4, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )[0].codes == ("PS3",)


@pytest.mark.parametrize(
    ("max_cardinality", "expected_codes"),
    [
        # Six supporting criteria are worth six points together and nothing less
        # will do, so a bound of four means there is no sufficient set at all.
        (4, []),
        (6, [("PP1", "PP2", "PP3", "PP4", "PP5", "PP6")]),
    ],
)
def test_a_gap_no_admissible_set_can_close_is_reported_as_no_set_rather_than_a_partial_one(
    max_cardinality, expected_codes
):
    """A set that does not reach the threshold is not a plan, it is a distraction.

    ``max_cardinality`` is a report-legibility bound, so exceeding it has to
    yield silence, never a truncated set the reader would take for sufficient
    evidence.
    """
    spec = _supporting_spec(6, gap=GapConfig(max_cardinality=max_cardinality))
    candidates = _candidates(spec)

    sets = minimum_sufficient_sets(
        candidates, 6, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )

    assert [s.codes for s in sets] == expected_codes


def test_evidence_the_specification_prices_at_zero_is_never_proposed():
    """Acquiring it could not move the total, so it is not a route to anything.

    A point system may legitimately price a strength at zero -- that is what a
    VCEP does when it retires a criterion without deleting it -- and a set built
    from such criteria would read as a plan while being arithmetically inert.

    Expectation corrected while writing this: candidates come back sorted by
    code, not in declaration order, so the pairs are asserted by code.
    """
    zero_supporting = PointSystem(
        pathogenic_points={
            Strength.SUPPORTING: 0,
            Strength.MODERATE: 2,
            Strength.STRONG: 4,
            Strength.VERY_STRONG: 8,
        }
    )
    spec = _spec(
        _criterion("PP3"),
        _criterion("PM1", strength=Strength.MODERATE),
        point_system=zero_supporting,
    )
    candidates = _candidates(spec)

    sets = minimum_sufficient_sets(
        candidates, 2, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )

    assert [(c.code, c.points) for c in candidates] == [("PM1", 2), ("PP3", 0)]
    assert [s.codes for s in sets] == [("PM1",)]


def test_only_the_cheapest_sets_are_kept_when_more_than_the_cap_would_qualify():
    """The cap is what keeps one row of the map readable by a human.

    Three interchangeable strong criteria give three equally valid routes; the
    report keeps the declared number of them, chosen by cost and then by code so
    that two runs of the same input pick the same two.
    """
    spec = _spec(
        _criterion("PS1", strength=Strength.STRONG),
        _criterion("PS2", strength=Strength.STRONG),
        _criterion("PS3", strength=Strength.STRONG),
        gap=GapConfig(max_sets_per_target=2),
    )
    candidates = _candidates(spec)

    sets = minimum_sufficient_sets(
        candidates, 4, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )

    assert [s.codes for s in sets] == [("PS1",), ("PS2",)]


def test_a_set_requiring_the_unobtainable_is_no_set_at_all_unless_the_spec_asks_for_it():
    """An intractable requirement is a research plan nobody can execute.

    Reporting "no route" is the honest answer for a variant whose only route is
    a de novo event in an adult-onset disease; the specification can still ask
    to see those routes explicitly, and then they must appear.
    """
    spec = _spec(_extrinsic("PS2", feasibility=FeasibilityTag.INTRACTABLE))
    permissive = _spec(
        _extrinsic("PS2", feasibility=FeasibilityTag.INTRACTABLE),
        gap=GapConfig(include_intractable=True),
    )

    blocked = minimum_sufficient_sets(
        _candidates(spec), 4, ACMGClass.LIKELY_PATHOGENIC, Direction.PATHOGENIC, spec
    )
    allowed = minimum_sufficient_sets(
        _candidates(permissive),
        4,
        ACMGClass.LIKELY_PATHOGENIC,
        Direction.PATHOGENIC,
        permissive,
    )

    assert blocked == ()
    assert [s.codes for s in allowed] == [("PS2",)]


# --------------------------------------------------------------------------
# Why the variant is stuck.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("priority", "expected"),
    [
        (
            (BlockingReason.MISSING_FUNCTIONAL, BlockingReason.PREDICTOR_INDETERMINATE),
            BlockingReason.MISSING_FUNCTIONAL,
        ),
        # Flipping the specification's list flips the answer, and nothing else
        # about the variant changed: the ranking is data, not code.
        (
            (BlockingReason.PREDICTOR_INDETERMINATE, BlockingReason.MISSING_FUNCTIONAL),
            BlockingReason.PREDICTOR_INDETERMINATE,
        ),
        # A signal the specification never ranked cannot win by accident.
        ((BlockingReason.MISSING_SEGREGATION,), BlockingReason.FREQUENCY_UNINFORMATIVE),
        ((), BlockingReason.FREQUENCY_UNINFORMATIVE),
    ],
)
def test_the_reported_block_is_ranked_by_the_specification_and_never_by_criterion_order(
    priority, expected
):
    """This label is the ``GROUP BY`` the whole map aggregates on.

    Two blocked criteria of different kinds are the normal case, so *which* one
    is reported has to come from the specification's own list. Deriving it from
    declaration order would make a cosmetic edit to a YAML file move thousands
    of variants between buckets.
    """
    spec = _spec(
        _criterion("PP3", blocks_as=BlockingReason.PREDICTOR_INDETERMINATE),
        _criterion("PS3", blocks_as=BlockingReason.MISSING_FUNCTIONAL),
        blocking=BlockingConfig(priority=priority, default=BlockingReason.FREQUENCY_UNINFORMATIVE),
    )
    evaluation = evaluate_variant(_variant(), _context(), spec)

    reason = derive_blocking_reason(evaluation, spec, [], ())

    assert {s.outcome for s in evaluation.skipped} == {CriterionOutcome.NOT_EVALUABLE}
    assert reason is expected


def test_a_criterion_answered_no_names_the_block_only_when_nothing_else_can():
    """Last resort, and it has to stay last.

    A criterion whose rule was answered is not a gap, so it must not outrank a
    criterion that is genuinely waiting for data -- here the not-met one carries
    the higher-priority reason and still loses. When nothing is waiting, though,
    the alternative is a default that names no evidence at all, which tells a
    curator nothing.
    """
    spec = _spec(
        _criterion("PS3", blocks_as=BlockingReason.MISSING_FUNCTIONAL),
        _criterion("PP3", blocks_as=BlockingReason.PREDICTOR_INDETERMINATE),
    )
    answered = {"evidence": {"ps3": False}}
    all_answered = {"evidence": {"ps3": False, "pp3": False}}

    with_a_waiting_criterion = derive_blocking_reason(
        evaluate_variant(_variant(), _context(answered), spec), spec, [], ()
    )
    with_nothing_waiting = derive_blocking_reason(
        evaluate_variant(_variant(), _context(all_answered), spec), spec, [], ()
    )

    assert evaluate_variant(_variant(), _context(all_answered), spec).skipped[0].reason is (
        SkipReason.RULE_NOT_MET
    )
    assert with_a_waiting_criterion is BlockingReason.PREDICTOR_INDETERMINATE
    assert with_nothing_waiting is BlockingReason.MISSING_FUNCTIONAL


def test_a_variant_no_criterion_can_explain_falls_back_to_the_declared_default():
    """Silence is a possible state of the trace, and the column is not nullable.

    Every criterion here was answered and none of them declares a blocking
    reason, so the engine has nothing to report and must still report something.
    The default is the specification's own, which keeps that row aggregable
    instead of dropping it out of the ``GROUP BY`` entirely.
    """
    spec = _spec(
        _criterion("PM1"),
        blocking=BlockingConfig(default=BlockingReason.FREQUENCY_UNINFORMATIVE),
    )
    evaluation = evaluate_variant(_variant(), _context({"evidence": {"pm1": False}}), spec)

    reason = derive_blocking_reason(evaluation, spec, [], ())

    assert evaluation.skipped[0].reason is SkipReason.RULE_NOT_MET
    assert reason is BlockingReason.FREQUENCY_UNINFORMATIVE


def test_contradictory_evidence_outranks_every_missing_piece_of_evidence():
    """A variant pulled both ways is not waiting for data, it is disputed.

    PS3 below is genuinely blocked and carries the highest-priority reason in
    the shipped list, and it still loses. Labelling this variant
    ``MISSING_FUNCTIONAL`` would put it in the queue of gaps an assay could
    close, when what it actually needs is adjudication of evidence already in
    hand.
    """
    spec = _spec(
        _criterion("PM1", strength=Strength.MODERATE),
        _criterion("BS1", direction=Direction.BENIGN, strength=Strength.MODERATE),
        _criterion("PS3", blocks_as=BlockingReason.MISSING_FUNCTIONAL),
    )
    context = _context({"evidence": {"pm1": True, "bs1": True}})

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    assert [c.code for c in analysis.candidates] == ["PS3"]
    assert analysis.blocking_reason is BlockingReason.CONFLICTING_INTRINSIC


def test_a_criterion_the_engine_can_never_evaluate_speaks_through_the_proposed_set():
    """Extrinsic criteria are never "not evaluable" -- they are never evaluated.

    They are skipped as out of the engine's reach, so the only place they can
    name a block is the set the report proposes. Losing that would leave every
    variant whose only route is a case-control study labelled by whatever
    intrinsic criterion happened to be answered no.
    """
    spec = _spec(
        _criterion("PS3", blocks_as=BlockingReason.MISSING_FUNCTIONAL),
        _extrinsic(
            "PS4",
            max_plausible_strength=Strength.VERY_STRONG,
            blocks_as=BlockingReason.MISSING_CASE_CONTROL,
        ),
    )
    context = _context({"evidence": {"ps3": False}})

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    assert [s.codes for s in analysis.minimum_sufficient_sets] == [("PS4",)]
    assert analysis.blocking_reason is BlockingReason.MISSING_CASE_CONTROL


def test_a_variant_already_at_the_benign_threshold_is_offered_no_further_benign_evidence():
    """The two directions are asked independently, and only one is still open.

    At exactly the likely-benign threshold there is nothing left to prove on
    that side, while the pathogenic side is still seven points away. Proposing
    benign evidence here would be work with no verdict attached to it; dropping
    the pathogenic route would hide the only reclassification still possible.
    """
    spec = _spec(
        _criterion("BP4", direction=Direction.BENIGN),
        _criterion("PVS1", strength=Strength.VERY_STRONG),
        _criterion("BS3", direction=Direction.BENIGN, strength=Strength.STRONG),
    )
    context = _context({"evidence": {"bp4": True}})

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    assert (analysis.gap_to_lp, analysis.gap_to_lb) == (7, None)
    assert [(s.target, s.codes) for s in analysis.minimum_sufficient_sets] == [
        (ACMGClass.LIKELY_PATHOGENIC, ("PVS1",))
    ]


def test_a_variant_that_left_uncertainty_is_reported_as_resolved_not_as_blocked():
    """The enum has to distinguish "no longer a question" from "still stuck".

    Counting a resolved variant under any missing-evidence reason would inflate
    the very aggregation the map is read for.
    """
    spec = _spec(_criterion("PVS1", strength=Strength.VERY_STRONG))
    context = _context({"evidence": {"pvs1": True}})

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    assert analysis.class_ceiling_intrinsic is ACMGClass.LIKELY_PATHOGENIC
    assert analysis.blocking_reason is BlockingReason.RESOLVED_NOT_BLOCKED


# --------------------------------------------------------------------------
# The ceiling and the headline report.
# --------------------------------------------------------------------------


def test_the_ceiling_counts_one_member_per_group_and_ignores_what_needs_a_patient():
    """The ceiling answers "how far can this go with nobody in the room".

    Adding a second member of a mutex group would count the same evidence twice,
    and adding an extrinsic criterion would answer a different question
    entirely -- one whose answer is not reachable without a patient.
    """
    spec = _spec(
        _criterion("PS3", strength=Strength.STRONG, mutex_group="functional"),
        _criterion("PM1", strength=Strength.MODERATE, mutex_group="functional"),
        _criterion("PP3"),
        _criterion("BS1", direction=Direction.BENIGN, strength=Strength.STRONG),
        _extrinsic("PS4"),
    )
    context = _context()

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    # 4 for the strongest of the functional group, 1 for the ungrouped PP3.
    assert analysis.points_ceiling_intrinsic == 5
    assert analysis.class_ceiling_intrinsic is ACMGClass.UNCERTAIN


def test_the_actionable_report_lists_only_evidence_that_already_exists_somewhere():
    """Engineering work and bench work must not share a column.

    ``available_uningested`` is read as a work queue: everything in it can be
    closed by wiring in a published dataset. One criterion that actually needs
    an assay or a cohort turns that queue into an estimate nobody can plan
    against. The order is alphabetical rather than declaration order because the
    tuple is written into a Parquet file whose bytes are compared between runs.
    """
    spec = _spec(
        _criterion("PS1", feasibility=FeasibilityTag.AVAILABLE_UNINGESTED),
        _criterion("PS3", feasibility=FeasibilityTag.ASSAY_FEASIBLE),
        _criterion("BA1", direction=Direction.BENIGN, strength=Strength.STAND_ALONE),
        _extrinsic("PS4"),
    )
    context = _context()

    analysis = analyse_gap(evaluate_variant(_variant(), context, spec), spec, context)

    assert [c.code for c in analysis.candidates] == ["BA1", "PS1", "PS3", "PS4"]
    assert analysis.available_uningested == ("BA1", "PS1")
