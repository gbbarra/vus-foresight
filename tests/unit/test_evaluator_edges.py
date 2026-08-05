"""The edges of criterion evaluation: the trace, not the verdict.

Oracle for this file: the evaluator's own stated contract, exercised against
specifications built inside the test rather than against ``config/specs``. Every
assertion names a concrete outcome, reason, detail string, strength or point
value, so a shipped threshold moving cannot make these tests lie -- and a change
to what the engine *means* by NOT_MET, NOT_EVALUABLE or NOT_APPLICABLE cannot
pass unnoticed.

How these tests could nevertheless be wrong: by describing a configuration no
YAML file could produce. The guard is that everything here is constructed
through the very Pydantic models ``load_spec`` validates into, with their
validators running, so anything a test builds is something a specification could
legitimately say. What is *not* guarded is realism -- ``config/specs/toy`` gives
every criterion a ``requires`` list covering its own rule, so several paths
below are unreachable from the shipped files and reachable only from a
specification that omits one. That is precisely why they are tested here: an
untested fallback is a wrong number waiting for the first spec that needs it.
"""

from __future__ import annotations

import pytest

from vus_foresight.acmg import (
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    SkipReason,
    Strength,
)
from vus_foresight.engine.context import EvidenceContext
from vus_foresight.engine.evaluator import Evaluation, evaluate_variant, render_evidence
from vus_foresight.engine.predicates import Leaf, RuleNode
from vus_foresight.engine.spec import CriterionSpec, StrengthRung, VCEPSpec
from vus_foresight.gapmap import AppliedCriterion
from vus_foresight.variant import Consequence, Variant, VariantKind

# --------------------------------------------------------------------------
# Builders. Deliberately thin: they set defaults, never behaviour.
# --------------------------------------------------------------------------


def _leaf(field: str, op: str, value=None) -> RuleNode:
    return RuleNode(leaf=Leaf(field=field, op=op, value=value))


#: A rule that fires whenever ``evidence.flag`` is present and true.
_FLAG_IS_TRUE = _leaf("evidence.flag", "is_true")

#: The context the flag rule expects. Nothing else lives in it, so any other
#: path a test names is genuinely absent rather than accidentally present.
_FLAG_SET = {"evidence": {"flag": True}}


def _criterion(code: str, **overrides) -> CriterionSpec:
    fields = {
        "code": code,
        "direction": Direction.PATHOGENIC,
        "evidence_class": EvidenceClass.INTRINSIC,
        "rule": _FLAG_IS_TRUE,
    }
    fields.update(overrides)
    return CriterionSpec(**fields)


def _spec(*criteria: CriterionSpec) -> VCEPSpec:
    return VCEPSpec(
        spec_id="edges",
        name="Specification built for the evaluator edge tests",
        version="0.0.0",
        verified=True,
        criteria=criteria,
    )


def _variant(
    consequence: Consequence = Consequence.MISSENSE,
    terms: tuple[Consequence, ...] = (),
) -> Variant:
    return Variant(
        gene="TOY1",
        transcript="NM_999003.1",
        kind=VariantKind.SNV,
        hgvs_c="c.10A>G",
        consequence=consequence,
        consequence_terms=terms,
    )


def _context(data: dict | None = None, sources: dict | None = None) -> EvidenceContext:
    return EvidenceContext(data if data is not None else {}, sources or {})


def _applied(code: str, direction: Direction, points: int) -> AppliedCriterion:
    return AppliedCriterion(
        code=code,
        direction=direction,
        strength=Strength.SUPPORTING,
        points=points,
        evidence_class=EvidenceClass.INTRINSIC,
        evidence="built by the test",
        source="test",
    )


# --------------------------------------------------------------------------
# Querying the trace. Principle 1: the trace is the map, so every criterion
# the specification names has to be answerable by code.
# --------------------------------------------------------------------------


@pytest.fixture
def mixed_trace() -> Evaluation:
    """One evaluation containing all four outcomes at once."""
    spec = _spec(
        _criterion("PM1"),
        _criterion("BA1", direction=Direction.BENIGN, rule=_leaf("evidence.flag", "is_false")),
        _criterion(
            "PS3",
            requires=("functional.classification",),
            rule=_leaf("functional.classification", "eq", "abnormal"),
        ),
        _criterion(
            "PS4",
            evidence_class=EvidenceClass.EXTRINSIC,
            rule=None,
            feasibility=FeasibilityTag.NEEDS_CASES,
        ),
    )

    return evaluate_variant(_variant(), _context(_FLAG_SET), spec)


@pytest.mark.parametrize(
    ("code", "expected"),
    [
        ("PM1", CriterionOutcome.APPLIED),
        ("BA1", CriterionOutcome.NOT_MET),
        ("PS3", CriterionOutcome.NOT_EVALUABLE),
        ("PS4", CriterionOutcome.NOT_APPLICABLE),
        ("PP1", None),
    ],
)
def test_the_trace_answers_for_every_code_and_admits_when_it_has_no_answer(
    mixed_trace, code, expected
):
    """``None`` for an unknown code must not be confused with NOT_MET.

    A caller that reads a missing criterion as "tested and rejected" would
    report a gap the specification never even asked about.
    """
    assert mixed_trace.outcome_of(code) is expected


def test_only_criteria_stopped_by_absent_data_are_indexed_as_gaps(mixed_trace):
    """``missing_by_code`` becomes ``available_uningested``.

    A criterion whose rule was answered "no", or one that is out of scope, is
    not a gap: acquiring more data would change nothing. Letting either leak in
    would inflate the single number the whole project exists to report.
    """
    assert mixed_trace.missing_by_code == {"PS3": ("functional.classification",)}


# --------------------------------------------------------------------------
# NOT_MET against NOT_EVALUABLE -- absence is not falsehood.
# --------------------------------------------------------------------------


def test_a_rule_blind_to_an_undeclared_field_is_not_evaluable_rather_than_unmet():
    """``requires`` is optional, so the rule's own lookup log is the second net.

    Every criterion in ``config/specs/toy`` declares the fields its rule reads,
    which means the prerequisite check catches absence first and this branch
    never runs from a shipped file. A specification that omits ``requires`` is
    still valid, and there the difference is the whole product: NOT_MET is
    excluded from the candidate sets in ``engine/gap.py``, NOT_EVALUABLE is not.
    """
    spec = _spec(_criterion("PP3", rule=_leaf("predictor.bayesdel", "ge", 0.30)))

    evaluation = evaluate_variant(_variant(), _context(_FLAG_SET), spec)

    skipped = evaluation.skipped[0]
    assert (skipped.code, skipped.outcome, skipped.reason) == (
        "PP3",
        CriterionOutcome.NOT_EVALUABLE,
        SkipReason.DATA_MISSING,
    )
    assert skipped.missing_fields == ("predictor.bayesdel",)
    assert evaluation.missing_by_code == {"PP3": ("predictor.bayesdel",)}


def test_a_rule_that_already_fired_is_applied_even_though_one_branch_was_blind():
    """Answerable beats complete: an ``any`` needs only one satisfied branch.

    The absent branch is still recorded in the lookup log, so the guard is
    ``not fired and log.has_missing`` rather than ``log.has_missing`` alone.
    Reporting the untouched path as a gap here would advertise data whose
    acquisition could not change the outcome by a single point.
    """
    rule = RuleNode(any=(_FLAG_IS_TRUE, _leaf("predictor.bayesdel", "ge", 0.30)))
    spec = _spec(_criterion("PP3", rule=rule))

    evaluation = evaluate_variant(_variant(), _context(_FLAG_SET), spec)

    assert evaluation.applied_codes == ("PP3",)
    assert evaluation.missing_by_code == {}
    assert evaluation.points == 1


# --------------------------------------------------------------------------
# The structured reason attached to a blocked criterion.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("code", "data", "expected_detail"),
    [
        (
            "PVS1",
            {"pvs1": {"undetermined_reason": "no splice prediction ingested"}},
            "no splice prediction ingested",
        ),
        # Specifications name modulated strengths as PVS1_Strong, PVS1_Moderate:
        # the test is a prefix, and dropping the suffixed forms would silently
        # blank the explanation for exactly the hardest PVS1 calls.
        (
            "PVS1_Moderate",
            {"pvs1": {"undetermined_reason": "NMD status unknown"}},
            "NMD status unknown",
        ),
        # The tree ran and had nothing structured to say.
        ("PVS1", {}, None),
        # A published null is absence, not an explanation reading "None".
        ("PVS1", {"pvs1": {"undetermined_reason": None}}, None),
        # No other criterion may borrow PVS1's reason as its own.
        ("PM2_Supporting", {"pvs1": {"undetermined_reason": "NMD status unknown"}}, None),
    ],
)
def test_a_blocked_criterion_carries_the_pvs1_reason_and_no_other_criterion_does(
    code, data, expected_detail
):
    spec = _spec(
        _criterion(code, requires=("splice.ds_max",), rule=_leaf("splice.ds_max", "ge", 0.5))
    )

    evaluation = evaluate_variant(_variant(), _context(data), spec)

    skipped = evaluation.skipped[0]
    assert skipped.detail == expected_detail
    assert skipped.missing_fields == ("splice.ds_max",)


# --------------------------------------------------------------------------
# The consequence gate.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("applies_to", "consequence", "terms", "expected"),
    [
        # No gate declared: the criterion is offered to every consequence.
        ((), Consequence.SYNONYMOUS, (), CriterionOutcome.APPLIED),
        # No secondary terms: the primary consequence stands in for the set.
        ((Consequence.MISSENSE,), Consequence.MISSENSE, (), CriterionOutcome.APPLIED),
        ((Consequence.MISSENSE,), Consequence.SYNONYMOUS, (), CriterionOutcome.NOT_APPLICABLE),
        # A secondary term opens the gate: a missense sitting in the splice
        # region is in scope for a splice criterion even though it is not its
        # primary consequence.
        (
            (Consequence.SPLICE_REGION,),
            Consequence.MISSENSE,
            (Consequence.MISSENSE, Consequence.SPLICE_REGION),
            CriterionOutcome.APPLIED,
        ),
        # ... and without that term the same variant is out of scope, so the
        # gate reads the annotation rather than guessing from the position.
        ((Consequence.SPLICE_REGION,), Consequence.MISSENSE, (), CriterionOutcome.NOT_APPLICABLE),
    ],
)
def test_the_consequence_gate_reads_every_term_and_falls_back_to_the_primary_one(
    applies_to, consequence, terms, expected
):
    spec = _spec(_criterion("PM1", applies_to=applies_to))

    evaluation = evaluate_variant(_variant(consequence, terms), _context(_FLAG_SET), spec)

    assert evaluation.outcome_of("PM1") is expected


def test_an_out_of_scope_skip_names_the_declared_scope_and_the_variant():
    """The detail is what a curator reads when asking why nothing was tried.

    It is built from the spec's own ``applies_to``, in declaration order, so it
    stays a byte-for-byte function of the specification file.
    """
    spec = _spec(
        _criterion(
            "PM4",
            applies_to=(Consequence.INFRAME_DELETION, Consequence.INFRAME_INSERTION),
        )
    )

    evaluation = evaluate_variant(_variant(Consequence.MISSENSE), _context(_FLAG_SET), spec)

    skipped = evaluation.skipped[0]
    assert skipped.detail == "applies to inframe_deletion, inframe_insertion; variant is missense"
    assert skipped.reason is SkipReason.CONSEQUENCE_OUT_OF_SCOPE


def test_the_consequence_gate_is_checked_before_the_extrinsic_gate():
    """An out-of-scope extrinsic criterion is out of scope, not merely extrinsic.

    ``engine/gap.py`` turns EXTRINSIC_EVIDENCE_REQUIRED into a candidate for the
    minimum sufficient sets and CONSEQUENCE_OUT_OF_SCOPE into nothing. The order
    of these two guards therefore decides whether the report proposes recruiting
    families to establish a criterion the specification says cannot apply to
    this consequence at all.
    """
    spec = _spec(
        _criterion(
            "BP2",
            direction=Direction.BENIGN,
            evidence_class=EvidenceClass.EXTRINSIC,
            rule=None,
            feasibility=FeasibilityTag.NEEDS_CASES,
            applies_to=(Consequence.NONSENSE,),
        )
    )

    evaluation = evaluate_variant(_variant(Consequence.MISSENSE), _context(_FLAG_SET), spec)

    assert evaluation.skipped[0].reason is SkipReason.CONSEQUENCE_OUT_OF_SCOPE


# --------------------------------------------------------------------------
# The strength ladder.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("bayesdel", "expected_strength", "expected_points"),
    [
        # Rule satisfied, no rung reached: the declared strength stands.
        (0.35, Strength.SUPPORTING, 1),
        (0.60, Strength.MODERATE, 2),
        # Both rungs match. The first one wins even though the second is
        # stronger -- a ladder is an ordered decision list, not a maximum, and
        # the obligation to write it strongest-first is the spec author's.
        (0.95, Strength.MODERATE, 2),
    ],
)
def test_the_first_matching_rung_wins_even_when_a_later_rung_is_stronger(
    bayesdel, expected_strength, expected_points
):
    spec = _spec(
        _criterion(
            "PP3",
            requires=("predictor.bayesdel",),
            rule=_leaf("predictor.bayesdel", "ge", 0.30),
            strength_ladder=(
                StrengthRung(
                    strength=Strength.MODERATE,
                    when=_leaf("predictor.bayesdel", "ge", 0.50),
                ),
                StrengthRung(
                    strength=Strength.STRONG,
                    when=_leaf("predictor.bayesdel", "ge", 0.90),
                ),
            ),
        )
    )

    evaluation = evaluate_variant(_variant(), _context({"predictor": {"bayesdel": bayesdel}}), spec)

    applied = evaluation.applied[0]
    assert applied.strength is expected_strength
    assert applied.points == expected_points


# --------------------------------------------------------------------------
# Provenance stamped onto an applied criterion.
# --------------------------------------------------------------------------


def test_the_source_names_every_adapter_consulted_sorted_and_deduplicated():
    """Provenance has to cover the ladder too, not just the rule.

    The strength a criterion contributes is decided by the ladder, so an adapter
    that only the ladder read is still an input to the point total. The string
    is sorted and deduplicated because it is written verbatim into the Parquet
    output, where byte-for-byte reproducibility is a guarantee.
    """
    spec = _spec(
        _criterion(
            "PS3",
            rule=RuleNode(
                all=(
                    _leaf("predictor.bayesdel", "ge", 0.30),
                    _leaf("frequency.gnomad.faf95_popmax", "le", 0.001),
                    _leaf("splice.ds_max", "ge", 0.50),
                )
            ),
            strength_ladder=(
                StrengthRung(
                    strength=Strength.STRONG,
                    when=_leaf("functional.classification", "eq", "abnormal"),
                ),
            ),
        )
    )
    context = _context(
        {
            "predictor": {"bayesdel": 0.90},
            "frequency": {"gnomad": {"faf95_popmax": 0.0}},
            "splice": {"ds_max": 0.80},
            "functional": {"classification": "abnormal"},
        },
        # "splice" and "frequency" deliberately share one release, and the
        # alphabetical order of the sources is not the order of the namespaces.
        {
            "predictor": "zeta@1",
            "frequency": "alpha@2",
            "splice": "alpha@2",
            "functional": "mid@3",
        },
    )

    evaluation = evaluate_variant(_variant(), context, spec)

    applied = evaluation.applied[0]
    assert applied.source == "alpha@2, mid@3, zeta@1"
    assert applied.strength is Strength.STRONG


def test_a_rule_that_consults_no_adapter_is_attributed_to_the_engine():
    """An empty ``all`` is vacuously true and reads nothing.

    ``source`` is a non-optional column of every applied criterion, so the
    degenerate case still has to name something. "engine" says the point came
    from the specification alone, which is honest; an empty string would be a
    hole in the schema that no downstream consumer could interpret.
    """
    spec = _spec(_criterion("PM4", rule=RuleNode(all=())))

    evaluation = evaluate_variant(_variant(), _context(), spec)

    assert evaluation.applied[0].source == "engine"


# --------------------------------------------------------------------------
# Mutually exclusive groups.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("pathogenic_first", [True, False])
def test_a_mutex_tie_is_broken_by_code_and_not_by_declaration_order(pathogenic_first):
    """A tie here flips the sign of the variant's entire score.

    PP3 and BP4 share the ``computational`` group and are both supporting, so
    overlapping thresholds make them tie at one point each in opposite
    directions. Whichever survives decides whether the variant scores +1 or -1.
    Resolving by list position would make that depend on the order the criteria
    happen to sit in the YAML, which byte-for-byte determinism forbids.
    """
    pp3 = _criterion(
        "PP3",
        mutex_group="computational",
        requires=("predictor.bayesdel",),
        rule=_leaf("predictor.bayesdel", "ge", 0.30),
    )
    bp4 = _criterion(
        "BP4",
        direction=Direction.BENIGN,
        mutex_group="computational",
        requires=("predictor.bayesdel",),
        rule=_leaf("predictor.bayesdel", "le", 0.90),
    )
    spec = _spec(*((pp3, bp4) if pathogenic_first else (bp4, pp3)))

    evaluation = evaluate_variant(_variant(), _context({"predictor": {"bayesdel": 0.50}}), spec)

    assert evaluation.applied_codes == ("PP3",)
    assert evaluation.points == 1
    superseded = evaluation.skipped[0]
    assert superseded.reason is SkipReason.SUPERSEDED_BY_MUTEX
    assert superseded.detail == "superseded by PP3 in mutex group 'computational'"


# --------------------------------------------------------------------------
# Directional conflict.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("pathogenic_points", "benign_points", "min_each", "expected"),
    [
        # One supporting criterion each way is ordinary curation at the shipped
        # threshold of 2, and a conflict only if a spec lowers it to 1.
        (1, -1, 1, True),
        (1, -1, 2, False),
        # The comparison is >=, so the threshold itself conflicts.
        (2, -2, 2, True),
        (2, -1, 2, False),
        # Evidence in one direction only is never a conflict, however large.
        (8, None, 1, False),
        (None, -8, 1, False),
    ],
)
def test_opposing_evidence_is_a_conflict_only_at_or_above_the_declared_threshold(
    pathogenic_points, benign_points, min_each, expected
):
    applied = []
    if pathogenic_points is not None:
        applied.append(_applied("PP3", Direction.PATHOGENIC, pathogenic_points))
    if benign_points is not None:
        applied.append(_applied("BP4", Direction.BENIGN, benign_points))
    evaluation = Evaluation(applied=applied)

    assert evaluation.has_directional_conflict(min_each) is expected


# --------------------------------------------------------------------------
# The evidence string.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("template", "expected"),
    [
        (None, "the fallback"),
        ("no tokens here at all", "no tokens here at all"),
        ("BayesDel {predictor.bayesdel} >= 0.30", "BayesDel 0.42 >= 0.30"),
        # A typo in a template must not abort a map-wide run.
        ("{predictor.absent}", "?"),
        ("{nowhere.at.all} and {predictor.bayesdel}", "? and 0.42"),
        # A published null is absence, and renders as absence.
        ("{predictor.published_null}", "?"),
    ],
)
def test_an_evidence_template_marks_what_it_could_not_resolve_instead_of_raising(
    template, expected
):
    context = _context({"predictor": {"bayesdel": 0.42, "published_null": None}})

    assert render_evidence(template, context, fallback="the fallback") == expected


@pytest.mark.parametrize(
    ("description", "expected"),
    [
        (
            "Computational evidence supports a deleterious effect",
            "Computational evidence supports a deleterious effect",
        ),
        (None, "PP3 rule satisfied"),
    ],
)
def test_a_criterion_without_a_template_falls_back_to_its_description_then_its_code(
    description, expected
):
    """The engine never writes prose of its own beyond this last resort.

    Everything a curator reads in the ``evidence`` column has to be traceable to
    the specification file; the bare code is the honest minimum when the file
    supplies nothing.
    """
    spec = _spec(_criterion("PP3", description=description))

    evaluation = evaluate_variant(_variant(), _context(_FLAG_SET), spec)

    assert evaluation.applied[0].evidence == expected
