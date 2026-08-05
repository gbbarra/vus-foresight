"""The closed predicate DSL, which is what a specification is allowed to say.

Local to one module, so it lives here rather than in a layer. The DSL is the
boundary between rule-as-data and code: everything a VCEP specification can
express passes through it, and everything it cannot express is a rule that
cannot be written. That makes its refusals part of the contract -- a
specification containing a typo must fail loudly at load or evaluation, never
quietly evaluate to false, because a criterion that silently never fires
removes points from every variant in the gene and looks exactly like "the
evidence is absent".
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from vus_foresight.engine.context import MISSING, EvidenceContext
from vus_foresight.engine.predicates import OPERATORS, Leaf, RuleNode, evaluate_rule


def _context(**values):
    return EvidenceContext(values)


# --------------------------------------------------------------------------
# What the DSL refuses to accept at all.
# --------------------------------------------------------------------------


def test_an_unknown_operator_is_rejected_when_the_rule_is_built():
    """At load, not at evaluation: a typo must not survive to run time."""
    with pytest.raises(ValidationError, match="unknown operator"):
        Leaf(field="gnomad.faf95_popmax", op="greater_than", value=0.001)


def test_the_error_lists_the_operators_that_do_exist():
    with pytest.raises(ValidationError) as exc:
        Leaf(field="x", op="=~", value=1)
    for name in ("eq", "gt", "lt"):
        assert name in str(exc.value)


@pytest.mark.parametrize(
    "node",
    [
        {},
        {"all": [{"leaf": {"field": "x", "op": "eq", "value": 1}}], "any": []},
    ],
)
def test_a_node_needs_exactly_one_of_all_any_none_leaf(node):
    with pytest.raises(ValidationError, match="exactly one"):
        RuleNode(**node)


def test_an_unknown_key_in_a_node_is_refused_rather_than_ignored():
    """extra=forbid: a misspelt 'nne' must not silently become no constraint."""
    with pytest.raises(ValidationError):
        RuleNode(nne=[{"leaf": {"field": "x", "op": "eq", "value": 1}}])


# --------------------------------------------------------------------------
# MISSING is not False.
# --------------------------------------------------------------------------


def test_a_missing_field_makes_a_leaf_false_without_raising():
    leaf = Leaf(field="gnomad.faf95_popmax", op="gt", value=0.001)
    verdict, log = evaluate_rule(RuleNode(leaf=leaf), _context())

    assert verdict is False
    assert "gnomad.faf95_popmax" in log.missing


def test_a_missing_field_is_recorded_even_inside_a_negation():
    """``none`` inverts the verdict; it must not invert the bookkeeping.

    A criterion whose absent field satisfies a ``none`` clause would otherwise
    apply on evidence nobody has, and the gap report would show nothing
    missing.
    """
    rule = RuleNode(none=[RuleNode(leaf=Leaf(field="region.tags", op="contains", value="repeat"))])
    verdict, log = evaluate_rule(rule, _context())

    assert verdict is True
    assert "region.tags" in log.missing


def test_comparing_incompatible_types_raises_rather_than_returning_false():
    """A specification bug must not be able to make a rule quietly not fire."""
    leaf = Leaf(field="gnomad.faf95_popmax", op="gt", value=0.001)
    # Nested, because the field is a dotted *path*: a literal key containing a
    # dot resolves to nothing and would leave this asserting on MISSING.
    context = _context(gnomad={"faf95_popmax": "not a number"})
    with pytest.raises(ValueError, match="cannot apply"):
        evaluate_rule(RuleNode(leaf=leaf), context)


# --------------------------------------------------------------------------
# The three combinators, including the one no shipped rule uses.
# --------------------------------------------------------------------------


def test_any_is_true_when_one_branch_is_true():
    rule = RuleNode(
        any=[
            RuleNode(leaf=Leaf(field="a", op="eq", value=1)),
            RuleNode(leaf=Leaf(field="b", op="eq", value=2)),
        ]
    )
    assert evaluate_rule(rule, _context(a=99, b=2))[0] is True
    assert evaluate_rule(rule, _context(a=99, b=99))[0] is False


def test_every_branch_is_evaluated_even_once_the_verdict_is_settled():
    """Not an optimisation left on the table -- the reason the gap map works.

    ``blocking_reason`` is derived from which fields were missing, so the
    lookup log has to be complete. Short-circuiting would make the recorded
    cause depend on the order the branches happen to sit in the YAML.
    """
    rule = RuleNode(
        all=[
            RuleNode(leaf=Leaf(field="known", op="eq", value=0)),
            RuleNode(leaf=Leaf(field="absent.one", op="eq", value=1)),
            RuleNode(leaf=Leaf(field="absent.two", op="eq", value=2)),
        ]
    )
    verdict, log = evaluate_rule(rule, _context(known=999))

    assert verdict is False
    assert {"absent.one", "absent.two"} <= set(log.missing)

    rule_any = RuleNode(
        any=[
            RuleNode(leaf=Leaf(field="satisfied", op="eq", value=1)),
            RuleNode(leaf=Leaf(field="absent.three", op="eq", value=1)),
        ]
    )
    verdict_any, log_any = evaluate_rule(rule_any, _context(satisfied=1))
    assert verdict_any is True
    assert "absent.three" in log_any.missing


def test_a_criterion_with_no_rule_never_applies():
    """``None`` means "this criterion has no machine-checkable condition"."""
    verdict, log = evaluate_rule(None, _context(anything=1))

    assert verdict is False
    assert log.missing == []


# --------------------------------------------------------------------------
# The operators themselves, at their boundaries.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("op", "observed", "value", "expected"),
    [
        ("gt", 0.001, 0.001, False),
        ("gt", 0.0010001, 0.001, True),
        ("ge", 0.001, 0.001, True),
        ("lt", 0.001, 0.001, False),
        ("le", 0.001, 0.001, True),
        ("eq", 0, 0, True),
        ("ne", 0, 0, False),
        ("in", "missense", ["missense", "nonsense"], True),
        ("in", "synonymous", ["missense", "nonsense"], False),
        ("not_in", "synonymous", ["missense"], True),
        ("contains", ["critical", "domain"], "critical", True),
        ("contains", ["domain"], "critical", False),
        ("contains_any", ["domain", "repeat"], ["repeat"], True),
        ("contains_any", ["domain"], ["repeat"], False),
        ("contains_all", ["critical", "domain"], ["critical", "domain"], True),
        ("contains_all", ["critical"], ["critical", "domain"], False),
        ("between", 5, [1, 10], True),
        ("between", 1, [1, 10], True),
        ("between", 10, [1, 10], True),
        ("between", 11, [1, 10], False),
        ("is_true", True, None, True),
        ("is_false", False, None, True),
    ],
)
def test_each_operator_at_its_boundary(op, observed, value, expected):
    """Every numeric comparison is pinned at equality, where off-by-one lives."""
    # No skip guard here on purpose. The first version had one, and it turned
    # two wrong operator names into two silent passes -- which is exactly the
    # decorative test the standard forbids.
    leaf = Leaf(field="x", op=op, value=value)
    assert evaluate_rule(RuleNode(leaf=leaf), _context(x=observed))[0] is expected


def test_the_operator_set_is_closed():
    """Adding one is a deliberate act: it widens what a specification may say.

    Pinned by name so that a new operator arrives with a test rather than
    arriving because somebody needed it once.
    """
    assert sorted(OPERATORS) == [
        "between",
        "contains",
        "contains_all",
        "contains_any",
        "eq",
        "ge",
        "gt",
        "in",
        "is_false",
        "is_true",
        "le",
        "lt",
        "ne",
        "not_in",
    ]


def test_a_published_null_counts_as_missing_rather_than_as_a_value():
    """Characterisation, and the expectation was what changed.

    This was written expecting a published ``None`` to be distinguishable from
    an absent key -- "consulted, no result" against "nobody looked". It is not:
    ``EvidenceContext.get`` coerces None to MISSING deliberately, and that is
    the right call. An adapter reporting None means the evidence is not
    available for this variant, which is precisely what the gap report should
    say. The alternative is worse: a present None would let ``ne`` fire, so a
    criterion could apply on the strength of a source having nothing to say.
    """
    leaf = Leaf(field="functional.classification", op="ne", value="abnormal")
    verdict, log = evaluate_rule(RuleNode(leaf=leaf), _context(functional={"classification": None}))

    assert verdict is False
    assert log.missing == ["functional.classification"]
    assert MISSING is not None
