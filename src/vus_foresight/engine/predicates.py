"""The rule DSL: a small, closed set of comparison operators.

The set is deliberately small and deliberately not extensible from YAML. A
specification says *what* to compare, never *how*; anything that needs a new
operator needs a code change plus a test, which is the point. A DSL that can
express arbitrary computation is a programming language with no debugger, and
VCEP specifications do not need one.

Absence is never true. Every operator applied to :data:`MISSING` returns
``False`` and the evaluator turns that into ``NOT_EVALUABLE`` -- so a rule can
never accidentally fire because data was not loaded.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from typing import Any

from pydantic import BaseModel, ConfigDict, model_validator

from .context import MISSING, EvidenceContext, LookupLog

__all__ = ["OPERATORS", "Leaf", "Operator", "RuleNode", "evaluate_rule"]


def _as_set(value: Any) -> set[Any]:
    if isinstance(value, (list, tuple, set, frozenset)):
        return {_scalar(v) for v in value}
    return {_scalar(value)}


def _scalar(value: Any) -> Any:
    """Normalise enums to their values so YAML can name them as plain strings."""
    return getattr(value, "value", value)


OPERATORS: dict[str, Callable[[Any, Any], bool]] = {
    "eq": lambda a, b: _scalar(a) == _scalar(b),
    "ne": lambda a, b: _scalar(a) != _scalar(b),
    "lt": lambda a, b: a < b,
    "le": lambda a, b: a <= b,
    "gt": lambda a, b: a > b,
    "ge": lambda a, b: a >= b,
    "in": lambda a, b: _scalar(a) in _as_set(b),
    "not_in": lambda a, b: _scalar(a) not in _as_set(b),
    # Collection-valued fields, e.g. consequence_terms or region tags.
    "contains": lambda a, b: _scalar(b) in _as_set(a),
    "contains_any": lambda a, b: bool(_as_set(a) & _as_set(b)),
    "contains_all": lambda a, b: _as_set(b) <= _as_set(a),
    "is_true": lambda a, b: a is True,
    "is_false": lambda a, b: a is False,
    "between": lambda a, b: b[0] <= a <= b[1],
}

Operator = str


class Leaf(BaseModel):
    """``{field: ..., op: ..., value: ...}``"""

    model_config = ConfigDict(frozen=True, extra="forbid")

    field: str
    op: Operator
    value: Any = None

    @model_validator(mode="after")
    def _known_operator(self) -> Leaf:
        if self.op not in OPERATORS:
            raise ValueError(f"unknown operator {self.op!r}; available: {sorted(OPERATORS)}")
        return self

    def evaluate(self, context: EvidenceContext, log: LookupLog) -> bool:
        observed = context.get(self.field, log)
        if observed is MISSING:
            return False
        try:
            return bool(OPERATORS[self.op](observed, self.value))
        except TypeError:
            # A comparison between incompatible types is a specification bug,
            # but it must not be able to make a rule fire.
            raise ValueError(
                f"cannot apply {self.op!r} to {self.field}={observed!r} and value={self.value!r}"
            ) from None


class RuleNode(BaseModel):
    """A boolean combination of leaves. Exactly one field may be set."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    all: tuple[RuleNode, ...] | None = None
    any: tuple[RuleNode, ...] | None = None
    none: tuple[RuleNode, ...] | None = None
    leaf: Leaf | None = None

    @model_validator(mode="before")
    @classmethod
    def _accept_bare_leaf(cls, value: Any) -> Any:
        """Allow ``{field: x, op: y, value: z}`` without a ``leaf:`` wrapper."""
        if isinstance(value, dict) and "field" in value and "op" in value:
            return {"leaf": value}
        return value

    @model_validator(mode="after")
    def _exactly_one(self) -> RuleNode:
        set_fields = [
            name for name in ("all", "any", "none", "leaf") if getattr(self, name) is not None
        ]
        if len(set_fields) != 1:
            raise ValueError(
                f"a rule node needs exactly one of all/any/none/leaf, got {set_fields}"
            )
        return self

    def iter_fields(self) -> Iterator[str]:
        """Every context path this subtree can read.

        Used to compute a specification's field footprint, which is what makes
        class-level evaluation sound: two variants agreeing on every path the
        spec can read must evaluate identically.
        """
        if self.leaf is not None:
            yield self.leaf.field
            return
        for child in self.all or self.any or self.none or ():
            yield from child.iter_fields()

    def evaluate(self, context: EvidenceContext, log: LookupLog) -> bool:
        if self.leaf is not None:
            return self.leaf.evaluate(context, log)
        if self.all is not None:
            # Every branch is evaluated, not short-circuited: the missing-field
            # log has to be complete for the gap report to be complete.
            return all([node.evaluate(context, log) for node in self.all])
        if self.any is not None:
            return any([node.evaluate(context, log) for node in self.any])
        assert self.none is not None
        return not any([node.evaluate(context, log) for node in self.none])


RuleNode.model_rebuild()


def evaluate_rule(rule: RuleNode | None, context: EvidenceContext) -> tuple[bool, LookupLog]:
    """Evaluate a rule and return the verdict together with its lookup log."""
    log = LookupLog()
    if rule is None:
        return False, log
    return rule.evaluate(context, log), log
