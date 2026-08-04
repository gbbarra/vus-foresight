"""Evidence signatures: the computational half of spec section 5.

Section 5 promises a double gain from equivalence classes -- a map legible by
region *and* less computation. The legibility half falls out of grouping the
finished rows. The computation half does not: a class derived from the
evaluation trace can only be known after paying for the evaluation.

The signature closes that. It is a canonical digest of every context path the
specification is able to read, plus the consequence terms the criterion gates
use. Two variants with the same signature cannot evaluate differently, because
there is nothing else for the engine to look at -- so the evaluation and the gap
analysis are computed once and reused.

Note what is *not* in the signature: the variant's identity. Two missense
variants with an identical evidence profile share an evaluation but still get
distinct ``equivalence_class_id`` values, because section 5 forbids collapsing
missense in the *report*. Sharing the computation and separating the reporting
are different questions, and conflating them would either lose the saving or
lose the distinction.
"""

from __future__ import annotations

from hashlib import blake2b
from typing import Any, Sequence

from ..variant import Variant
from .context import MISSING, EvidenceContext

__all__ = ["canonical", "evidence_signature", "SignatureCache"]

_MISSING_TOKEN = ("\x00missing",)


def canonical(value: Any) -> Any:
    """Reduce a context value to something hashable and order-independent.

    Lists keep their order -- ``region.tags`` is emitted sorted by its adapter,
    and two different orders would mean two different adapter outputs, which is
    a difference worth separating rather than smoothing over.
    """
    if value is MISSING:
        return _MISSING_TOKEN
    if isinstance(value, dict):
        return tuple((key, canonical(value[key])) for key in sorted(value))
    if isinstance(value, (list, tuple)):
        return tuple(canonical(item) for item in value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def evidence_signature(variant: Variant, context: EvidenceContext, footprint: Sequence[str]) -> str:
    """A digest of everything that can influence this variant's evaluation.

    ``footprint`` must be :attr:`VCEPSpec.field_footprint`. Consequence terms are
    included separately because the ``applies_to`` gate reads them from the
    variant rather than from the context.
    """
    parts: list[Any] = [tuple(term.value for term in variant.consequence_terms)]
    for path in footprint:
        parts.append((path, canonical(context.get(path))))
    return blake2b(repr(tuple(parts)).encode("utf-8"), digest_size=16).hexdigest()


class SignatureCache:
    """Memoises whatever a caller computes per signature, and counts the hits.

    The counters are not instrumentation for its own sake: a run whose hit rate
    collapses has usually gained a per-variant field in its footprint, and that
    is worth noticing at the point it happens rather than in a profiler later.
    """

    __slots__ = ("_entries", "hits", "misses")

    def __init__(self) -> None:
        self._entries: dict[str, Any] = {}
        self.hits = 0
        self.misses = 0

    def get(self, signature: str) -> Any | None:
        entry = self._entries.get(signature)
        if entry is None:
            self.misses += 1
        else:
            self.hits += 1
        return entry

    def put(self, signature: str, value: Any) -> None:
        self._entries[signature] = value

    def clear(self) -> None:
        self._entries.clear()
        self.hits = 0
        self.misses = 0

    @property
    def size(self) -> int:
        return len(self._entries)

    @property
    def lookups(self) -> int:
        return self.hits + self.misses

    @property
    def hit_rate(self) -> float:
        return self.hits / self.lookups if self.lookups else 0.0

    def summary(self) -> str:
        return (
            f"{self.size:,} distinct evidence profiles across {self.lookups:,} "
            f"variants ({self.hit_rate:.1%} reused)"
        )
