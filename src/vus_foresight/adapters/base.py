"""Data source adapters (spec section 8).

Every external source enters through one interface, declares a version, and is
read from a local snapshot. **No adapter performs a network call at engine
runtime.** Reproducibility is the reason: "same input + same spec version + same
database snapshot = same output, always" is not achievable if a source can move
under a running job.

An adapter that has nothing to say about a variant returns an empty mapping. It
must never return a placeholder value -- absent and zero are different evidence
states, and :mod:`vus_foresight.engine.context` keeps them different.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import Any, ClassVar

from ..genome.transcript import Transcript
from ..variant import Variant

__all__ = ["Adapter", "AdapterRegistry", "NullAdapter"]


class Adapter(ABC):
    """One namespace of the evidence context, from one versioned source."""

    #: Top-level context namespace this adapter owns, e.g. ``frequency``.
    namespace: ClassVar[str]
    #: Short source name, e.g. ``gnomad``.
    name: ClassVar[str] = "unnamed"

    def __init__(self, version: str) -> None:
        if not version:
            raise ValueError(
                f"{type(self).__name__} must declare a version; an undated source "
                "cannot appear in a reproducible provenance record"
            )
        self.version = version

    @property
    def source_id(self) -> str:
        """``name@version``, copied into every criterion this adapter feeds."""
        return f"{self.name}@{self.version}"

    @abstractmethod
    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        """Evidence for one variant. Empty mapping when nothing is known."""

    def describe(self) -> dict[str, str]:
        return {"namespace": self.namespace, "source": self.source_id}


class NullAdapter(Adapter):
    """An adapter that knows nothing, on purpose.

    Used to run the pipeline with a source deliberately switched off, so the map
    shows what the gap looks like *without* that source. That is the cheapest
    way to answer "how much would ingesting this dataset buy us".
    """

    def __init__(self, namespace: str, *, name: str = "null", version: str = "0") -> None:
        super().__init__(version)
        self.namespace = namespace  # type: ignore[misc]
        self.name = name  # type: ignore[misc]

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        return {}


class AdapterRegistry:
    """An ordered, namespace-unique collection of adapters."""

    def __init__(self, adapters: Iterable[Adapter] = ()) -> None:
        self._adapters: list[Adapter] = []
        for adapter in adapters:
            self.add(adapter)

    def add(self, adapter: Adapter) -> None:
        existing = {a.namespace for a in self._adapters}
        if adapter.namespace in existing:
            raise ValueError(
                f"namespace {adapter.namespace!r} is already claimed; two adapters "
                "writing one namespace would make provenance ambiguous"
            )
        self._adapters.append(adapter)

    def __iter__(self):
        return iter(self._adapters)

    def __len__(self) -> int:
        return len(self._adapters)

    @property
    def versions(self) -> dict[str, str]:
        """``namespace -> source@version``, for the provenance block of each row."""
        return {a.namespace: a.source_id for a in self._adapters}
