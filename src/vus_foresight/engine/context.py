"""The evidence context: everything known about one variant, at one moment.

A context is a nested mapping addressed by dotted path (``frequency.gnomad.af``).
Two things distinguish it from a plain dict and both matter:

* it separates *absent* from *false*. "No population frequency data" and
  "frequency is zero" are completely different evidence states, and collapsing
  them is how a system silently invents evidence;
* every lookup a rule performs is recorded, so a criterion that could not be
  evaluated can say exactly which fields were missing. That record is what turns
  into ``available_uningested`` in the report.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterator

__all__ = ["MISSING", "Missing", "EvidenceContext", "LookupLog"]


class Missing:
    """Sentinel for an absent path. Distinct from ``None`` and from ``False``."""

    _instance: "Missing | None" = None

    def __new__(cls) -> "Missing":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:  # pragma: no cover - debug aid
        return "MISSING"

    def __bool__(self) -> bool:
        return False


MISSING = Missing()


@dataclass(slots=True)
class LookupLog:
    """Which paths a single criterion evaluation touched, and which were absent."""

    touched: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)

    def record(self, path: str, present: bool) -> None:
        if path not in self.touched:
            self.touched.append(path)
        if not present and path not in self.missing:
            self.missing.append(path)

    @property
    def has_missing(self) -> bool:
        return bool(self.missing)


@dataclass(slots=True)
class EvidenceContext:
    """Nested evidence for one variant, with provenance per top-level namespace."""

    data: dict[str, Any] = field(default_factory=dict)
    #: namespace -> "adapter@version", copied into every applied criterion so a
    #: row can be traced back to the exact dataset that produced it.
    sources: dict[str, str] = field(default_factory=dict)

    def merge(self, namespace: str, payload: dict[str, Any], source: str) -> None:
        """Attach one adapter's output under its own namespace."""
        if namespace in self.data:
            raise ValueError(f"namespace {namespace!r} is already populated")
        self.data[namespace] = payload
        self.sources[namespace] = source

    def get(self, path: str, log: LookupLog | None = None) -> Any:
        """Resolve a dotted path, returning :data:`MISSING` if any segment is absent."""
        node: Any = self.data
        for segment in path.split("."):
            if isinstance(node, dict) and segment in node:
                node = node[segment]
            else:
                node = MISSING
                break
        if node is None:
            node = MISSING
        if log is not None:
            log.record(path, node is not MISSING)
        return node

    def source_for(self, path: str) -> str:
        """The adapter that supplied the namespace a path lives in."""
        namespace = path.split(".", 1)[0]
        return self.sources.get(namespace, "unknown")

    def namespaces(self) -> Iterator[str]:
        yield from sorted(self.data)

    def flat(self) -> dict[str, Any]:
        """Flatten to dotted keys. Used for evidence-string rendering only."""
        out: dict[str, Any] = {}

        def walk(prefix: str, node: Any) -> None:
            if isinstance(node, dict):
                for key, value in node.items():
                    walk(f"{prefix}.{key}" if prefix else str(key), value)
            else:
                out[prefix] = node

        walk("", self.data)
        return out
