"""Inspect what has actually been materialised, before computing anything.

Acquiring the external sources (``docs/data-acquisition.md``) is a sequence of
large downloads and region slices, and the failure modes are quiet: a snapshot
built against the wrong transcript version yields an empty table, a region slice
with the wrong coordinates yields a table that covers nothing, and both look
exactly like "this variant simply has no data" once they reach the engine.

So the check is separate from the run, and it reports **coverage** rather than
mere presence. A frequency snapshot with ten thousand rows that intersects the
transcript in none of them is worse than an absent one, because the absent one
is honest.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .adapters.base import Adapter
from .genome.reference import GeneConfig, ReferenceUnavailable
from .genome.transcript import Transcript
from .variant import Variant

__all__ = ["SourceStatus", "DataReport", "check_reference", "check_sources"]


@dataclass(frozen=True, slots=True)
class SourceStatus:
    """One source, as it currently stands on disk."""

    namespace: str
    source: str
    present: bool
    records: int | None = None
    #: Fraction of the sampled variants for which the source returns anything.
    coverage: float | None = None
    detail: str | None = None

    @property
    def symbol(self) -> str:
        if not self.present:
            return "--"
        if self.coverage is None:
            return "ok"
        if self.coverage == 0.0:
            return "!!"
        return "ok" if self.coverage > 0.5 else " ~"

    def line(self) -> str:
        records = "" if self.records is None else f"{self.records:>10,} records"
        coverage = "" if self.coverage is None else f"  {self.coverage:6.1%} of sample"
        detail = f"  {self.detail}" if self.detail else ""
        return f"  {self.symbol}  {self.namespace:<12} {self.source:<28}{records}{coverage}{detail}"


@dataclass(slots=True)
class DataReport:
    gene: str
    reference: SourceStatus
    sources: list[SourceStatus] = field(default_factory=list)
    sampled: int = 0

    @property
    def usable(self) -> bool:
        """Whether a map can be computed at all. Sources may legitimately be absent."""
        return self.reference.present

    @property
    def empty_but_present(self) -> list[SourceStatus]:
        """Loaded and covering nothing -- the quiet failure worth shouting about."""
        return [s for s in self.sources if s.present and s.coverage == 0.0]

    def summary(self) -> str:
        lines = [f"{self.gene}", self.reference.line()]
        lines.extend(status.line() for status in self.sources)
        if self.sampled:
            lines.append(f"  (coverage measured over {self.sampled:,} sampled variants)")
        for status in self.empty_but_present:
            lines.append(
                f"  !! {status.namespace} is loaded but matches none of the sampled "
                "variants. Check the region slice and the coordinate convention "
                "(chrom-pos-ref-alt on the genomic plus strand) before running a map."
            )
        return "\n".join(lines)


def check_reference(config: GeneConfig, *, data_root: str | Path) -> SourceStatus:
    """Whether the MANE Select resources are present and self-consistent."""
    from .genome.reference import build_transcript

    if not config.transcript.has_coordinates:
        return SourceStatus(
            namespace="reference",
            source=config.transcript.id,
            present=False,
            detail="no exon coordinates; run 'vus-foresight reference import'",
        )
    try:
        transcript = build_transcript(config, data_root=data_root)
    except ReferenceUnavailable as exc:
        return SourceStatus(
            namespace="reference",
            source=config.transcript.id,
            present=False,
            detail=str(exc).split(".")[0],
        )
    return SourceStatus(
        namespace="reference",
        source=config.transcript.id,
        present=True,
        records=len(transcript.exons),
        detail=(
            f"CDS {transcript.cds_length} nt, {transcript.protein_length} residues, "
            f"strand {transcript.strand}"
        ),
    )


def check_sources(
    adapters: Iterable[Adapter],
    transcript: Transcript,
    sample: list[Variant],
) -> list[SourceStatus]:
    """Measure how much of a sample each adapter actually answers for.

    The sample is drawn from the enumeration rather than from observed variants,
    which is the only way to notice that a slice covers the wrong region: every
    possible variant in the transcript is in scope, so a source that answers for
    none of them is answering for nothing.
    """
    statuses: list[SourceStatus] = []
    for adapter in adapters:
        answered = sum(1 for v in sample if adapter.lookup(v, transcript))
        described = adapter.describe()
        records = described.get("records")
        min_stars = described.get("min_stars")
        statuses.append(
            SourceStatus(
                namespace=adapter.namespace,
                source=adapter.source_id,
                present=True,
                records=int(records) if records is not None else _len_or_none(adapter),
                coverage=answered / len(sample) if sample else None,
                detail=None if min_stars is None else f"min_stars={min_stars}",
            )
        )
    return sorted(statuses, key=lambda s: s.namespace)


def _len_or_none(adapter: Adapter) -> int | None:
    try:
        return len(adapter)  # type: ignore[arg-type]
    except TypeError:
        return None
