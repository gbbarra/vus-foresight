"""Build a dated snapshot from ClinVar's ``variant_summary`` release.

Run by a curator, once per snapshot date, from a file already on disk. Nothing
here is fetched at runtime.

The parser reads columns **by header name**, never by position: ClinVar has
added columns between releases, and a positional parser silently reads the
wrong field the first time that happens -- which would look like a data change
rather than a bug.
"""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import IO, Iterable, Iterator

from .clinvar import ClinVarRecord, normalise_classification, review_stars

__all__ = ["ParseStats", "parse_variant_summary", "write_snapshot", "build_snapshot"]

#: ``NM_007294.4(BRCA1):c.5074G>A (p.Asp1692Asn)`` -- the ``Name`` column.
_NAME = re.compile(
    r"^(?P<transcript>[A-Z_0-9.]+)"
    r"(?:\((?P<gene>[^)]+)\))?"
    r":(?P<hgvs_c>[cn]\.[^ ]+)"
    r"(?:\s+\((?P<hgvs_p>p\.[^)]+)\))?"
)

#: Residue number inside a three-letter ``p.`` description.
_RESIDUE = re.compile(r"^p\.\(?[A-Z][a-z]{2}(?P<residue>\d+)")


@dataclass(slots=True)
class ParseStats:
    """What the import saw, so a thin snapshot is noticed rather than shipped."""

    rows_read: int = 0
    wrong_assembly: int = 0
    other_transcript: int = 0
    unparsable_name: int = 0
    unmapped_classification: int = 0
    kept: int = 0

    def summary(self) -> str:
        return (
            f"{self.rows_read:,} rows read, {self.kept:,} kept "
            f"({self.other_transcript:,} other transcript, "
            f"{self.wrong_assembly:,} other assembly, "
            f"{self.unparsable_name:,} unparsable name, "
            f"{self.unmapped_classification:,} unmapped significance)"
        )


def _open(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _residue(hgvs_p: str | None) -> int | None:
    if not hgvs_p:
        return None
    match = _RESIDUE.match(hgvs_p)
    return int(match.group("residue")) if match else None


def parse_variant_summary(
    path: str | Path,
    *,
    transcript: str,
    assembly: str = "GRCh38",
    stats: ParseStats | None = None,
) -> Iterator[ClinVarRecord]:
    """Yield records for one transcript from a ``variant_summary`` file.

    ``transcript`` is matched with its version, because a classification made
    against ``NM_007294.3`` is not evidence about a coordinate in
    ``NM_007294.4`` unless somebody has checked that the two agree. Matching
    loosely would quietly mix transcripts, which is exactly what the MANE-only
    rule exists to prevent.
    """
    stats = stats if stats is not None else ParseStats()
    with _open(Path(path)) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: no header row")
        fields = {name.lstrip("#").strip(): name for name in reader.fieldnames}
        for required in ("Name", "ClinicalSignificance", "ReviewStatus"):
            if required not in fields:
                raise ValueError(
                    f"{path}: missing required column {required!r}. Columns seen: {sorted(fields)}"
                )
        assembly_column = fields.get("Assembly")
        evaluated_column = fields.get("LastEvaluated")

        for row in reader:
            stats.rows_read += 1
            if assembly_column and row.get(assembly_column, "").strip() not in (assembly, ""):
                stats.wrong_assembly += 1
                continue

            match = _NAME.match(row[fields["Name"]].strip())
            if match is None:
                stats.unparsable_name += 1
                continue
            if match.group("transcript") != transcript:
                stats.other_transcript += 1
                continue

            classification = normalise_classification(row[fields["ClinicalSignificance"]])
            if classification is None:
                stats.unmapped_classification += 1
                continue

            hgvs_p = match.group("hgvs_p")
            evaluated = (row.get(evaluated_column) or "").strip() if evaluated_column else ""
            stats.kept += 1
            yield ClinVarRecord(
                hgvs_c=match.group("hgvs_c"),
                hgvs_p=hgvs_p,
                classification=classification,
                stars=review_stars(row[fields["ReviewStatus"]]),
                codon=_residue(hgvs_p),
                last_evaluated=_parse_date(evaluated),
            )


def _parse_date(text: str) -> date | None:
    """ClinVar writes ``Mar 21, 2024``; a missing value is ``-``."""
    if not text or text == "-":
        return None
    for fmt in ("%b %d, %Y", "%Y-%m-%d"):
        try:
            from datetime import datetime

            return datetime.strptime(text, fmt).date()
        except ValueError:
            continue
    return None


def write_snapshot(records: Iterable[ClinVarRecord], path: str | Path) -> int:
    """Write the normalised snapshot TSV, sorted, so two builds match byte for byte."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda r: r.hgvs_c)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["hgvs_c", "hgvs_p", "codon", "classification", "stars", "last_evaluated"])
        for record in ordered:
            writer.writerow(
                [
                    record.hgvs_c,
                    record.hgvs_p or "",
                    "" if record.codon is None else record.codon,
                    record.classification,
                    record.stars,
                    record.last_evaluated.isoformat() if record.last_evaluated else "",
                ]
            )
    return len(ordered)


def build_snapshot(
    source: str | Path,
    destination: str | Path,
    *,
    transcript: str,
    assembly: str = "GRCh38",
) -> ParseStats:
    """Parse a ``variant_summary`` release and write a normalised snapshot."""
    stats = ParseStats()
    records = list(
        parse_variant_summary(source, transcript=transcript, assembly=assembly, stats=stats)
    )
    write_snapshot(records, destination)
    return stats
