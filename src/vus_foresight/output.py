"""Parquet serialisation and the aggregate reports.

Two requirements shape this module.

*Determinism* (spec section 11): running the pipeline twice must produce
byte-identical Parquet. So the schema is written out explicitly rather than
inferred, dictionaries become sorted key/value lists rather than maps, and no
value is derived from iteration order or from the clock.

*Queryability*: the output has to answer "how many BRCA2 VUS are blocked for
want of functional data" as one ``GROUP BY``, from DuckDB, without a Python
process in the loop.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import polars as pl

from .acmg import (
    ACMGClass,
    BlockingReason,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    SkipReason,
    Strength,
)
from .gapmap import (
    GAP_MAP_COLUMNS,
    AppliedCriterion,
    EvidenceRequirement,
    EvidenceSet,
    GapMapRow,
    SkippedCriterion,
)
from .variant import Consequence, VariantKind

__all__ = [
    "GAP_MAP_SCHEMA",
    "PARQUET_BATCH_SIZE",
    "rows_to_frame",
    "rows_from_frame",
    "read_rows",
    "write_parquet",
    "read_parquet",
    "blocking_summary",
    "available_uningested_report",
    "equivalence_summary",
]

_APPLIED = pl.Struct(
    {
        "code": pl.Utf8,
        "direction": pl.Utf8,
        "strength": pl.Utf8,
        "points": pl.Int32,
        "evidence_class": pl.Utf8,
        "evidence": pl.Utf8,
        "source": pl.Utf8,
    }
)

_SKIPPED = pl.Struct(
    {
        "code": pl.Utf8,
        "direction": pl.Utf8,
        "evidence_class": pl.Utf8,
        "outcome": pl.Utf8,
        "reason": pl.Utf8,
        "missing_fields": pl.List(pl.Utf8),
        "detail": pl.Utf8,
    }
)

_REQUIREMENT = pl.Struct(
    {
        "code": pl.Utf8,
        "direction": pl.Utf8,
        "strength": pl.Utf8,
        "points": pl.Int32,
        "evidence_class": pl.Utf8,
        "feasibility": pl.Utf8,
        "required_observations": pl.Int32,
        "description": pl.Utf8,
    }
)

_EVIDENCE_SET = pl.Struct(
    {
        "target": pl.Utf8,
        "requirements": pl.List(_REQUIREMENT),
        "total_points": pl.Int32,
        "feasibility": pl.Utf8,
        "acquisition_cost": pl.Int32,
    }
)

#: Explicit schema. Inference would let an all-null column in one shard become a
#: different type from the same column in another, which silently breaks a
#: multi-gene read.
GAP_MAP_SCHEMA: dict[str, Any] = {
    "gene": pl.Utf8,
    "transcript": pl.Utf8,
    "hgvs_c": pl.Utf8,
    "hgvs_p": pl.Utf8,
    "grch38_pos": pl.Utf8,
    "consequence": pl.Utf8,
    "variant_kind": pl.Utf8,
    "equivalence_class_id": pl.Utf8,
    "mutational_distance": pl.Int32,
    "criteria_applied": pl.List(_APPLIED),
    "criteria_evaluated_not_applied": pl.List(_SKIPPED),
    "points_current": pl.Int32,
    "class_current": pl.Utf8,
    "points_ceiling_intrinsic": pl.Int32,
    "class_ceiling_intrinsic": pl.Utf8,
    "gap_to_LP": pl.Int32,
    "gap_to_LB": pl.Int32,
    "minimum_sufficient_sets": pl.List(_EVIDENCE_SET),
    "blocking_reason": pl.Utf8,
    "spec_version": pl.Utf8,
    "clinvar_snapshot": pl.Date,
    "gnomad_version": pl.Utf8,
    "source_versions": pl.List(pl.Struct({"key": pl.Utf8, "value": pl.Utf8})),
    "computed_at": pl.Datetime("us"),
}


def _row_to_record(row: GapMapRow) -> dict[str, Any]:
    return {
        "gene": row.gene,
        "transcript": row.transcript,
        "hgvs_c": row.hgvs_c,
        "hgvs_p": row.hgvs_p,
        "grch38_pos": row.grch38_pos,
        "consequence": row.consequence.value,
        "variant_kind": row.variant_kind.value,
        "equivalence_class_id": row.equivalence_class_id,
        "mutational_distance": row.mutational_distance,
        "criteria_applied": [
            {
                "code": c.code,
                "direction": c.direction.value,
                "strength": c.strength.value,
                "points": c.points,
                "evidence_class": c.evidence_class.value,
                "evidence": c.evidence,
                "source": c.source,
            }
            for c in row.criteria_applied
        ],
        "criteria_evaluated_not_applied": [
            {
                "code": s.code,
                "direction": s.direction.value,
                "evidence_class": s.evidence_class.value,
                "outcome": s.outcome.value,
                "reason": s.reason.value,
                "missing_fields": list(s.missing_fields),
                "detail": s.detail,
            }
            for s in row.criteria_evaluated_not_applied
        ],
        "points_current": row.points_current,
        "class_current": row.class_current.value,
        "points_ceiling_intrinsic": row.points_ceiling_intrinsic,
        "class_ceiling_intrinsic": row.class_ceiling_intrinsic.value,
        "gap_to_LP": row.gap_to_LP,
        "gap_to_LB": row.gap_to_LB,
        "minimum_sufficient_sets": [
            {
                "target": s.target.value,
                "requirements": [
                    {
                        "code": r.code,
                        "direction": r.direction.value,
                        "strength": r.strength.value,
                        "points": r.points,
                        "evidence_class": r.evidence_class.value,
                        "feasibility": r.feasibility.value,
                        "required_observations": r.required_observations,
                        "description": r.description,
                    }
                    for r in s.requirements
                ],
                "total_points": s.total_points,
                "feasibility": s.feasibility.value,
                "acquisition_cost": s.acquisition_cost,
            }
            for s in row.minimum_sufficient_sets
        ],
        "blocking_reason": row.blocking_reason.value,
        "spec_version": row.spec_version,
        "clinvar_snapshot": row.clinvar_snapshot,
        "gnomad_version": row.gnomad_version,
        "source_versions": [
            {"key": key, "value": value} for key, value in sorted(row.source_versions.items())
        ],
        "computed_at": row.computed_at,
    }


def rows_to_frame(rows: Iterable[GapMapRow]) -> pl.DataFrame:
    """Materialise rows into a frame with the fixed schema and column order."""
    records = [_row_to_record(row) for row in rows]
    frame = pl.DataFrame(records, schema=GAP_MAP_SCHEMA)
    return frame.select(list(GAP_MAP_COLUMNS))


#: Rows converted per Parquet row group.
#:
#: Every row carries its complete evaluation trace -- two dozen nested criterion
#: structs plus the sufficient sets -- which is the design (spec section 1: the
#: trace *is* the map) and also makes rows expensive to materialise. Measured on
#: real BRCA1: 16,776 rows cost 3.1 GB to convert in one shot, so the full tier
#: 1+2 enumeration would need roughly 25 GB and simply dies. Batching caps the
#: conversion at a few hundred megabytes regardless of gene size.
#:
#: Fixed rather than tuned, because the batch boundary decides the row-group
#: boundary, and the determinism guarantee compares bytes.
PARQUET_BATCH_SIZE = 2000


def _batched(rows: Iterable[GapMapRow], size: int) -> Iterator[list[GapMapRow]]:
    batch: list[GapMapRow] = []
    for row in rows:
        batch.append(row)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def write_parquet(
    rows: Iterable[GapMapRow], path: str | Path, *, batch_size: int = PARQUET_BATCH_SIZE
) -> int:
    """Stream a gap map partition to Parquet, returning the number of rows.

    Takes an iterable and never holds more than one batch of converted rows, so
    a whole-gene map does not have to fit in memory at once.

    ``write_statistics=False`` is not an optimisation: min/max statistics embed
    per-row-group values whose encoding has varied between writer versions, and
    the determinism test compares bytes.
    """
    import pyarrow.parquet as pq

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    writer: "pq.ParquetWriter | None" = None
    written = 0
    try:
        for batch in _batched(rows, batch_size):
            table = rows_to_frame(batch).to_arrow()
            if writer is None:
                writer = pq.ParquetWriter(
                    path, table.schema, compression="zstd", write_statistics=False
                )
            writer.write_table(table)
            written += len(batch)
        if writer is None:
            # An empty map is still a valid answer, and it must carry the schema
            # so a downstream read of several genes does not change type.
            table = rows_to_frame([]).to_arrow()
            writer = pq.ParquetWriter(
                path, table.schema, compression="zstd", write_statistics=False
            )
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    return written


def read_parquet(path: str | Path) -> pl.DataFrame:
    return pl.read_parquet(path)


def rows_from_frame(frame: pl.DataFrame) -> list[GapMapRow]:
    """Rehydrate full :class:`GapMapRow` objects, criteria and all.

    The timeline diff attributes a transition by looking at *which* criteria
    moved and what evidence class they belong to, so a partial rehydration that
    drops ``criteria_applied`` would silently make every transition look like a
    specification change.
    """
    rows: list[GapMapRow] = []
    for record in frame.iter_rows(named=True):
        rows.append(
            GapMapRow(
                gene=record["gene"],
                transcript=record["transcript"],
                hgvs_c=record["hgvs_c"],
                hgvs_p=record["hgvs_p"],
                grch38_pos=record["grch38_pos"],
                consequence=Consequence(record["consequence"]),
                variant_kind=VariantKind(record["variant_kind"]),
                equivalence_class_id=record["equivalence_class_id"],
                mutational_distance=record["mutational_distance"],
                criteria_applied=tuple(
                    AppliedCriterion(
                        code=c["code"],
                        direction=Direction(c["direction"]),
                        strength=Strength(c["strength"]),
                        points=c["points"],
                        evidence_class=EvidenceClass(c["evidence_class"]),
                        evidence=c["evidence"],
                        source=c["source"],
                    )
                    for c in record["criteria_applied"]
                ),
                criteria_evaluated_not_applied=tuple(
                    SkippedCriterion(
                        code=s["code"],
                        direction=Direction(s["direction"]),
                        evidence_class=EvidenceClass(s["evidence_class"]),
                        outcome=CriterionOutcome(s["outcome"]),
                        reason=SkipReason(s["reason"]),
                        missing_fields=tuple(s["missing_fields"] or ()),
                        detail=s["detail"],
                    )
                    for s in record["criteria_evaluated_not_applied"]
                ),
                points_current=record["points_current"],
                class_current=ACMGClass(record["class_current"]),
                points_ceiling_intrinsic=record["points_ceiling_intrinsic"],
                class_ceiling_intrinsic=ACMGClass(record["class_ceiling_intrinsic"]),
                gap_to_LP=record["gap_to_LP"],
                gap_to_LB=record["gap_to_LB"],
                minimum_sufficient_sets=tuple(
                    EvidenceSet(
                        target=ACMGClass(s["target"]),
                        requirements=tuple(
                            EvidenceRequirement(
                                code=r["code"],
                                direction=Direction(r["direction"]),
                                strength=Strength(r["strength"]),
                                points=r["points"],
                                evidence_class=EvidenceClass(r["evidence_class"]),
                                feasibility=FeasibilityTag(r["feasibility"]),
                                required_observations=r["required_observations"],
                                description=r["description"],
                            )
                            for r in s["requirements"]
                        ),
                        total_points=s["total_points"],
                        feasibility=FeasibilityTag(s["feasibility"]),
                        acquisition_cost=s["acquisition_cost"],
                    )
                    for s in record["minimum_sufficient_sets"]
                ),
                blocking_reason=BlockingReason(record["blocking_reason"]),
                spec_version=record["spec_version"],
                clinvar_snapshot=record["clinvar_snapshot"],
                gnomad_version=record["gnomad_version"],
                source_versions={
                    entry["key"]: entry["value"] for entry in record["source_versions"]
                },
                computed_at=record["computed_at"],
            )
        )
    return rows


def read_rows(path: str | Path) -> list[GapMapRow]:
    """Read a gap map partition straight back into model objects."""
    return rows_from_frame(read_parquet(path))


def blocking_summary(frame: pl.DataFrame) -> pl.DataFrame:
    """``GROUP BY blocking_reason`` -- the aggregation the enum exists for."""
    return (
        frame.group_by("gene", "blocking_reason")
        .agg(pl.len().alias("variants"))
        .sort("gene", "blocking_reason")
    )


def equivalence_summary(frame: pl.DataFrame) -> pl.DataFrame:
    """Class-level view: the map read by region rather than by variant."""
    return (
        frame.group_by("gene", "equivalence_class_id")
        .agg(
            pl.len().alias("variants"),
            pl.first("consequence").alias("consequence"),
            pl.first("class_current").alias("class_current"),
            pl.first("blocking_reason").alias("blocking_reason"),
            pl.first("gap_to_LP").alias("gap_to_LP"),
        )
        .sort("variants", "equivalence_class_id", descending=[True, False])
    )


def available_uningested_report(frame: pl.DataFrame) -> pl.DataFrame:
    """Variants one publicly available, un-ingested dataset away from resolution.

    The single most actionable output the system produces (spec section 7): this
    is engineering work, not bench work, so it belongs at the top of any report.
    """
    tag = FeasibilityTag.AVAILABLE_UNINGESTED.value
    exploded = (
        frame.filter(pl.col("blocking_reason") != BlockingReason.RESOLVED_NOT_BLOCKED.value)
        .explode("minimum_sufficient_sets")
        .drop_nulls("minimum_sufficient_sets")
    )
    if exploded.height == 0:
        return pl.DataFrame(
            schema={"gene": pl.Utf8, "codes": pl.Utf8, "target": pl.Utf8, "variants": pl.UInt32}
        )
    return (
        exploded.filter(pl.col("minimum_sufficient_sets").struct.field("feasibility") == tag)
        .with_columns(
            pl.col("minimum_sufficient_sets")
            .struct.field("requirements")
            .list.eval(pl.element().struct.field("code"))
            .list.sort()
            .list.join("+")
            .alias("codes"),
            pl.col("minimum_sufficient_sets").struct.field("target").alias("target"),
        )
        .group_by("gene", "codes", "target")
        .agg(pl.len().alias("variants"))
        .sort("variants", "codes", descending=[True, False])
    )


def summarise_counts(frame: pl.DataFrame, columns: Sequence[str]) -> pl.DataFrame:
    """Generic count-by, used by the CLI's report command."""
    return frame.group_by(list(columns)).agg(pl.len().alias("variants")).sort(list(columns))
