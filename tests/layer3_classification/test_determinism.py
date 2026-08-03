"""Layer 3: determinism and the golden-snapshot regression.

"Same input + same spec version + same database snapshot = same output, always."
Anything that breaks it -- dictionary iteration order, parallelism, a timestamp
leaking into the data -- fails here rather than in a quarterly comparison
nobody can explain.

The golden snapshot is the audit trail for specification changes: a diff in the
committed output has to be justified in the pull request that causes it.
"""

from __future__ import annotations

import hashlib
import os
from datetime import datetime
from pathlib import Path

import pytest

from vus_foresight.engine.pipeline import MapRunner, default_registry
from vus_foresight.enumeration import EnumerationClass, enumerate_all
from vus_foresight.output import rows_to_frame, write_parquet

GOLDEN = Path(__file__).resolve().parents[1] / "fixtures" / "golden"
GOLDEN_FILE = GOLDEN / "toy_minus_gap_map.tsv"
REGENERATE = os.environ.get("VUS_FORESIGHT_REGENERATE_GOLDEN") == "1"

FIXED_STAMP = datetime(1970, 1, 1)


def _rows(minus_gene, minus_config, toy_spec, limit: int | None = None):
    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=default_registry(minus_config),
        computed_at=FIXED_STAMP,
    )
    variants = list(
        enumerate_all(
            minus_gene.transcript,
            classes=(
                EnumerationClass.CODING_SNV,
                EnumerationClass.FRAMESHIFT_CLASS,
                EnumerationClass.INFRAME_DELETION,
            ),
            flanks=minus_gene.flanks,
        )
    )
    if limit is not None:
        variants = variants[:limit]
    return [result.row for result in runner.run(variants)]


def test_two_runs_produce_byte_identical_parquet(
    tmp_path, minus_gene, minus_config, toy_spec
):
    first = write_parquet(_rows(minus_gene, minus_config, toy_spec, 300), tmp_path / "a.parquet")
    second = write_parquet(_rows(minus_gene, minus_config, toy_spec, 300), tmp_path / "b.parquet")
    digest_a = hashlib.sha256(first.read_bytes()).hexdigest()
    digest_b = hashlib.sha256(second.read_bytes()).hexdigest()
    assert digest_a == digest_b


def test_row_order_follows_the_enumerator_not_a_set(minus_gene, minus_config, toy_spec):
    first = [row.hgvs_c for row in _rows(minus_gene, minus_config, toy_spec, 200)]
    second = [row.hgvs_c for row in _rows(minus_gene, minus_config, toy_spec, 200)]
    assert first == second
    assert first[:3] == ["c.1A>C", "c.1A>G", "c.1A>T"]


def test_criteria_within_a_row_are_deterministically_ordered(
    minus_gene, minus_config, toy_spec
):
    for row in _rows(minus_gene, minus_config, toy_spec, 100):
        applied = [c.code for c in row.criteria_applied]
        skipped = [s.code for s in row.criteria_evaluated_not_applied]
        assert applied == sorted(applied)
        assert skipped == sorted(skipped)


def test_source_versions_are_stored_sorted(minus_gene, minus_config, toy_spec):
    for row in _rows(minus_gene, minus_config, toy_spec, 20):
        keys = list(row.source_versions)
        assert keys == sorted(keys)


def test_equivalence_class_ids_are_stable_across_runs(
    minus_gene, minus_config, toy_spec
):
    first = {r.hgvs_c: r.equivalence_class_id for r in _rows(minus_gene, minus_config, toy_spec, 200)}
    second = {r.hgvs_c: r.equivalence_class_id for r in _rows(minus_gene, minus_config, toy_spec, 200)}
    assert first == second


def test_variants_with_the_same_intrinsic_profile_share_a_class(
    minus_gene, minus_config, toy_spec
):
    """And missense never collapses, per spec section 5."""
    from vus_foresight.variant import Consequence

    rows = _rows(minus_gene, minus_config, toy_spec)
    by_class: dict[str, list] = {}
    for row in rows:
        by_class.setdefault(row.equivalence_class_id, []).append(row)

    missense_classes = [
        members
        for members in by_class.values()
        if any(r.consequence is Consequence.MISSENSE for r in members)
    ]
    assert missense_classes, "the fixture should produce missense variants"
    assert all(len(members) == 1 for members in missense_classes)

    synonymous_classes = [
        members
        for members in by_class.values()
        if all(r.consequence is Consequence.SYNONYMOUS for r in members)
    ]
    assert any(len(members) > 1 for members in synonymous_classes), (
        "synonymous variants with an identical evidence profile should collapse"
    )


def _golden_lines(rows) -> list[str]:
    header = "\t".join(
        [
            "hgvs_c",
            "hgvs_p",
            "consequence",
            "equivalence_class_id",
            "points_current",
            "class_current",
            "points_ceiling_intrinsic",
            "gap_to_LP",
            "gap_to_LB",
            "blocking_reason",
            "applied",
            "minimum_sufficient_sets",
        ]
    )
    lines = [header]
    for row in rows:
        applied = ";".join(f"{c.code}:{c.strength.value}:{c.points}" for c in row.criteria_applied)
        sets = "|".join(
            f"{s.target.value}[{'+'.join(s.codes)}]={s.total_points}:{s.feasibility.value}"
            for s in row.minimum_sufficient_sets
        )
        lines.append(
            "\t".join(
                [
                    row.hgvs_c,
                    row.hgvs_p or "",
                    row.consequence.value,
                    row.equivalence_class_id,
                    str(row.points_current),
                    row.class_current.value,
                    str(row.points_ceiling_intrinsic),
                    "" if row.gap_to_LP is None else str(row.gap_to_LP),
                    "" if row.gap_to_LB is None else str(row.gap_to_LB),
                    row.blocking_reason.value,
                    applied,
                    sets,
                ]
            )
        )
    return lines


def test_golden_snapshot(minus_gene, minus_config, toy_spec):
    """A committed output for a fixed subset. Any diff needs a reason in the PR.

    Regenerate deliberately::

        VUS_FORESIGHT_REGENERATE_GOLDEN=1 pytest tests/layer3_classification -k golden
    """
    rows = _rows(minus_gene, minus_config, toy_spec, 400)
    produced = _golden_lines(rows)

    if REGENERATE or not GOLDEN_FILE.exists():
        GOLDEN.mkdir(parents=True, exist_ok=True)
        GOLDEN_FILE.write_text("\n".join(produced) + "\n", encoding="utf-8")
        if not REGENERATE:
            pytest.skip(f"created {GOLDEN_FILE.name}; commit it and re-run")
        return

    expected = GOLDEN_FILE.read_text(encoding="utf-8").splitlines()
    if produced != expected:
        differences = [
            f"line {i}: expected {e!r}, produced {p!r}"
            for i, (e, p) in enumerate(zip(expected, produced), start=1)
            if e != p
        ]
        pytest.fail(
            f"{len(differences)} golden line(s) changed (and "
            f"{abs(len(expected) - len(produced))} added or removed). "
            "If the change is intended, regenerate the snapshot and justify the diff "
            f"in the pull request. First five:\n" + "\n".join(differences[:5])
        )


def test_frame_schema_is_fixed_not_inferred(minus_gene, minus_config, toy_spec):
    from vus_foresight.gapmap import GAP_MAP_COLUMNS

    frame = rows_to_frame(_rows(minus_gene, minus_config, toy_spec, 50))
    assert tuple(frame.columns) == GAP_MAP_COLUMNS
