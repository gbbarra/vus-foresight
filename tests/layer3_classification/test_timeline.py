"""Layer 3: recomputing the map over time and attributing what moved.

The claim being tested is the one that justifies separating semi-intrinsic
evidence from intrinsic evidence at all: **a variant can leave VUS without any
new evidence about itself appearing**, because a neighbour in the same codon was
classified. If the diff cannot tell that apart from a new assay result, the
distinction is decorative.
"""

from __future__ import annotations

from datetime import date, datetime

import pytest

from vus_foresight.acmg import ACMGClass
from vus_foresight.adapters import AdapterRegistry, ClinVarSnapshotAdapter
from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from vus_foresight.adapters.clinvar import ClinVarRecord
from vus_foresight.adapters.clinvar_import import write_snapshot
from vus_foresight.adapters.tabular import PredictorAdapter
from vus_foresight.engine.pipeline import MapRunner
from vus_foresight.engine.timeline import TransitionCause, compare_maps, compare_series
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.variant import Consequence


@pytest.fixture()
def missense_variants(minus_gene):
    return [
        v
        for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE and v.codon_index and v.codon_index > 3
    ]


def _map(
    minus_gene,
    minus_config,
    toy_spec,
    variants,
    *,
    clinvar_path=None,
    predictor_path=None,
    snapshot_date=None,
):
    registry = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(minus_config)])
    if clinvar_path is not None:
        registry.add(ClinVarSnapshotAdapter.from_path(clinvar_path, snapshot_date=snapshot_date))
    if predictor_path is not None:
        registry.add(PredictorAdapter(predictor_path, "dbnsfp-test"))
    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=registry,
        computed_at=datetime(1970, 1, 1),
        clinvar_snapshot=snapshot_date,
    )
    return [r.row for r in runner.run(variants)]


def test_a_neighbour_classification_moves_a_variant_with_no_evidence_about_itself(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """The section 4 observation, end to end."""
    target = missense_variants[0]
    variants = missense_variants[:40]

    empty = tmp_path / "clinvar_2018.tsv"
    write_snapshot([], empty)

    later = tmp_path / "clinvar_2024.tsv"
    write_snapshot(
        [
            # A different nucleotide change reaching the same protein change.
            ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)
        ],
        later,
    )

    before = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=empty,
        snapshot_date=date(2018, 1, 1),
    )
    after = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=later,
        snapshot_date=date(2024, 1, 1),
    )

    diff = compare_maps(before, after, label_before="2018-01-01", label_after="2024-01-01")
    assert diff.compared == len(variants)
    assert not diff.appeared and not diff.disappeared

    transition = next(t for t in diff.transitions if t.variant_id == target.variant_id)
    assert transition.criteria_gained == (("PS1", "strong"),)
    assert transition.criteria_lost == ()
    assert transition.cause is TransitionCause.NEIGHBOUR_EVIDENCE
    assert transition.points_after == transition.points_before + 4

    # One classification does not move one variant -- it moves the whole codon.
    # The variants reaching the same protein change gain PS1; the rest of the
    # codon gains PM5. Not one of them had any evidence generated about itself.
    by_codon = {v.variant_id: v.codon_index for v in variants if v.codon_index is not None}
    assert {by_codon[t.variant_id] for t in diff.transitions} == {target.codon_index}
    assert all(t.cause is TransitionCause.NEIGHBOUR_EVIDENCE for t in diff.transitions)
    gained = {t.variant_id: {code for code, _ in t.criteria_gained} for t in diff.transitions}
    assert all(codes <= {"PS1", "PM5"} for codes in gained.values())
    assert any(codes == {"PS1"} for codes in gained.values())
    assert any(codes == {"PM5"} for codes in gained.values())

    # And nothing outside that codon moved at all.
    untouched = {v.variant_id for v in variants} - set(gained)
    assert untouched
    assert all(by_codon[vid] != target.codon_index for vid in untouched)


def test_new_intrinsic_data_is_attributed_differently(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    variants = missense_variants[:40]
    target = variants[0]

    empty_predictor = tmp_path / "pred_empty.tsv"
    empty_predictor.write_text("grch38_pos\tbayesdel\n", encoding="utf-8")
    scored = tmp_path / "pred_scored.tsv"
    scored.write_text(f"grch38_pos\tbayesdel\n{target.grch38_pos}\t0.62\n", encoding="utf-8")

    before = _map(minus_gene, minus_config, toy_spec, variants, predictor_path=empty_predictor)
    after = _map(minus_gene, minus_config, toy_spec, variants, predictor_path=scored)

    diff = compare_maps(before, after)
    transition = next(t for t in diff.transitions if t.variant_id == target.variant_id)
    assert transition.criteria_gained == (("PP3", "moderate"),)
    assert transition.cause is TransitionCause.INTRINSIC_EVIDENCE
    assert not transition.resolved_without_self_evidence


def test_both_kinds_moving_in_one_interval_is_reported_as_mixed(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    variants = missense_variants[:40]
    target = variants[0]

    empty_cv = tmp_path / "cv_empty.tsv"
    write_snapshot([], empty_cv)
    later_cv = tmp_path / "cv_later.tsv"
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later_cv,
    )
    empty_pred = tmp_path / "pred_empty.tsv"
    empty_pred.write_text("grch38_pos\tbayesdel\n", encoding="utf-8")
    scored = tmp_path / "pred_scored.tsv"
    scored.write_text(f"grch38_pos\tbayesdel\n{target.grch38_pos}\t0.62\n", encoding="utf-8")

    before = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=empty_cv,
        predictor_path=empty_pred,
    )
    after = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=later_cv,
        predictor_path=scored,
    )

    transition = next(
        t for t in compare_maps(before, after).transitions if t.variant_id == target.variant_id
    )
    assert {code for code, _ in transition.criteria_gained} == {"PS1", "PP3"}
    assert transition.cause is TransitionCause.MIXED_EVIDENCE


def test_resolution_without_self_evidence_is_singled_out(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """A variant that crosses into LP purely on a neighbour's classification."""
    variants = missense_variants[:40]
    target = variants[0]

    # Start with enough intrinsic weight that one strong criterion tips it over.
    predictor = tmp_path / "pred.tsv"
    predictor.write_text(f"grch38_pos\tbayesdel\n{target.grch38_pos}\t0.62\n", encoding="utf-8")
    empty_cv = tmp_path / "cv_empty.tsv"
    write_snapshot([], empty_cv)
    later_cv = tmp_path / "cv_later.tsv"
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later_cv,
    )

    before = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=empty_cv,
        predictor_path=predictor,
    )
    after = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        clinvar_path=later_cv,
        predictor_path=predictor,
    )

    diff = compare_maps(before, after)
    transition = next(t for t in diff.transitions if t.variant_id == target.variant_id)
    assert transition.class_before is ACMGClass.UNCERTAIN
    assert transition.class_after is ACMGClass.LIKELY_PATHOGENIC
    assert transition.resolved
    assert transition.resolved_without_self_evidence
    assert diff.neighbour_driven_resolutions == [transition]
    assert "no new evidence about the variant itself" in diff.summary()


def test_rows_present_in_only_one_snapshot_are_reported_separately(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """Enumeration moving is a different event from evidence moving."""
    before = _map(minus_gene, minus_config, toy_spec, missense_variants[:20])
    after = _map(minus_gene, minus_config, toy_spec, missense_variants[5:25])
    diff = compare_maps(before, after)
    assert len(diff.appeared) == 5
    assert len(diff.disappeared) == 5
    assert diff.compared == 15
    assert not diff.transitions


def test_an_unchanged_interval_produces_no_transitions(
    minus_gene, minus_config, toy_spec, missense_variants
):
    rows = _map(minus_gene, minus_config, toy_spec, missense_variants[:30])
    diff = compare_maps(rows, rows)
    assert diff.transitions == []
    assert diff.class_changes == []
    assert diff.counts_by_cause() == {}
    assert compare_maps(rows, rows, include_unchanged=True).compared == 30


def test_a_series_is_diffed_consecutively(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """Consecutive, not all-against-the-first: *when* a variant moved is the point."""
    variants = missense_variants[:30]
    first = variants[0]
    # A different codon, so the two intervals cannot bleed into each other: one
    # classification moves an entire codon, not a single variant.
    second = next(v for v in variants if v.codon_index != first.codon_index)

    snapshots = []
    for index, records in enumerate(
        [
            [],
            [ClinVarRecord("c.999999A>G", first.hgvs_p, "pathogenic", 3, first.codon_index)],
            [
                ClinVarRecord("c.999999A>G", first.hgvs_p, "pathogenic", 3, first.codon_index),
                ClinVarRecord("c.999998A>G", second.hgvs_p, "pathogenic", 3, second.codon_index),
            ],
        ]
    ):
        path = tmp_path / f"cv_{index}.tsv"
        write_snapshot(records, path)
        snapshots.append(
            (
                f"T{index}",
                _map(
                    minus_gene,
                    minus_config,
                    toy_spec,
                    variants,
                    clinvar_path=path,
                    snapshot_date=date(2018 + index, 1, 1),
                ),
            )
        )

    diffs = compare_series(snapshots)
    assert [d.label_before for d in diffs] == ["T0", "T1"]
    codon_of = {v.variant_id: v.codon_index for v in variants}
    # Each interval moves exactly the codon whose neighbour was classified in it.
    assert {codon_of[t.variant_id] for t in diffs[0].transitions} == {first.codon_index}
    assert {codon_of[t.variant_id] for t in diffs[1].transitions} == {second.codon_index}
    assert all(
        t.cause is TransitionCause.NEIGHBOUR_EVIDENCE for diff in diffs for t in diff.transitions
    )


def test_a_series_needs_at_least_two_snapshots(minus_gene, minus_config, toy_spec):
    with pytest.raises(ValueError, match="at least two"):
        compare_series([("only", [])])


def test_transitions_survive_a_parquet_round_trip(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """The CLI diffs written maps, so the rehydration has to keep the criteria.

    A partial rehydration that dropped ``criteria_applied`` would make every
    transition look like a specification change -- silently, and only in the
    command anyone would actually run.
    """
    from vus_foresight.output import read_rows, write_parquet

    variants = missense_variants[:20]
    target = variants[0]
    empty = tmp_path / "cv_empty.tsv"
    write_snapshot([], empty)
    later = tmp_path / "cv_later.tsv"
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later,
    )

    before_path = tmp_path / "before.parquet"
    after_path = tmp_path / "after.parquet"
    assert write_parquet(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=empty), before_path
    ) == len(variants)
    assert write_parquet(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=later), after_path
    ) == len(variants)

    in_memory = compare_maps(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=empty),
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=later),
    )
    round_tripped = compare_maps(read_rows(before_path), read_rows(after_path))

    assert round_tripped.transitions == in_memory.transitions
    assert round_tripped.transitions
    assert all(t.cause is TransitionCause.NEIGHBOUR_EVIDENCE for t in round_tripped.transitions)
    assert any(t.criteria_gained == (("PS1", "strong"),) for t in round_tripped.transitions)


def test_streaming_a_map_from_disk_gives_the_identical_diff(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """The whole point of GapMapSource: same answer, without holding two maps.

    A real gene's map is ~9 GB once rehydrated and the diff needs two of them,
    so the streaming path is the only one that runs at full scale -- which makes
    it the one that has to be pinned to the reference implementation.
    """
    from vus_foresight.output import GapMapSource, read_rows, write_parquet

    variants = missense_variants[:30]
    target = variants[0]
    empty = tmp_path / "cv_empty.tsv"
    write_snapshot([], empty)
    later = tmp_path / "cv_later.tsv"
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later,
    )
    before_path, after_path = tmp_path / "b.parquet", tmp_path / "a.parquet"
    write_parquet(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=empty), before_path
    )
    write_parquet(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=later), after_path
    )

    materialised = compare_maps(read_rows(before_path), read_rows(after_path))
    # A batch size below the row count forces several Parquet batches, so the
    # test covers the boundary rather than a single-batch shortcut.
    streamed = compare_maps(
        GapMapSource(before_path, batch_size=7), GapMapSource(after_path, batch_size=7)
    )

    assert streamed.transitions == materialised.transitions
    assert streamed.compared == materialised.compared == len(variants)
    assert streamed.appeared == materialised.appeared
    assert streamed.disappeared == materialised.disappeared
    assert streamed.transitions


def test_transitions_come_out_in_enumeration_order(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """Streaming fixed the order to the later map's, which is the enumerator's."""
    from vus_foresight.output import GapMapSource, write_parquet

    variants = missense_variants[:30]
    target = variants[0]
    empty, later = tmp_path / "e.tsv", tmp_path / "l.tsv"
    write_snapshot([], empty)
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later,
    )
    before_path, after_path = tmp_path / "b.parquet", tmp_path / "a.parquet"
    write_parquet(
        _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=empty), before_path
    )
    after_rows = _map(minus_gene, minus_config, toy_spec, variants, clinvar_path=later)
    write_parquet(after_rows, after_path)

    diff = compare_maps(GapMapSource(before_path), GapMapSource(after_path))
    moved = [t.variant_id for t in diff.transitions]
    order = [row.variant_id for row in after_rows]
    assert moved == [vid for vid in order if vid in set(moved)]


def test_a_one_shot_iterator_is_refused_rather_than_silently_emptied(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """A generator would make the second interval report everything as dropped."""
    rows = _map(minus_gene, minus_config, toy_spec, missense_variants[:10])
    with pytest.raises(TypeError, match="one-shot iterator"):
        compare_series([("T0", iter(rows)), ("T1", iter(rows))])


def test_parquet_round_trip_preserves_every_field(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    from vus_foresight.output import read_rows, write_parquet

    original = _map(minus_gene, minus_config, toy_spec, missense_variants[:25])
    path = tmp_path / "roundtrip.parquet"
    assert write_parquet(original, path) == len(original)
    assert read_rows(path) == original
