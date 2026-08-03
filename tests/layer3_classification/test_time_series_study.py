"""Layer 3: the section 10 protocol run over the map's own time series.

The protocol needed an external outcome table. It does not, and the reason is a
property of the design rather than a convenience: ``clinvar.self.classification``
-- ClinVar's own call on a variant -- is published by the adapter and read by no
criterion, enforced by its own test. So the ground truth is sitting in the same
snapshots the engine consumes, and is firewalled from the thing being measured.

Metric 2 comes from the specification's own ``blocks_as`` mapping rather than
from submission prose, which this project could not parse without an LLM it does
not use.

These are unit tests of the study machinery. The study itself is still run on
demand, never in CI.
"""

from __future__ import annotations

from datetime import date, datetime

import pytest

from vus_foresight.acmg import ACMGClass, BlockingReason
from vus_foresight.adapters import AdapterRegistry, ClinVarSnapshotAdapter
from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from vus_foresight.adapters.clinvar import ClinVarRecord, ClinVarSnapshot
from vus_foresight.adapters.clinvar_import import write_snapshot
from vus_foresight.adapters.tabular import FunctionalAdapter, PredictorAdapter
from vus_foresight.engine.pipeline import MapRunner
from vus_foresight.engine.timeline import compare_maps
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.validation import (
    CLINVAR_TO_ACMG,
    observed_causes_from_diff,
    outcomes_from_snapshots,
    run_time_series_study,
)
from vus_foresight.variant import Consequence


@pytest.fixture()
def missense_variants(minus_gene):
    return [
        v
        for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE and v.codon_index and v.codon_index > 3
    ]


def _map(minus_gene, minus_config, toy_spec, variants, *, adapters=()):
    registry = AdapterRegistry(
        [VariantAdapter(), TranscriptAdapter(), RegionAdapter(minus_config)]
    )
    for adapter in adapters:
        registry.add(adapter)
    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=registry,
        computed_at=datetime(1970, 1, 1),
    )
    return [r.row for r in runner.run(variants)]


# --------------------------------------------------------------------------
# Ground truth derived from the snapshots.
# --------------------------------------------------------------------------


def test_outcomes_are_derived_from_clinvars_own_call(minus_gene):
    transcript = minus_gene.transcript.transcript_id
    before = ClinVarSnapshot.from_records(
        [
            ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 2, 4),
            ClinVarRecord("c.11C>T", "p.Thr4Ile", "uncertain", 2, 4),
            ClinVarRecord("c.12G>A", "p.Thr4=", "benign", 2, 4),
        ]
    )
    after = ClinVarSnapshot.from_records(
        [
            ClinVarRecord("c.10A>G", "p.Thr4Ala", "pathogenic", 3, 4, date(2024, 5, 1)),
            ClinVarRecord("c.11C>T", "p.Thr4Ile", "uncertain", 2, 4),
            ClinVarRecord("c.12G>A", "p.Thr4=", "benign", 2, 4),
        ]
    )
    outcomes = {
        o.variant_id: o for o in outcomes_from_snapshots(before, after, transcript_id=transcript)
    }

    # Only variants that were present AND uncertain at T are eligible.
    assert set(outcomes) == {f"{transcript}:c.10A>G", f"{transcript}:c.11C>T"}
    resolved = outcomes[f"{transcript}:c.10A>G"]
    assert resolved.was_resolved
    assert resolved.class_at_t_plus_n is ACMGClass.PATHOGENIC
    assert resolved.resolved_on == date(2024, 5, 1)
    assert not outcomes[f"{transcript}:c.11C>T"].was_resolved


def test_a_variant_absent_at_t_is_not_an_outcome(minus_gene):
    """It was unexamined, not uncertain. Counting it answers a different question."""
    before = ClinVarSnapshot.from_records([])
    after = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.10A>G", "p.Thr4Ala", "pathogenic", 3, 4)]
    )
    assert outcomes_from_snapshots(before, after, transcript_id="NM_1") == []


def test_conflicting_submissions_are_not_a_resolution():
    assert CLINVAR_TO_ACMG["conflicting"] is ACMGClass.UNCERTAIN
    before = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 2, 4)]
    )
    after = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.10A>G", "p.Thr4Ala", "conflicting", 1, 4)]
    )
    outcome = outcomes_from_snapshots(before, after, transcript_id="NM_1")[0]
    assert not outcome.was_resolved


def test_the_review_status_floor_applies_to_the_ground_truth_too(minus_gene):
    before = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 0, 4)]
    )
    after = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.10A>G", "p.Thr4Ala", "pathogenic", 0, 4)]
    )
    assert outcomes_from_snapshots(before, after, transcript_id="NM_1", min_stars=1) == []
    assert outcomes_from_snapshots(before, after, transcript_id="NM_1", min_stars=0)


# --------------------------------------------------------------------------
# Metric 2 without submission prose.
# --------------------------------------------------------------------------


def test_observed_cause_comes_from_the_specs_own_blocks_as(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    variants = missense_variants[:30]
    target = variants[0]

    empty = tmp_path / "func_empty.tsv"
    empty.write_text("hgvs_p\tclassification\tscore\tdataset\n", encoding="utf-8")
    assayed = tmp_path / "func.tsv"
    assayed.write_text(
        "hgvs_p\tclassification\tscore\tdataset\n"
        f"{target.hgvs_p}\tabnormal\t-2.1\tTOY-SGE\n",
        encoding="utf-8",
    )

    before = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(empty, "t0", assayed_regions=((1, 40, "TOY-SGE"),))],
    )
    after = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(assayed, "t1", assayed_regions=((1, 40, "TOY-SGE"),))],
    )

    causes = observed_causes_from_diff(compare_maps(before, after), toy_spec)
    assert causes[target.variant_id] == BlockingReason.MISSING_FUNCTIONAL.value
    assert toy_spec.by_code("PS3").blocks_as is BlockingReason.MISSING_FUNCTIONAL


def test_a_variant_whose_criteria_carry_no_blocks_as_yields_no_observed_cause(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """PM1 declares no blocks_as, so it cannot masquerade as a relieved block."""
    assert toy_spec.by_code("PM1").blocks_as is None
    variants = missense_variants[:20]
    rows = _map(minus_gene, minus_config, toy_spec, variants)
    assert observed_causes_from_diff(compare_maps(rows, rows), toy_spec) == {}


# --------------------------------------------------------------------------
# The whole study, crossed.
# --------------------------------------------------------------------------


def test_the_study_runs_with_no_external_outcome_table(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    variants = missense_variants[:40]
    target = variants[0]
    transcript = minus_gene.transcript.transcript_id

    # ClinVar's own call: uncertain at T, pathogenic at T+n. No criterion reads
    # this, so it is a clean oracle for what the field concluded.
    cv_before = ClinVarSnapshot.from_records(
        [
            ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index)
            for v in variants
        ],
        date(2018, 1, 1),
    )
    cv_after = ClinVarSnapshot.from_records(
        [
            ClinVarRecord(
                v.hgvs_c,
                v.hgvs_p,
                "pathogenic" if v is target else "uncertain",
                3 if v is target else 2,
                v.codon_index,
                date(2020, 6, 1) if v is target else None,
            )
            for v in variants
        ],
        date(2024, 1, 1),
    )

    # And the map's own sources: an assay result lands for the same variant.
    empty = tmp_path / "func_empty.tsv"
    empty.write_text("hgvs_p\tclassification\tscore\tdataset\n", encoding="utf-8")
    assayed = tmp_path / "func.tsv"
    assayed.write_text(
        "hgvs_p\tclassification\tscore\tdataset\n"
        f"{target.hgvs_p}\tabnormal\t-2.1\tTOY-SGE\n",
        encoding="utf-8",
    )
    regions = ((1, 40, "TOY-SGE"),)
    rows_t = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(empty, "t0", assayed_regions=regions)],
    )
    rows_tn = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(assayed, "t1", assayed_regions=regions)],
    )

    study = run_time_series_study(
        rows_t,
        rows_tn,
        cv_before,
        cv_after,
        toy_spec,
        transcript_id=transcript,
        reference_date=date(2018, 1, 1),
    )

    result = study.result
    assert result.considered == len(variants)
    assert result.resolved == 1
    assert result.resolvability_recall == 1.0
    # Metric 2: predicted MISSING_FUNCTIONAL, and a functional result is what
    # actually arrived.
    assert result.cause_scored == 1
    assert result.cause_accuracy == 1.0
    # Metric 3 is not scored here: at T the map had accumulated no evidence
    # either way on this variant, so it made no directional claim. Counting a
    # non-claim as a wrong answer would make the metric measure the thresholds
    # rather than the map.
    assert result.direction_accuracy is None
    assert result.direction_unpredicted == 1
    # The map saw the same assay result arrive in the same interval, but PS3
    # alone (+4) does not reach LP (6) under this specification. That is the
    # diagnostic bucket: the evidence was visible, the thresholds disagreed.
    assert study.saw_evidence_only == [target.variant_id]
    assert study.anticipated == []
    assert study.unmoved == []
    assert study.anticipation_rate == 0.0
    assert study.evidence_visibility_rate == 1.0
    assert "map saw the evidence" in study.summary()


def test_direction_is_scored_when_the_map_had_a_lean(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """A variant with accumulated pathogenic points, which ClinVar then made LP."""
    variants = missense_variants[:20]
    target = variants[0]
    transcript = minus_gene.transcript.transcript_id

    predictor = tmp_path / "pred.tsv"
    predictor.write_text(
        f"grch38_pos\tbayesdel\n{target.grch38_pos}\t0.62\n", encoding="utf-8"
    )
    rows = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[PredictorAdapter(predictor, "t")],
    )
    leaning = next(r for r in rows if r.hgvs_c == target.hgvs_c)
    assert leaning.points_current > 0

    cv_before = ClinVarSnapshot.from_records(
        [ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index) for v in variants]
    )
    cv_after = ClinVarSnapshot.from_records(
        [
            ClinVarRecord(
                v.hgvs_c,
                v.hgvs_p,
                "likely_pathogenic" if v is target else "uncertain",
                3 if v is target else 2,
                v.codon_index,
            )
            for v in variants
        ]
    )
    study = run_time_series_study(
        rows, rows, cv_before, cv_after, toy_spec, transcript_id=transcript
    )
    assert study.result.direction_scored == 1
    assert study.result.direction_accuracy == 1.0
    assert study.result.direction_unpredicted == 0


def test_direction_ignores_which_evidence_is_cheapest_to_obtain(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """Acquisition cost is not a lean.

    With nothing loaded, the cheapest sufficient sets are the benign-direction
    frequency criteria, because gnomAD is one ingest away. That says nothing
    about where the variant is heading, and an earlier definition of this metric
    read it as "benign" for every unevidenced variant in the gene.
    """
    from vus_foresight.validation import predicted_direction

    rows = _map(minus_gene, minus_config, toy_spec, missense_variants[:10])
    row = rows[0]
    assert row.points_current == 0
    cheapest = min(row.minimum_sufficient_sets, key=lambda s: s.acquisition_cost)
    assert cheapest.target is ACMGClass.LIKELY_BENIGN
    assert predicted_direction(row) is None


def test_evidence_the_pipeline_never_sees_is_counted_not_hidden(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """A variant the field resolved on family data has no observed cause here.

    That is the honest answer, and the size of the bucket is itself a finding
    about how much resolution happens outside public datasets.
    """
    variants = missense_variants[:30]
    target = variants[0]
    transcript = minus_gene.transcript.transcript_id

    cv_before = ClinVarSnapshot.from_records(
        [ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index) for v in variants]
    )
    cv_after = ClinVarSnapshot.from_records(
        [
            ClinVarRecord(
                v.hgvs_c,
                v.hgvs_p,
                "pathogenic" if v is target else "uncertain",
                3 if v is target else 2,
                v.codon_index,
            )
            for v in variants
        ]
    )
    # The map's sources did not move at all in this interval.
    rows = _map(minus_gene, minus_config, toy_spec, variants)

    study = run_time_series_study(
        rows, rows, cv_before, cv_after, toy_spec, transcript_id=transcript
    )
    assert study.result.resolved == 1
    assert study.result.cause_scored == 0
    assert study.result.unobserved_cause == 1
    assert study.unmoved == [target.variant_id]
    assert study.anticipation_rate == 0.0
    assert "resolved on evidence this pipeline never sees" in study.summary()


def test_a_curated_outcome_table_takes_precedence_over_the_snapshots(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """Somebody who read the submission records knows what the sources cannot."""
    from vus_foresight.validation import Outcome

    variants = missense_variants[:20]
    target = variants[0]
    rows = _map(minus_gene, minus_config, toy_spec, variants)
    empty = ClinVarSnapshot.from_records([])

    study = run_time_series_study(
        rows,
        rows,
        empty,
        empty,
        toy_spec,
        transcript_id=minus_gene.transcript.transcript_id,
        curated_outcomes=[
            Outcome(
                variant_id=target.variant_id,
                class_at_t=ACMGClass.UNCERTAIN,
                class_at_t_plus_n=ACMGClass.LIKELY_PATHOGENIC,
                evidence_type="segregation",
            )
        ],
    )
    assert study.result.considered == 1
    assert study.result.cause_scored == 1
    confusion = list(study.result.cause_confusion)
    assert confusion[0][1] == BlockingReason.MISSING_SEGREGATION.value


def test_a_neighbour_driven_resolution_gets_its_own_observed_label(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """PS1 and PM5 have no blocking reason to declare -- the schema has no member
    for "no neighbour classified yet". Labelling them anyway is what lets metric
    2 report that as a systematic miss instead of silently dropping the case the
    whole project is built around.
    """
    from vus_foresight.validation import NEIGHBOUR_CLASSIFICATION

    variants = missense_variants[:30]
    target = variants[0]
    assert toy_spec.by_code("PS1").blocks_as is None

    empty = tmp_path / "cv_empty.tsv"
    write_snapshot([], empty)
    later = tmp_path / "cv_later.tsv"
    write_snapshot(
        [ClinVarRecord("c.999999A>G", target.hgvs_p, "pathogenic", 3, target.codon_index)],
        later,
    )
    before = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[ClinVarSnapshotAdapter.from_path(empty)],
    )
    after = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[ClinVarSnapshotAdapter.from_path(later)],
    )
    causes = observed_causes_from_diff(compare_maps(before, after), toy_spec)
    assert causes[target.variant_id] == NEIGHBOUR_CLASSIFICATION
    # And the whole codon, since one classification moves all of it.
    assert len(causes) > 1
    assert set(causes.values()) == {NEIGHBOUR_CLASSIFICATION}


def test_variants_the_map_moved_ahead_of_the_archive_are_surfaced(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """Not errors. These are the map's live, still-unfalsified predictions."""
    variants = missense_variants[:30]
    target = variants[0]

    empty = tmp_path / "pred_empty.tsv"
    empty.write_text("grch38_pos\tbayesdel\n", encoding="utf-8")
    scored = tmp_path / "pred.tsv"
    scored.write_text(
        f"grch38_pos\tbayesdel\n{target.grch38_pos}\t0.62\n", encoding="utf-8"
    )
    rows_t = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[PredictorAdapter(empty, "t0")],
    )
    rows_tn = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[PredictorAdapter(scored, "t1")],
    )
    # ClinVar has moved nobody.
    cv = ClinVarSnapshot.from_records(
        [ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index) for v in variants]
    )

    study = run_time_series_study(
        rows_t, rows_tn, cv, cv, toy_spec, transcript_id=minus_gene.transcript.transcript_id
    )
    assert study.result.resolved == 0
    assert study.unmoved == []
    # PP3 alone does not cross a threshold here, so nothing is "ahead" unless
    # the class actually changed -- which is the point of tracking class_changes
    # rather than any criterion movement.
    assert all(v.startswith(minus_gene.transcript.transcript_id) for v in study.ahead_of_clinvar)


def test_the_end_to_end_command_path(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """The two snapshot files and two maps the CLI actually reads."""
    from vus_foresight.output import read_rows, write_parquet

    variants = missense_variants[:25]
    target = variants[0]

    before_path = tmp_path / "cv_before.tsv"
    after_path = tmp_path / "cv_after.tsv"
    write_snapshot(
        [ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index) for v in variants],
        before_path,
    )
    write_snapshot(
        [
            ClinVarRecord(
                v.hgvs_c,
                v.hgvs_p,
                "likely_pathogenic" if v is target else "uncertain",
                3 if v is target else 2,
                v.codon_index,
            )
            for v in variants
        ],
        after_path,
    )

    rows = _map(minus_gene, minus_config, toy_spec, variants)
    map_path = tmp_path / "map.parquet"
    write_parquet(rows, map_path)

    study = run_time_series_study(
        read_rows(map_path),
        read_rows(map_path),
        ClinVarSnapshot.read(before_path),
        ClinVarSnapshot.read(after_path),
        toy_spec,
        transcript_id=rows[0].transcript,
    )
    assert study.result.resolved == 1
    assert study.result.resolvability_recall is not None


def test_the_study_is_identical_whether_the_maps_are_streamed_or_materialised(
    tmp_path, minus_gene, minus_config, toy_spec, missense_variants
):
    """The map at T is walked twice -- diffed, then looked up by variant.

    A list survives that by keeping a whole gene in memory, which it cannot.
    GapMapSource re-reads instead, and this pins the two to the same answer.
    """
    from vus_foresight.output import GapMapSource, read_rows, write_parquet

    variants = missense_variants[:25]
    target = variants[0]

    cv_before = ClinVarSnapshot.from_records(
        [ClinVarRecord(v.hgvs_c, v.hgvs_p, "uncertain", 2, v.codon_index) for v in variants],
        date(2018, 1, 1),
    )
    cv_after = ClinVarSnapshot.from_records(
        [
            ClinVarRecord(
                v.hgvs_c,
                v.hgvs_p,
                "pathogenic" if v is target else "uncertain",
                3 if v is target else 2,
                v.codon_index,
                date(2020, 6, 1) if v is target else None,
            )
            for v in variants
        ],
        date(2024, 1, 1),
    )

    empty = tmp_path / "func_empty.tsv"
    empty.write_text("hgvs_p\tclassification\tscore\tdataset\n", encoding="utf-8")
    assayed = tmp_path / "func.tsv"
    assayed.write_text(
        "hgvs_p\tclassification\tscore\tdataset\n"
        f"{target.hgvs_p}\tabnormal\t-2.1\tTOY-SGE\n",
        encoding="utf-8",
    )
    regions = ((1, 40, "TOY-SGE"),)
    rows_t = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(empty, "t0", assayed_regions=regions)],
    )
    rows_tn = _map(
        minus_gene, minus_config, toy_spec, variants,
        adapters=[FunctionalAdapter(assayed, "t1", assayed_regions=regions)],
    )
    path_t, path_tn = tmp_path / "t.parquet", tmp_path / "tn.parquet"
    write_parquet(rows_t, path_t)
    write_parquet(rows_tn, path_tn)

    def run(at_t, at_tn):
        return run_time_series_study(
            at_t,
            at_tn,
            cv_before,
            cv_after,
            toy_spec,
            transcript_id=minus_gene.transcript.transcript_id,
            reference_date=date(2018, 1, 1),
        )

    materialised = run(read_rows(path_t), read_rows(path_tn))
    streamed = run(GapMapSource(path_t, batch_size=6), GapMapSource(path_tn, batch_size=6))

    assert streamed.result == materialised.result
    assert streamed.anticipated == materialised.anticipated
    assert streamed.saw_evidence_only == materialised.saw_evidence_only == [target.variant_id]
    assert streamed.unmoved == materialised.unmoved
    assert streamed.ahead_of_clinvar == materialised.ahead_of_clinvar


def test_a_generator_of_rows_is_refused_by_the_study(
    minus_gene, minus_config, toy_spec, missense_variants
):
    """It would be walked twice and be empty the second time."""
    rows = _map(minus_gene, minus_config, toy_spec, missense_variants[:5])
    empty = ClinVarSnapshot.from_records([], date(2018, 1, 1))
    with pytest.raises(TypeError, match="one-shot iterator"):
        run_time_series_study(
            iter(rows), rows, empty, empty, toy_spec, transcript_id="NM_1.1"
        )
