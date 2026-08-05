"""The edges of the section 10 protocol: the inputs a real study actually hits.

The oracle of every test here is a *hand-computed* metric or an exact message,
never a value read back from the code. Metric 2 of the protocol is the
publishable one, and a validation module that silently miscounts is worse than
no validation at all: it produces a number that looks like a result. So each
test states the count it expects and why that count is the right one.

How these tests could be wrong, and what is done about it:

* **The rows could be unrepresentative.** Every ``GapMapRow`` used here comes
  out of the real ``MapRunner`` over the synthetic gene -- none is hand-built --
  so a row shape that the engine cannot produce cannot be asserted about by
  accident. Where an edge needs a row the default fixtures never yield (no
  feasible evidence set, a benign lean, a gap already closed), the *rules* are
  varied rather than the rows: the repository's own toy specification is loaded,
  mutated as data, and re-loaded through ``load_spec``.
* **A spec mutation could pass vacuously.** Two of the three cases below expect
  *no* observed cause, which is also what a broken fixture would produce. The
  third spells the surviving label out in full, so if the toy specification ever
  stops mapping PS3 onto ``MISSING_FUNCTIONAL`` that case fails loudly and the
  premise of the other two is not silently lost.
* **A branch could be reached without being tested.** Reaching a line and
  distinguishing its two arms are different things, and the guards below are the
  cases where they came apart: a closed gap only proves anything on a row whose
  sufficient sets are all intractable, since a feasible row is resolvable either
  way; and an absent timing gap only proves anything in a cohort of at least
  three, since ``spearman`` returns ``None`` below that whatever the code did
  with it. Each such test states the value the wrong behaviour would produce.
* **A negative could pass because the premise collapsed.** Every "nothing came
  out" assertion over the snapshots carries a second variant that *does* come
  out, so a lost star rating or a mangled classification column fails the test
  instead of satisfying it.
* **A metric could be checked only where it is defined.** Half of this module is
  the undefined cases -- no resolutions, no scored directions, no observed
  causes -- because "n/a" and ``0.0`` are different claims about a study, and
  only one of them is honest when the denominator is empty.

The distinction between ``observed_causes=None`` and ``observed_causes={}`` gets
its own test: they select different arms of the same branch, and collapsing them
would report zero unobserved causes exactly when every cause is unobserved.
"""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

import pytest
import yaml

from vus_foresight.acmg import ACMGClass, BlockingReason, FeasibilityTag
from vus_foresight.adapters import AdapterRegistry
from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from vus_foresight.adapters.clinvar import ClinVarRecord, ClinVarSnapshot
from vus_foresight.adapters.clinvar_import import write_snapshot
from vus_foresight.adapters.tabular import FunctionalAdapter, PredictorAdapter
from vus_foresight.engine.pipeline import MapRunner
from vus_foresight.engine.spec import load_spec
from vus_foresight.engine.timeline import compare_maps
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.validation import (
    Outcome,
    TimeSeriesStudy,
    ValidationResult,
    is_resolvable,
    observed_causes_from_diff,
    outcomes_from_snapshots,
    predicted_direction,
    read_outcomes,
    validate,
)
from vus_foresight.variant import Consequence

TOY_SPEC_PATH = Path(__file__).resolve().parents[2] / "config" / "specs" / "toy_v0.1.0.yaml"

#: The assayed interval the functional fixture below declares, in codons.
ASSAYED_REGIONS = ((1, 40, "TOY-SGE"),)


def _map(gene, config, spec, variants, *, adapters=()):
    registry = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(config)])
    for adapter in adapters:
        registry.add(adapter)
    runner = MapRunner(
        transcript=gene.transcript,
        gene=config,
        spec=spec,
        adapters=registry,
        computed_at=datetime(1970, 1, 1),
    )
    return [result.row for result in runner.run(variants)]


def _spec_variant(destination: Path, mutate):
    """Load the repository's toy specification with one thing changed.

    Rule-as-data is the project's own claim, so the way to reach a rule edge is
    to vary the rules through the real loader rather than to fake a spec object.
    """
    raw = yaml.safe_load(TOY_SPEC_PATH.read_text(encoding="utf-8"))
    mutate(raw)
    destination.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return load_spec(destination)


@pytest.fixture(scope="module")
def missense_variants(minus_gene):
    return [
        v
        for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE and v.codon_index and v.codon_index > 3
    ]


def _make_everything_intractable(raw):
    """No criterion in this specification can be obtained by anybody.

    Not a contrived state: it is what a gene with no assay, no cohort and no
    predictor coverage looks like. It is also the way to take *feasibility* out
    of an answer, so that a test about a closed gap cannot pass through the
    feasible-set fallback it is not testing.
    """
    for criterion in raw["criteria"]:
        criterion["feasibility"] = "intractable"
        criterion.pop("feasibility_ladder", None)
    # Keep the sets in the output so the "no feasible route" case is the one
    # being tested, rather than the empty-set case that hides it.
    raw["gap"]["include_intractable"] = True


def _nonsense_variants(gene):
    return [
        v for v in enumerate_coding_snvs(gene.transcript) if v.consequence is Consequence.NONSENSE
    ][:6]


def _benign_predictor_table(directory: Path, variant) -> Path:
    table = directory / "bayesdel_low.tsv"
    table.write_text(f"grch38_pos\tbayesdel\n{variant.grch38_pos}\t0.05\n", encoding="utf-8")
    return table


@pytest.fixture(scope="module")
def predictor_map(tmp_path_factory, minus_gene, minus_config, toy_spec, missense_variants):
    """A real map in which exactly one variant carries a benign predictor call.

    The whole map is returned, not just that row: metric 4 needs several
    variants with *different* gaps before a rank correlation exists at all.
    """
    target = missense_variants[0]
    table = _benign_predictor_table(tmp_path_factory.mktemp("predictor"), target)
    rows = _map(
        minus_gene,
        minus_config,
        toy_spec,
        missense_variants[:20],
        adapters=[PredictorAdapter(table, "t")],
    )
    return rows, target.hgvs_c


@pytest.fixture(scope="module")
def leaning_benign_row(predictor_map):
    """A real row on which the map leans benign: BP4 fired and nothing else did.

    Its ``gap_to_LB`` is ``None`` -- one supporting benign point already reaches
    the toy threshold -- which makes it the row for both "which way did the map
    lean" and "a gap already closed contributes no timing point".
    """
    rows, hgvs_c = predictor_map
    return next(row for row in rows if row.hgvs_c == hgvs_c)


@pytest.fixture(scope="module")
def rows_past_a_verdict(tmp_path_factory, minus_gene, minus_config, missense_variants):
    """One row past LP and one past LB, under a spec where nothing is obtainable.

    Both directions are needed and the intractable specification is the point:
    with every sufficient set intractable, the feasible-set fallback answers
    *no*, so anything these rows still get scored as resolvable comes from the
    closed gap and from nothing else.
    """
    directory = tmp_path_factory.mktemp("past-verdict")
    spec = _spec_variant(directory / "intractable.yaml", _make_everything_intractable)
    target = missense_variants[0]
    table = _benign_predictor_table(directory, target)

    leaning = _map(
        minus_gene,
        minus_config,
        spec,
        missense_variants[:12],
        adapters=[PredictorAdapter(table, "t")],
    )
    truncating = _map(minus_gene, minus_config, spec, _nonsense_variants(minus_gene))
    return {
        "gap_to_LB": next(row for row in leaning if row.hgvs_c == target.hgvs_c),
        "gap_to_LP": next(row for row in truncating if row.gap_to_LP is None),
    }


@pytest.fixture(scope="module")
def no_feasible_route_rows(tmp_path_factory, minus_gene, minus_config, missense_variants):
    """A map computed under a specification whose every criterion is intractable.

    The only state in which the map reports a finite gap it cannot tell anyone
    how to close.
    """
    directory = tmp_path_factory.mktemp("intractable")
    spec = _spec_variant(directory / "intractable.yaml", _make_everything_intractable)
    return _map(minus_gene, minus_config, spec, missense_variants[:12])


@pytest.fixture(scope="module")
def assay_diff(tmp_path_factory, minus_gene, minus_config, toy_spec, missense_variants):
    """A real diff in which exactly one variant gained PS3 between T and T+n."""
    directory = tmp_path_factory.mktemp("functional")
    variants = missense_variants[:30]
    target = variants[0]

    empty = directory / "functional_t0.tsv"
    empty.write_text("hgvs_p\tclassification\tscore\tdataset\n", encoding="utf-8")
    assayed = directory / "functional_t1.tsv"
    assayed.write_text(
        f"hgvs_p\tclassification\tscore\tdataset\n{target.hgvs_p}\tabnormal\t-2.1\tTOY-SGE\n",
        encoding="utf-8",
    )

    before = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        adapters=[FunctionalAdapter(empty, "t0", assayed_regions=ASSAYED_REGIONS)],
    )
    after = _map(
        minus_gene,
        minus_config,
        toy_spec,
        variants,
        adapters=[FunctionalAdapter(assayed, "t1", assayed_regions=ASSAYED_REGIONS)],
    )
    return compare_maps(before, after), target.variant_id


# --------------------------------------------------------------------------
# read_outcomes: the curated table a human hands the study.
# --------------------------------------------------------------------------


def test_an_outcome_table_with_only_the_required_columns_records_no_evidence_claim(tmp_path):
    """Absent columns must read as "unrecorded", never as a default cause.

    An outcome that silently acquired an ``evidence_type`` would be scored by
    metric 2 against a claim nobody made.
    """
    path = tmp_path / "outcomes.tsv"
    path.write_text(
        "variant_id\tclass_at_t\tclass_at_t_plus_n\nNM_1.1:c.10A>G\tVUS\tLP\n",
        encoding="utf-8",
    )

    outcomes = read_outcomes(path)

    assert outcomes == [
        Outcome(
            variant_id="NM_1.1:c.10A>G",
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.LIKELY_PATHOGENIC,
            evidence_type="",
            resolved_on=None,
        )
    ]


def test_the_optional_columns_are_read_and_the_evidence_type_is_trimmed(tmp_path):
    """``evidence_type`` is a key into the mapping, so stray whitespace in a
    hand-maintained table would silently drop the outcome out of metric 2."""
    path = tmp_path / "outcomes.tsv"
    path.write_text(
        "variant_id\tclass_at_t\tclass_at_t_plus_n\tevidence_type\tresolved_on\n"
        "NM_1.1:c.10A>G\tVUS\tB\t  functional \t2021-03-04\n",
        encoding="utf-8",
    )

    outcomes = read_outcomes(path)

    assert outcomes == [
        Outcome(
            variant_id="NM_1.1:c.10A>G",
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.BENIGN,
            evidence_type="functional",
            resolved_on=date(2021, 3, 4),
        )
    ]


def test_a_row_that_stops_short_of_the_header_leaves_the_optional_fields_empty(tmp_path):
    """A ragged TSV is the normal state of a hand-curated table.

    The short row must land in the same state as an absent column rather than
    carrying a ``None`` through into ``evidence_type`` and blowing up in the
    mapping lookup later, where the cause would be much harder to see.
    """
    path = tmp_path / "outcomes.tsv"
    path.write_text(
        "variant_id\tclass_at_t\tclass_at_t_plus_n\tevidence_type\tresolved_on\n"
        "NM_1.1:c.10A>G\tVUS\tLB\n",
        encoding="utf-8",
    )

    outcomes = read_outcomes(path)

    assert outcomes == [
        Outcome(
            variant_id="NM_1.1:c.10A>G",
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.LIKELY_BENIGN,
            evidence_type="",
            resolved_on=None,
        )
    ]


def test_a_table_with_a_header_and_no_rows_yields_no_outcomes(tmp_path):
    path = tmp_path / "outcomes.tsv"
    path.write_text("variant_id\tclass_at_t\tclass_at_t_plus_n\n", encoding="utf-8")

    assert read_outcomes(path) == []


@pytest.mark.parametrize(
    ("content", "error", "message"),
    [
        pytest.param(
            "class_at_t\tclass_at_t_plus_n\nVUS\tLP\n",
            KeyError,
            "variant_id",
            id="identifier-column-missing",
        ),
        pytest.param(
            "variant_id\tclass_at_t\tclass_at_t_plus_n\nNM_1.1:c.10A>G\tVUS\tPROBABLY\n",
            ValueError,
            r"'PROBABLY' is not a valid ACMGClass",
            id="class-outside-the-five-tier-vocabulary",
        ),
        pytest.param(
            "variant_id\tclass_at_t\tclass_at_t_plus_n\tresolved_on\n"
            "NM_1.1:c.10A>G\tVUS\tLP\t04/03/2021\n",
            ValueError,
            r"Invalid isoformat string: '04/03/2021'",
            id="date-not-in-iso-form",
        ),
    ],
)
def test_a_malformed_outcome_row_is_refused_rather_than_guessed(tmp_path, content, error, message):
    """Every one of these could be coerced into something. None should be.

    The outcome table is the study's ground truth; a row nobody can read is a
    row nobody should score, and guessing at it would corrupt all four metrics
    at once with no trace in the output.
    """
    path = tmp_path / "outcomes.tsv"
    path.write_text(content, encoding="utf-8")

    with pytest.raises(error, match=message):
        read_outcomes(path)


# --------------------------------------------------------------------------
# Outcome: which way the field moved.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "class_at_t_plus_n", [ACMGClass.BENIGN, ACMGClass.LIKELY_BENIGN], ids=["B", "LB"]
)
def test_a_resolution_towards_either_benign_tier_reports_one_direction(class_at_t_plus_n):
    """Metric 3 asks about direction, not tier, so B and LB are one answer.

    Keeping them apart here would split the denominator of a metric whose whole
    question is "did the map point the right way".
    """
    outcome = Outcome(
        variant_id="NM_1.1:c.10A>G",
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=class_at_t_plus_n,
    )

    assert outcome.direction == "benign"
    assert outcome.was_resolved is True


# --------------------------------------------------------------------------
# The metrics when their denominators are empty.
# --------------------------------------------------------------------------


def test_a_study_that_resolved_nothing_reports_undefined_rather_than_zero():
    """An empty denominator is not a score of zero.

    ``0.0`` says the map got everything wrong; ``None`` says nothing was asked.
    Reporting the first for the second is how a study with no data becomes a
    published claim that the map failed.
    """
    result = ValidationResult(considered=7)

    assert result.resolvability_recall is None
    assert result.cause_accuracy is None
    assert result.direction_accuracy is None


def test_the_summary_of_an_empty_study_prints_n_a_for_every_metric():
    """The summary is what a reader sees; the undefined cases must survive it."""
    study = TimeSeriesStudy(result=ValidationResult())

    lines = study.summary().splitlines()

    assert lines[2] == "1. resolvability recall  n/a"
    assert lines[5] == "4. temporal calibration  n/a"
    assert lines[6] == "5. map reached a verdict n/a (0 of 0)"


def test_a_time_series_with_no_resolutions_leaves_both_of_its_rates_undefined():
    """Verdict rate and visibility rate answer different questions, and both are
    unanswerable before the archive has resolved anything."""
    study = TimeSeriesStudy(result=ValidationResult(), ahead_of_clinvar=["NM_1.1:c.10A>G"])

    assert study.resolved_total == 0
    assert study.anticipation_rate is None
    assert study.evidence_visibility_rate is None


# --------------------------------------------------------------------------
# Ground truth from the snapshots: what may not become an outcome.
# --------------------------------------------------------------------------


def test_a_variant_that_left_the_archive_is_not_counted_as_a_resolution():
    """A withdrawn record is a variant nobody classified, not a variant resolved.

    Counting the disappearance would inflate metric 1's denominator with cases
    the map was never given a chance to be right about.

    The second variant is present in both snapshots and carries the test: an
    empty result would otherwise be indistinguishable from a fixture in which
    nothing was eligible in the first place, and the withdrawal would go
    unmeasured. What is asserted is that exactly one of the two survives.
    """
    before = ClinVarSnapshot.from_records(
        [
            ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 2, 4),
            ClinVarRecord("c.11C>T", "p.Thr4Ile", "uncertain", 2, 4),
        ]
    )
    after = ClinVarSnapshot.from_records(
        [ClinVarRecord("c.11C>T", "p.Thr4Ile", "benign", 2, 4, date(2024, 5, 1))]
    )

    assert outcomes_from_snapshots(before, after, transcript_id="NM_1.1") == [
        Outcome(
            variant_id="NM_1.1:c.11C>T",
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.BENIGN,
            evidence_type="",
            resolved_on=date(2024, 5, 1),
        )
    ]


def test_a_classification_the_study_cannot_map_is_dropped_rather_than_guessed(tmp_path):
    """Snapshot files carry the classification column verbatim.

    ``ClinVarSnapshot.read`` does not re-normalise, so a term this module has no
    entry for -- one a future build step learns before the study does -- reaches
    here as a raw string. Dropping it keeps an unrecognised word out of the
    numerator of every metric; the alternative is scoring against a class nobody
    has decided the meaning of.

    Only the first of the two variants is unreadable. The second is the control
    that keeps the assertion honest: this snapshot goes through a real
    ``write_snapshot``/``read`` round trip, so an empty result on its own would
    also be what a lost star rating or a mangled classification column looks
    like, and the drop would never be attributable to the term.
    """
    before_path, after_path = tmp_path / "cv_t.tsv", tmp_path / "cv_tn.tsv"
    write_snapshot(
        [
            ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 2, 4),
            ClinVarRecord("c.11C>T", "p.Thr4Ile", "uncertain", 2, 4),
        ],
        before_path,
    )
    write_snapshot(
        [
            ClinVarRecord("c.10A>G", "p.Thr4Ala", "drug_response", 2, 4),
            ClinVarRecord("c.11C>T", "p.Thr4Ile", "likely_benign", 2, 4, date(2024, 5, 1)),
        ],
        after_path,
    )

    outcomes = outcomes_from_snapshots(
        ClinVarSnapshot.read(before_path),
        ClinVarSnapshot.read(after_path),
        transcript_id="NM_1.1",
    )

    assert outcomes == [
        Outcome(
            variant_id="NM_1.1:c.11C>T",
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.LIKELY_BENIGN,
            evidence_type="",
            resolved_on=date(2024, 5, 1),
        )
    ]


# --------------------------------------------------------------------------
# Metric 2 when the specification and the diff disagree about the vocabulary.
# --------------------------------------------------------------------------


def _drop_ps3(raw):
    """The map at T was computed with a criterion this specification no longer has."""
    raw["criteria"] = [c for c in raw["criteria"] if c["code"] != "PS3"]


def _drop_ps3_blocks_as(raw):
    """PS3 stays, but stops declaring which block its absence causes."""
    for criterion in raw["criteria"]:
        if criterion["code"] == "PS3":
            del criterion["blocks_as"]


def _shorten_blocking_priority(raw):
    """A priority list that does not mention a reason the criteria still declare."""
    raw["blocking"]["priority"] = ["MISSING_SEGREGATION"]


@pytest.mark.parametrize(
    ("mutation", "expected_label"),
    [
        pytest.param(_drop_ps3, None, id="criterion-gone-from-the-spec"),
        pytest.param(_drop_ps3_blocks_as, None, id="criterion-declares-no-blocks-as"),
        pytest.param(
            _shorten_blocking_priority,
            BlockingReason.MISSING_FUNCTIONAL.value,
            id="label-absent-from-the-blocking-priority",
        ),
    ],
)
def test_an_observed_cause_survives_only_while_the_spec_can_explain_it(
    tmp_path, assay_diff, mutation, expected_label
):
    """The diff is fixed; the specification reading it is varied.

    Two of these are the honest silence the module documents -- a criterion the
    spec no longer has, and a criterion that declares no block -- and both must
    yield *no* cause rather than a plausible-looking one, because metric 2 is
    the publishable metric and a fabricated attribution is unfalsifiable.

    The third is the opposite risk: a specification whose ``priority`` list does
    not mention every reason its own criteria declare. There the label is known
    and only its ranking is missing, so dropping it would lose a real
    observation to a bookkeeping omission.
    """
    diff, target_id = assay_diff
    spec = _spec_variant(tmp_path / "variant.yaml", mutation)

    causes = observed_causes_from_diff(diff, spec)

    assert causes == ({} if expected_label is None else {target_id: expected_label})


# --------------------------------------------------------------------------
# What the map said, read back by the study.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "closed_gap", ["gap_to_LP", "gap_to_LB"], ids=["already-past-LP", "already-past-LB"]
)
def test_a_variant_already_past_a_verdict_counts_as_resolvable(rows_past_a_verdict, closed_gap):
    """No gap in one direction means nothing is left to obtain.

    Falling through to the sufficient-set check here would mark a variant the
    map had *already* resolved as unresolvable, which is metric 1 counting its
    own successes as failures.

    Both rows are computed under the specification in which every sufficient set
    is intractable, and that is what makes the answer mean something: the
    fallback says *not resolvable* for these rows, so ``True`` can only have
    come from the closed gap. A row with feasible sets would have been scored
    resolvable either way and would have tested nothing.
    """
    row = rows_past_a_verdict[closed_gap]

    assert getattr(row, closed_gap) is None
    assert {s.feasibility for s in row.minimum_sufficient_sets} == {FeasibilityTag.INTRACTABLE}
    assert is_resolvable(row) is True


def test_a_finite_gap_nobody_can_close_does_not_count_as_resolvable(no_feasible_route_rows):
    """A number is not a route.

    Every sufficient set here is intractable, so the map can say how many points
    are missing and still have no answer to "who would obtain them". Metric 1
    asks the second question.
    """
    row = next(r for r in no_feasible_route_rows if r.minimum_sufficient_sets)

    assert {s.feasibility for s in row.minimum_sufficient_sets} == {FeasibilityTag.INTRACTABLE}
    assert is_resolvable(row) is False


def test_the_map_leans_benign_when_its_accumulated_points_are_negative(leaning_benign_row):
    """The sign of the points, and nothing else.

    BP4 alone is one supporting benign point. That is a lean, and it is the only
    thing metric 3 is entitled to read.
    """
    assert leaning_benign_row.points_current == -1
    assert predicted_direction(leaning_benign_row) == "benign"


# --------------------------------------------------------------------------
# validate: the bookkeeping of the four metrics.
# --------------------------------------------------------------------------


def test_an_outcome_with_no_row_in_the_map_is_not_considered(leaning_benign_row):
    """The map is one gene; an outcome table need not be.

    A variant the map never enumerated cannot be scored, and counting it as
    considered would let a mismatched pair of inputs quietly depress metric 1.
    """
    outcome = Outcome(
        variant_id="NM_999002.1:c.999999A>G",
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.PATHOGENIC,
    )

    result = validate([leaning_benign_row], [outcome])

    assert result.considered == 0
    assert result.resolved == 0
    assert result.misses == []


def test_a_resolution_the_map_saw_no_way_to_close_is_written_down_as_a_miss(
    no_feasible_route_rows,
):
    """Metric 1's failures are listed by name, never summarised away.

    The map claimed there was no feasible evidence and the field resolved the
    variant anyway: that is the one bucket a reviewer has to read case by case,
    so the identifier and the predicted blocking reason both have to survive
    into the message.
    """
    row = next(
        r
        for r in no_feasible_route_rows
        if not is_resolvable(r) and r.blocking_reason is BlockingReason.MISSING_FUNCTIONAL
    )
    outcome = Outcome(
        variant_id=row.variant_id,
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.LIKELY_PATHOGENIC,
    )

    result = validate(no_feasible_route_rows, [outcome])

    assert result.resolved == 1
    assert result.resolvability_recall == 0.0
    assert result.misses == [
        f"{row.variant_id}: resolved to LP but the map offered no feasible "
        "evidence set (blocking_reason=MISSING_FUNCTIONAL)"
    ]


def test_no_time_series_at_all_is_not_the_same_as_a_time_series_that_moved_nothing(
    leaning_benign_row,
):
    """``None`` and ``{}`` are different claims and must stay apart.

    Without a time series there is nothing to say about where the evidence came
    from, so ``unobserved_cause`` must stay at zero. With a time series in which
    nothing moved, every resolution happened on evidence this pipeline never
    saw, and that count is the finding. Collapsing the two reports zero
    unobserved causes precisely when every cause is unobserved.
    """
    outcome = Outcome(
        variant_id=leaning_benign_row.variant_id,
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.BENIGN,
    )

    without_series = validate([leaning_benign_row], [outcome])
    empty_series = validate([leaning_benign_row], [outcome], observed_causes={})

    assert (without_series.resolved, without_series.unobserved_cause) == (1, 0)
    assert (empty_series.resolved, empty_series.unobserved_cause) == (1, 1)


def test_an_evidence_type_outside_the_mapping_scores_nothing_and_shadows_the_series(
    leaning_benign_row,
):
    """Characterisation, and a sharp edge worth knowing about.

    A curated ``evidence_type`` wins over the time series by design -- a human
    who read the submission record knows more than the pipeline's sources. But
    the curated value wins *before* it is looked up, so a typo or a vocabulary
    the mapping has not learned discards a cause the time series did observe,
    and the outcome lands in ``unobserved_cause`` instead of being scored.
    """
    outcome = Outcome(
        variant_id=leaning_benign_row.variant_id,
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.BENIGN,
        evidence_type="astrology",
    )

    result = validate(
        [leaning_benign_row],
        [outcome],
        observed_causes={leaning_benign_row.variant_id: BlockingReason.MISSING_FUNCTIONAL.value},
    )

    assert result.cause_scored == 0
    assert result.cause_confusion == {}
    assert result.unobserved_cause == 1


def test_a_lean_in_the_wrong_direction_is_scored_wrong_and_not_quietly_dropped(
    leaning_benign_row,
):
    """Metric 3's denominator is every claim the map made, right or wrong.

    The map leant benign and the field went pathogenic. Skipping the case would
    leave metric 3 measuring only the predictions that happened to come true.
    """
    outcome = Outcome(
        variant_id=leaning_benign_row.variant_id,
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.PATHOGENIC,
    )

    result = validate([leaning_benign_row], [outcome])

    assert result.direction_scored == 1
    assert result.direction_correct == 0
    assert result.direction_accuracy == 0.0
    assert result.direction_unpredicted == 0


def test_a_gap_already_closed_contributes_no_point_to_the_timing_correlation(predictor_map):
    """Metric 4 correlates gap size against delay, and there is no gap here.

    The cohort is built so that the answer distinguishes *skipped* from
    *counted as zero*, which a single-variant cohort cannot do -- ``spearman``
    returns ``None`` below three points, so any one-row study reports ``None``
    whatever the code does with the absent gap.

    Three variants carry a real gap and are resolved in the exact reverse order
    of its size, a rank correlation of exactly ``-1``:

    ====================  ====  =====
    gap the metric reads  gap   days
    ====================  ====  =====
    ``gap_to_LP`` (P)     6      100
    ``gap_to_LP`` (P)     4      500
    ``gap_to_LB`` (B)     3     1000
    ====================  ====  =====

    The fourth is the row whose ``gap_to_LB`` is ``None`` because the map had
    already reached LB: no distance, so no point. Were that absent gap read as
    zero, its 300-day delay would land out of rank order and rho would fall to
    ``-0.4``. So ``-1`` is the assertion that the made-up point stayed out.
    """
    rows, _ = predictor_map
    reference = date(2018, 1, 1)
    unevidenced = [row for row in rows if row.gap_to_LP == 6 and row.gap_to_LB == 1]
    leaning_pathogenic = [row for row in rows if row.gap_to_LP == 4 and row.gap_to_LB == 3]
    closed = next(row for row in rows if row.gap_to_LB is None)

    def outcome(row, tier, resolved_on):
        return Outcome(
            variant_id=row.variant_id,
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=tier,
            resolved_on=resolved_on,
        )

    result = validate(
        rows,
        [
            outcome(unevidenced[0], ACMGClass.PATHOGENIC, date(2018, 4, 11)),  # gap 6, 100 days
            outcome(leaning_pathogenic[0], ACMGClass.PATHOGENIC, date(2019, 5, 16)),  # gap 4, 500
            outcome(leaning_pathogenic[1], ACMGClass.BENIGN, date(2020, 9, 27)),  # gap 3, 1000
            outcome(closed, ACMGClass.LIKELY_BENIGN, date(2018, 10, 28)),  # no gap, 300 days
        ],
        reference_date=reference,
    )

    assert closed.gap_to_LB is None
    assert result.resolved == 4
    assert result.temporal_correlation == pytest.approx(-1.0)
