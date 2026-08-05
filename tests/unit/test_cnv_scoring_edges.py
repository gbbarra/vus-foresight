"""Edges of the copy-number framework: the scale that must never be mixed in.

Oracle: the ClinGen CNV framework (Riggs et al. 2020) as transcribed into
:class:`~vus_foresight.engine.cnv_scoring.CNVScoringConfig`, read as a table of
concrete numbers. Every expectation below is a literal a curator could check
against that table -- the score in hundredths, the section code, the sentence
of rationale, the resulting class.

How these tests could be *wrong*: the framework's numbers could have been
transcribed wrongly here, in the same way they could have been transcribed
wrongly into the module. They are not protected against that. What they *are*
protected against is the failure this module exists to prevent, which is a
number on the CNV scale being read as a number on the Tavtigian scale: every
assertion pins both the raw integer and the class it produces, so a change to
one that is not matched by the other shows up as a failure rather than as an
aggregation that silently means something else.

The other deliberate choice: the contrasts are always tested as *pairs* on one
axis (deletion vs duplication, whole-gene vs intragenic, frame-preserving vs
frame-disrupting, established vs unestablished mechanism), because a
copy-and-paste error between two arms of the branch tree is invisible when each
arm is asserted alone.
"""

from __future__ import annotations

from datetime import date, datetime

import pytest

from vus_foresight.acmg import ACMGClass, BlockingReason, Direction, EvidenceClass, Strength
from vus_foresight.engine.cnv_scoring import CNVScore, CNVScoringConfig, cnv_row, score_cnv
from vus_foresight.enumeration import enumerate_exon_cnvs
from vus_foresight.testing import synthetic_gene_config
from vus_foresight.variant import Consequence, Variant, VariantKind

FIXED_TIME = datetime(1970, 1, 1)


def pick_cnv(transcript, *, first: str, last: str, consequence: Consequence) -> Variant:
    """The one enumerated CNV spanning exons ``first``..``last``.

    Selecting by exon label rather than by position keeps the test readable
    against the fixture's labels, which skip ``4`` on purpose.
    """
    return next(
        v
        for v in enumerate_exon_cnvs(transcript)
        if v.attributes["exon_first"] == first
        and v.attributes["exon_last"] == last
        and v.consequence is consequence
    )


# ---------------------------------------------------------------------------
# The scale itself: where each class begins.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("score", "expected"),
    [
        (250, ACMGClass.PATHOGENIC),
        (100, ACMGClass.PATHOGENIC),
        (99, ACMGClass.PATHOGENIC),
        (98, ACMGClass.LIKELY_PATHOGENIC),
        (90, ACMGClass.LIKELY_PATHOGENIC),
        (89, ACMGClass.UNCERTAIN),
        (0, ACMGClass.UNCERTAIN),
        (-89, ACMGClass.UNCERTAIN),
        (-90, ACMGClass.LIKELY_BENIGN),
        (-98, ACMGClass.LIKELY_BENIGN),
        (-99, ACMGClass.BENIGN),
        (-250, ACMGClass.BENIGN),
    ],
)
def test_each_cnv_class_begins_exactly_at_its_own_threshold(score, expected):
    """The thresholds are inclusive on both signs, and symmetric about zero.

    A one-hundredth slip at any of the four boundaries reclassifies a whole band
    of copy-number calls, so each boundary is pinned together with the value
    immediately outside it.
    """
    config = CNVScoringConfig()

    result = config.classify(score)

    assert result is expected


@pytest.mark.parametrize(
    ("score", "expected"),
    [(100, 1.0), (90, 0.9), (15, 0.15), (0, 0.0), (-99, -0.99)],
)
def test_a_stored_score_reads_back_as_the_framework_s_own_decimal(score, expected):
    """The integer column is hundredths; the framework speaks in units of 1.00.

    ``as_fraction`` is the only place the two representations meet, so it is the
    only place a factor-of-100 error can be introduced without touching a
    threshold.
    """
    scored = CNVScore(score=score, section="2E", rationale="irrelevant here")

    assert scored.as_fraction == expected


@pytest.mark.parametrize(
    ("kind", "consequence"),
    [
        (VariantKind.SNV, Consequence.MISSENSE),
        (VariantKind.FRAMESHIFT_CLASS, Consequence.FRAMESHIFT),
        (VariantKind.INFRAME_DELETION, Consequence.INFRAME_DELETION),
    ],
)
def test_a_sequence_variant_offered_to_the_cnv_scorer_is_refused_by_name(
    minus_config, kind, consequence
):
    """The two scales share one schema, so this guard is the only thing that

    keeps a point-scored variant from acquiring a copy-number score. It must
    name the offending variant, otherwise a whole-gene run reports an
    unlocatable failure. The frameshift class is included because it is also a
    class-level row with ``mutational_distance == 0``, i.e. the row that most
    resembles a CNV without being one.
    """
    variant = Variant(
        gene="TOYM",
        transcript="NM_999002.1",
        kind=kind,
        hgvs_c="c.100A>T",
        consequence=consequence,
        cds_position=100,
    )

    with pytest.raises(ValueError, match=r"NM_999002\.1:c\.100A>T is not a copy number variant"):
        score_cnv(variant, minus_config, CNVScoringConfig())


# ---------------------------------------------------------------------------
# Whole-gene loss: the gene-level fact, not the variant, decides.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mechanism", "expected_score", "expected_class"),
    [
        ("established", 100, ACMGClass.PATHOGENIC),
        ("not_established", 15, ACMGClass.UNCERTAIN),
        ("unknown", 15, ACMGClass.UNCERTAIN),
    ],
)
def test_a_whole_gene_deletion_is_pathogenic_only_where_loss_of_function_is_established(
    minus_gene, mechanism, expected_score, expected_class
):
    """Identical deletion, three genes: section 2A's full value is reserved for

    an established haploinsufficient gene. Awarding it by default would call
    every whole-gene deletion in the genome Pathogenic.
    """
    gene = synthetic_gene_config(minus_gene, lof_mechanism=mechanism)
    config = CNVScoringConfig()
    whole = pick_cnv(
        minus_gene.transcript, first="1", last="6", consequence=Consequence.EXON_DELETION
    )

    scored = score_cnv(whole, gene, config)

    assert (scored.score, scored.section) == (expected_score, "2A")
    assert config.classify(scored.score) is expected_class


@pytest.mark.parametrize(
    ("mechanism", "expected_rationale"),
    [
        ("established", "complete deletion of TOYM, an established haploinsufficient gene"),
        (
            "not_established",
            "complete deletion of TOYM, whose loss-of-function mechanism is not_established",
        ),
        ("unknown", "complete deletion of TOYM, whose loss-of-function mechanism is unknown"),
    ],
)
def test_the_rationale_of_a_whole_gene_deletion_names_the_mechanism_it_relied_on(
    minus_gene, mechanism, expected_rationale
):
    """The rationale is the audit trail for a curator re-deriving the score.

    Two of the three arms share a score, so the sentence is the only thing that
    distinguishes "we know LoF is not the mechanism" from "nobody has looked".
    """
    gene = synthetic_gene_config(minus_gene, lof_mechanism=mechanism)
    whole = pick_cnv(
        minus_gene.transcript, first="1", last="6", consequence=Consequence.EXON_DELETION
    )

    scored = score_cnv(whole, gene, CNVScoringConfig())

    assert scored.rationale == expected_rationale


# ---------------------------------------------------------------------------
# Deletion vs duplication, whole-gene vs intragenic: the four corners.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("first", "last", "consequence", "expected"),
    [
        # Whole gene, both directions: a deletion removes the gene, a
        # duplication leaves both copies intact and is not LoF evidence.
        ("1", "6", Consequence.EXON_DELETION, CNVScore(100, "2A", "")),
        ("1", "6", Consequence.EXON_DUPLICATION, CNVScore(0, "3A", "")),
        # Intragenic, frame-disrupting (exon 2 alone is 40 coding nt).
        ("2", "2", Consequence.EXON_DELETION, CNVScore(90, "2E", "")),
        ("2", "2", Consequence.EXON_DUPLICATION, CNVScore(90, "2E", "")),
        # Intragenic, frame-preserving (exon 3 alone is 30 coding nt).
        ("3", "3", Consequence.EXON_DELETION, CNVScore(15, "2E", "")),
        ("3", "3", Consequence.EXON_DUPLICATION, CNVScore(0, "3A", "")),
    ],
)
def test_the_four_corners_of_direction_and_extent_score_on_their_own_terms(
    minus_gene, minus_config, first, last, consequence, expected
):
    """Deletion/duplication crossed with whole-gene/intragenic, all six cells.

    Read as a table this is where a copy-and-paste between the loss branch and
    the gain branch becomes visible: an intragenic duplication that preserves
    the frame must score 0 while the equivalent deletion scores 15, and a
    whole-gene duplication must score 0 while the equivalent deletion scores
    100. Any of those four numbers leaking into the wrong branch is a silent
    reclassification.
    """
    variant = pick_cnv(minus_gene.transcript, first=first, last=last, consequence=consequence)

    scored = score_cnv(variant, minus_config, CNVScoringConfig())

    assert (scored.score, scored.section) == (expected.score, expected.section)


@pytest.mark.parametrize(
    ("first", "last", "consequence", "expected"),
    [
        (
            "2",
            "2",
            Consequence.EXON_DELETION,
            "intragenic deletion of 1 exon(s) disrupting the reading frame",
        ),
        (
            "3",
            "3",
            Consequence.EXON_DELETION,
            "intragenic deletion of 1 exon(s) preserving the reading frame",
        ),
        (
            "2",
            "3",
            Consequence.EXON_DUPLICATION,
            "intragenic duplication of 2 exon(s) disrupting the reading frame",
        ),
        (
            "2",
            "5",
            Consequence.EXON_DUPLICATION,
            "intragenic duplication of 3 exon(s) preserving the reading frame",
        ),
    ],
)
def test_an_intragenic_rationale_reports_the_exon_count_and_the_frame_outcome(
    minus_gene, minus_config, first, last, consequence, expected
):
    """The count is the number of exons, not the span of exon labels.

    The fixture's labels skip ``4``, so a rationale built from label arithmetic
    instead of the enumerated count would report ``4`` exons for the 2..5 span.
    """
    variant = pick_cnv(minus_gene.transcript, first=first, last=last, consequence=consequence)

    scored = score_cnv(variant, minus_config, CNVScoringConfig())

    assert scored.rationale == expected


# ---------------------------------------------------------------------------
# The emitted row: what a query over the CNV file actually sees.
# ---------------------------------------------------------------------------


def test_a_frame_disrupting_intragenic_deletion_reports_no_remaining_gap_to_lp(
    minus_gene, minus_config
):
    """Score 90 is exactly the LP threshold, so the gap to LP is closed (None)

    while the distance to LB stays reportable. The boundary matters: an
    off-by-one here would advertise a gap of 0, which a downstream report would
    render as "one hundredth away" rather than "already there".
    """
    config = CNVScoringConfig()
    variant = pick_cnv(
        minus_gene.transcript, first="2", last="2", consequence=Consequence.EXON_DELETION
    )

    row = cnv_row(variant, minus_gene.transcript, minus_config, config, computed_at=FIXED_TIME)

    assert (row.points_current, row.class_current) == (90, ACMGClass.LIKELY_PATHOGENIC)
    assert (row.gap_to_LP, row.gap_to_LB) == (None, 180)
    assert row.blocking_reason is BlockingReason.RESOLVED_NOT_BLOCKED


def test_an_unresolved_whole_gene_duplication_reports_the_distance_to_both_verdicts(
    minus_gene, minus_config
):
    """A zero-scoring duplication is the archetypal CNV gap: equidistant from

    both thresholds, with nothing intrinsic left to add, so the row must say
    what kind of evidence is missing rather than leaving the reason blank.

    The direction assertion pins behaviour this test does *not* endorse. My
    expectation was that a criterion worth zero points carries no direction, or
    is not emitted as applied at all; the builder splits on ``score >= 0`` and
    so files it as pathogenic. That is characterised here, not corrected, so
    that a count of pathogenic-direction CNV criteria is at least known to
    include every whole-gene duplication in the gene.
    """
    config = CNVScoringConfig()
    variant = pick_cnv(
        minus_gene.transcript, first="1", last="6", consequence=Consequence.EXON_DUPLICATION
    )

    row = cnv_row(variant, minus_gene.transcript, minus_config, config, computed_at=FIXED_TIME)

    assert (row.points_current, row.class_current) == (0, ACMGClass.UNCERTAIN)
    assert (row.gap_to_LP, row.gap_to_LB) == (90, 90)
    assert row.blocking_reason is BlockingReason.MISSING_CASE_CONTROL
    assert row.minimum_sufficient_sets == ()
    assert (row.criteria_applied[0].points, row.criteria_applied[0].direction) == (
        0,
        Direction.PATHOGENIC,
    )


def test_a_curated_benign_copy_number_value_is_emitted_on_the_benign_side(minus_gene, minus_config):
    """A gene with a documented benign whole-gene duplication scores -1.00.

    Nothing in the default table is negative, so without a configured negative
    value the benign half of the row builder -- the criterion's direction and
    the two gap fields -- is never exercised at all.
    """
    config = CNVScoringConfig(gain_full_gene=-100)
    variant = pick_cnv(
        minus_gene.transcript, first="1", last="6", consequence=Consequence.EXON_DUPLICATION
    )

    row = cnv_row(variant, minus_gene.transcript, minus_config, config, computed_at=FIXED_TIME)

    assert (row.points_current, row.class_current) == (-100, ACMGClass.BENIGN)
    assert (row.gap_to_LP, row.gap_to_LB) == (190, None)
    assert row.criteria_applied[0].direction is Direction.BENIGN


def test_the_applied_criterion_carries_the_framework_section_and_its_own_source(
    minus_gene, minus_config
):
    """One criterion per CNV, coded by framework section and sourced to the

    framework version -- not to the VCEP spec. A row whose single criterion
    claimed a VCEP source would be indistinguishable from a point-scored row
    once the two files are concatenated.
    """
    config = CNVScoringConfig()
    variant = pick_cnv(
        minus_gene.transcript, first="2", last="2", consequence=Consequence.EXON_DELETION
    )

    row = cnv_row(variant, minus_gene.transcript, minus_config, config, computed_at=FIXED_TIME)

    assert len(row.criteria_applied) == 1
    criterion = row.criteria_applied[0]
    assert (criterion.code, criterion.points, criterion.source) == (
        "CNV_2E",
        90,
        "clingen_cnv@1.0",
    )
    assert (criterion.direction, criterion.strength, criterion.evidence_class) == (
        Direction.PATHOGENIC,
        Strength.SUPPORTING,
        EvidenceClass.INTRINSIC,
    )
    assert row.criteria_evaluated_not_applied == ()


@pytest.mark.parametrize(
    ("consequence", "expected_id"),
    [
        (Consequence.EXON_DELETION, "cnv_1_6_exon_deletion"),
        (Consequence.EXON_DUPLICATION, "cnv_1_6_exon_duplication"),
    ],
)
def test_a_deletion_and_a_duplication_of_one_span_are_not_the_same_class(
    minus_gene, minus_config, consequence, expected_id
):
    """The equivalence class keys on the span *and* the direction of change.

    Dropping the direction would collapse a 100-point deletion and a 0-point
    duplication of the same exons into one class, which is the CNV form of the
    scale-mixing this module exists to prevent.
    """
    variant = pick_cnv(minus_gene.transcript, first="1", last="6", consequence=consequence)

    row = cnv_row(
        variant, minus_gene.transcript, minus_config, CNVScoringConfig(), computed_at=FIXED_TIME
    )

    assert row.equivalence_class_id == expected_id


def test_a_cnv_row_has_no_protein_change_and_no_mutational_distance(minus_gene, minus_config):
    """A CNV row stands for a class of alleles, not one substitution.

    ``mutational_distance`` is a nucleotide count on the sequence scale; leaving
    it at the SNV default of 1 would make copy-number rows answer queries about
    single-nucleotide reachability.
    """
    variant = pick_cnv(
        minus_gene.transcript, first="2", last="5", consequence=Consequence.EXON_DELETION
    )

    row = cnv_row(
        variant, minus_gene.transcript, minus_config, CNVScoringConfig(), computed_at=FIXED_TIME
    )

    assert (row.hgvs_p, row.mutational_distance) == (None, 0)
    assert (row.variant_kind, row.consequence) == (VariantKind.CNV, Consequence.EXON_DELETION)


def test_the_provenance_of_a_cnv_row_records_the_snapshots_it_was_given(minus_gene, minus_config):
    """Provenance is passed through verbatim so the CNV file can be joined to

    the sequence file by snapshot. Defaulting these to None inside the builder
    would produce rows that look reproducible but name no input state.
    """
    config = CNVScoringConfig()
    variant = pick_cnv(
        minus_gene.transcript, first="3", last="3", consequence=Consequence.EXON_DELETION
    )

    row = cnv_row(
        variant,
        minus_gene.transcript,
        minus_config,
        config,
        computed_at=FIXED_TIME,
        clinvar_snapshot=date(2024, 3, 1),
        gnomad_version="4.1",
    )

    assert (row.clinvar_snapshot, row.gnomad_version) == (date(2024, 3, 1), "4.1")
    assert row.computed_at == FIXED_TIME
    assert row.source_versions == {"cnv_framework": "clingen_cnv@1.0"}


def test_a_cnv_score_is_already_at_its_ceiling_because_no_intrinsic_evidence_remains(
    minus_gene, minus_config
):
    """Current and ceiling coincide for every CNV: the framework consumes all

    the intrinsic evidence there is in one pass. If the two ever diverged, the
    gap columns would be describing evidence this module cannot produce.
    """
    config = CNVScoringConfig()
    variant = pick_cnv(
        minus_gene.transcript, first="3", last="3", consequence=Consequence.EXON_DELETION
    )

    row = cnv_row(variant, minus_gene.transcript, minus_config, config, computed_at=FIXED_TIME)

    assert (row.points_current, row.points_ceiling_intrinsic) == (15, 15)
    assert (row.class_current, row.class_ceiling_intrinsic) == (
        ACMGClass.UNCERTAIN,
        ACMGClass.UNCERTAIN,
    )
