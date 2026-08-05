"""Layer 3: the two decision procedures that are not the point engine.

PVS1 is a tree over transcript geometry; copy number is a different framework
entirely. Both are parameterised by configuration, and both have to be honest
about what they cannot determine.
"""

from __future__ import annotations

import pytest

from vus_foresight.acmg import ACMGClass, Strength
from vus_foresight.annotate.annotator import CodingEdit, annotate_coding_edits
from vus_foresight.engine.cnv_scoring import CNVScoringConfig, score_cnv
from vus_foresight.engine.context import MISSING, EvidenceContext
from vus_foresight.engine.pvs1 import compute_pvs1
from vus_foresight.enumeration import enumerate_exon_cnvs
from vus_foresight.variant import Consequence


def _nonsense_at(transcript, codon_index):
    """Find any single substitution in this codon that creates a terminator."""
    first, _, _ = transcript.codon_bounds_cds(codon_index)
    ref_codon = transcript.codon_sequence(codon_index)
    for offset in range(3):
        for alt in "ACGT":
            if alt == ref_codon[offset]:
                continue
            candidate = list(ref_codon)
            candidate[offset] = alt
            if "".join(candidate) in ("TAA", "TAG", "TGA"):
                return annotate_coding_edits(
                    transcript,
                    (
                        CodingEdit(
                            cds_position=first + offset,
                            ref=ref_codon[offset],
                            alt=alt,
                        ),
                    ),
                )
    return None


def test_pvs1_is_very_strong_when_nmd_is_triggered(minus_gene, minus_config, toy_spec):
    transcript = minus_gene.transcript
    early = next(
        v
        for codon in range(2, transcript.protein_length)
        if (v := _nonsense_at(transcript, codon)) is not None
        and not transcript.ptc_escapes_nmd(codon)
    )
    result = compute_pvs1(early, transcript, minus_config, toy_spec.pvs1, EvidenceContext())
    assert result["applicable"] is True
    assert result["nmd_escape"] is False
    assert result["strength"] == Strength.VERY_STRONG.value


def test_pvs1_drops_to_strong_when_lof_mechanism_is_not_established(
    minus_gene, minus_config, toy_spec
):
    from vus_foresight.testing import synthetic_gene_config

    transcript = minus_gene.transcript
    unestablished = synthetic_gene_config(minus_gene, lof_mechanism="not_established")
    early = next(
        v
        for codon in range(2, transcript.protein_length)
        if (v := _nonsense_at(transcript, codon)) is not None
        and not transcript.ptc_escapes_nmd(codon)
    )
    result = compute_pvs1(early, transcript, unestablished, toy_spec.pvs1, EvidenceContext())
    assert result["strength"] == Strength.STRONG.value
    assert "not_established" in result["rationale"]


def test_pvs1_is_undetermined_for_a_splice_site_without_a_prediction(
    minus_gene, minus_config, toy_spec
):
    """The un-ingested SpliceAI table becomes a visible gap, not an assumption."""
    from vus_foresight.enumeration.snv import enumerate_intronic_snvs

    transcript = minus_gene.transcript
    donor = next(
        v
        for v in enumerate_intronic_snvs(transcript, minus_gene.flanks)
        if v.consequence is Consequence.SPLICE_DONOR
    )
    result = compute_pvs1(donor, transcript, minus_config, toy_spec.pvs1, EvidenceContext())
    assert result["applicable"] is False
    assert result["undetermined_reason"] == "splice_prediction_missing"
    assert result["missing_field"] == toy_spec.pvs1.splice_prediction_field


def test_pvs1_fires_for_a_splice_site_once_the_prediction_is_present(
    minus_gene, minus_config, toy_spec
):
    from vus_foresight.enumeration.snv import enumerate_intronic_snvs

    transcript = minus_gene.transcript
    donor = next(
        v
        for v in enumerate_intronic_snvs(transcript, minus_gene.flanks)
        if v.consequence is Consequence.SPLICE_DONOR
    )
    context = EvidenceContext()
    context.merge("splice", {"ds_max": 0.98}, "spliceai@test")
    result = compute_pvs1(donor, transcript, minus_config, toy_spec.pvs1, context)
    assert result["applicable"] is True
    assert result["strength"] == Strength.VERY_STRONG.value


def test_pvs1_uses_the_critical_region_branch_when_nmd_is_escaped(plus_gene, plus_config, toy_spec):
    transcript = plus_gene.transcript
    late = [
        v
        for codon in range(2, transcript.protein_length + 1)
        if (v := _nonsense_at(transcript, codon)) is not None and transcript.ptc_escapes_nmd(codon)
    ]
    assert late, "the fixture should allow an NMD-escaping PTC"
    result = compute_pvs1(late[0], transcript, plus_config, toy_spec.pvs1, EvidenceContext())
    assert result["nmd_escape"] is True
    assert result["strength"] in {Strength.STRONG.value, Strength.MODERATE.value}
    assert 0 < result["fraction_removed"] <= 1


def test_start_loss_is_reported_as_unpredictable(minus_gene, minus_config, toy_spec):
    transcript = minus_gene.transcript
    variant = annotate_coding_edits(
        transcript, (CodingEdit(cds_position=1, ref=transcript.cds[0], alt="C"),)
    )
    result = compute_pvs1(variant, transcript, minus_config, toy_spec.pvs1, EvidenceContext())
    assert result["strength"] == toy_spec.pvs1.start_lost_strength.value
    assert "not predictable" in result["rationale"]


# --------------------------------------------------------------------------
# Copy number: a different framework, deliberately on a different scale.
# --------------------------------------------------------------------------


def test_whole_gene_deletion_of_an_established_hi_gene_is_pathogenic(minus_gene, minus_config):
    config = CNVScoringConfig()
    whole = next(
        v
        for v in enumerate_exon_cnvs(minus_gene.transcript)
        if v.attributes["spans_whole_gene"] == "true" and v.consequence is Consequence.EXON_DELETION
    )
    scored = score_cnv(whole, minus_config, config)
    assert scored.section == "2A"
    assert config.classify(scored.score) is ACMGClass.PATHOGENIC


def test_whole_gene_duplication_is_not_evidence_of_loss_of_function(minus_gene, minus_config):
    config = CNVScoringConfig()
    whole = next(
        v
        for v in enumerate_exon_cnvs(minus_gene.transcript)
        if v.attributes["spans_whole_gene"] == "true"
        and v.consequence is Consequence.EXON_DUPLICATION
    )
    scored = score_cnv(whole, minus_config, config)
    assert scored.score == 0
    assert config.classify(scored.score) is ACMGClass.UNCERTAIN


def test_in_frame_and_frame_disrupting_intragenic_deletions_score_differently(
    minus_gene, minus_config
):
    config = CNVScoringConfig()
    deletions = [
        v
        for v in enumerate_exon_cnvs(minus_gene.transcript)
        if v.consequence is Consequence.EXON_DELETION
        and v.attributes["spans_whole_gene"] == "false"
    ]
    in_frame = [v for v in deletions if v.attributes["in_frame"] == "true"]
    shifted = [v for v in deletions if v.attributes["in_frame"] == "false"]
    assert in_frame and shifted
    assert (
        score_cnv(in_frame[0], minus_config, config).score
        < score_cnv(shifted[0], minus_config, config).score
    )


def test_cnv_rows_declare_their_own_scoring_framework(minus_gene, minus_config):
    from datetime import datetime

    from vus_foresight.engine.cnv_scoring import cnv_row

    config = CNVScoringConfig()
    variant = next(iter(enumerate_exon_cnvs(minus_gene.transcript)))
    row = cnv_row(
        variant,
        minus_gene.transcript,
        minus_config,
        config,
        computed_at=datetime(1970, 1, 1),
    )
    assert row.spec_version == config.version_string
    assert row.spec_version != "toy@0.1.0"
    assert row.source_versions["cnv_framework"] == config.version_string


def test_the_pipeline_never_mixes_cnv_rows_into_the_point_scale(minus_gene, minus_config, toy_spec):
    from datetime import datetime

    from vus_foresight.engine.pipeline import MapRunner, default_registry

    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=default_registry(minus_config),
        computed_at=datetime(1970, 1, 1),
    )
    cnvs = list(enumerate_exon_cnvs(minus_gene.transcript))
    assert cnvs
    assert list(runner.run(cnvs)) == []


# --------------------------------------------------------------------------
# The context's absent/false distinction, which everything above rests on.
# --------------------------------------------------------------------------


def test_absent_is_not_false_and_not_zero():
    context = EvidenceContext()
    context.merge("frequency", {"gnomad": {"faf95_popmax": 0.0}}, "gnomad@v4")
    assert context.get("frequency.gnomad.faf95_popmax") == 0.0
    assert context.get("frequency.gnomad.missing_field") is MISSING
    assert context.get("nowhere.at.all") is MISSING


def test_a_namespace_may_only_be_claimed_once():
    context = EvidenceContext()
    context.merge("frequency", {}, "a@1")
    with pytest.raises(ValueError, match="already populated"):
        context.merge("frequency", {}, "b@1")


# --------------------------------------------------------------------------
# The branches where PVS1 declines to reach a verdict.
#
# Every one of these is the difference between "no evidence" and "the tree
# could not run", and the distinction has to survive into the trace: a
# criterion requiring pvs1.strength reports NOT_EVALUABLE with the reason
# attached rather than quietly failing to fire, which would silently remove
# eight points from every truncating variant in the gene.
# --------------------------------------------------------------------------


def test_a_disabled_pvs1_says_so_instead_of_returning_no_strength(
    minus_gene, minus_config, toy_spec
):
    """A specification that switches PVS1 off must be distinguishable from one
    whose tree ran and declined."""
    transcript = minus_gene.transcript
    variant = next(
        v
        for codon in range(2, transcript.protein_length)
        if (v := _nonsense_at(transcript, codon)) is not None
    )
    disabled = toy_spec.pvs1.model_copy(update={"enabled": False})

    result = compute_pvs1(variant, transcript, minus_config, disabled, EvidenceContext())

    assert result["applicable"] is False
    assert result["undetermined_reason"] == "pvs1_disabled_in_spec"
    assert "strength" not in result


def test_a_non_truncating_consequence_is_declined_with_its_reason(
    minus_gene, minus_config, toy_spec
):
    from vus_foresight.enumeration import enumerate_coding_snvs
    from vus_foresight.variant import Consequence

    transcript = minus_gene.transcript
    missense = next(
        v for v in enumerate_coding_snvs(transcript) if v.consequence is Consequence.MISSENSE
    )

    result = compute_pvs1(missense, transcript, minus_config, toy_spec.pvs1, EvidenceContext())

    assert result["applicable"] is False
    assert result["undetermined_reason"] == "consequence_not_truncating"
    assert "strength" not in result


def test_a_truncating_variant_with_no_locatable_ptc_is_declined_not_assumed(
    minus_gene, minus_config, toy_spec
):
    """Without a codon there is no NMD zone, so there is no tree to walk.

    Guessing here would be the worst available option: the NMD boundary is what
    separates very strong from strong, so a defaulted position would hand out
    the strongest criterion in the system on no information.
    """
    from vus_foresight.variant import Consequence, Variant, VariantKind

    transcript = minus_gene.transcript
    homeless = Variant(
        transcript=transcript.transcript_id,
        gene=transcript.gene,
        hgvs_c="c.100del",
        hgvs_p=None,
        consequence=Consequence.FRAMESHIFT,
        kind=VariantKind.FRAMESHIFT_CLASS,
        codon_index=None,
        ptc_codon=None,
    )

    result = compute_pvs1(homeless, transcript, minus_config, toy_spec.pvs1, EvidenceContext())

    assert result["applicable"] is False
    assert result["undetermined_reason"] == "ptc_position_unknown"
    assert "strength" not in result


def test_a_nonsense_variant_without_an_explicit_ptc_falls_back_to_its_own_codon(
    minus_gene, minus_config, toy_spec
):
    """A nonsense variant *is* its PTC, so the tree can run without being told.

    The enumerator already fills ptc_codon in, which is why this constructs the
    variant by hand: the fallback exists for anything that reaches the tree
    from another path, and an untested fallback in the NMD decision is a
    silent route to the strongest criterion in the system.
    """
    from vus_foresight.variant import Consequence, Variant, VariantKind

    transcript = minus_gene.transcript
    enumerated = next(
        v
        for codon in range(2, transcript.protein_length)
        if (v := _nonsense_at(transcript, codon)) is not None
    )
    assert enumerated.ptc_codon == enumerated.codon_index

    without_ptc = Variant(
        transcript=transcript.transcript_id,
        gene=transcript.gene,
        hgvs_c=enumerated.hgvs_c,
        hgvs_p=enumerated.hgvs_p,
        consequence=Consequence.NONSENSE,
        kind=VariantKind.SNV,
        codon_index=enumerated.codon_index,
        ptc_codon=None,
    )

    result = compute_pvs1(without_ptc, transcript, minus_config, toy_spec.pvs1, EvidenceContext())

    assert result["applicable"] is True
    assert result["ptc_codon"] == enumerated.codon_index
