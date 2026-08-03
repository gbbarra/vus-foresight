"""Layer 1: the enumeration engine.

Failures here produce variants that do not exist, and the oracle is *counting*.
Every assertion in this file is an exact closed form, never an approximation --
"about sixteen thousand" would pass with an off-by-one in the CDS bounds, which
is the error that matters.
"""

from __future__ import annotations

import pytest

from vus_foresight.enumeration import (
    MNV_PER_CODON,
    build_frameshift_classes,
    coding_snv_count,
    enumerate_coding_snvs,
    enumerate_exon_cnvs,
    enumerate_inframe_deletions,
    enumerate_intracodon_mnvs,
    enumerate_intronic_snvs,
    exon_interval_count,
    intronic_snv_positions,
    mnv_count,
)
from vus_foresight.enumeration.mnv import alternative_codons
from vus_foresight.genome.reference import load_gene_config
from vus_foresight.genome.sequence import BASES
from vus_foresight.variant import VariantKind

from ..conftest import CONFIG_DIR

STRANDS = ("+", "-")


@pytest.mark.parametrize("strand", STRANDS)
def test_coding_snv_cardinality_is_three_per_base(both_strands, strand):
    transcript = both_strands[strand].transcript
    variants = list(enumerate_coding_snvs(transcript))
    assert len(variants) == 3 * transcript.cds_length
    assert len(variants) == coding_snv_count(transcript)


@pytest.mark.parametrize("strand", STRANDS)
def test_coding_snv_set_equals_independent_cartesian_product(both_strands, strand):
    """Completeness, computed a second way and compared as sets.

    The independent construction deliberately does not reuse the enumerator's
    ordering, its annotation, or its alternative-base helper.
    """
    transcript = both_strands[strand].transcript
    expected = {
        (position, alt)
        for position, ref in enumerate(transcript.cds, start=1)
        for alt in "ACGT"
        if alt != ref
    }
    produced = {
        (variant.cds_position, variant.hgvs_c.split(">")[-1])
        for variant in enumerate_coding_snvs(transcript)
    }
    assert produced == expected


@pytest.mark.parametrize("strand", STRANDS)
def test_no_duplicates_and_no_no_ops(both_strands, strand):
    transcript = both_strands[strand].transcript
    seen: set[str] = set()
    for variant in enumerate_coding_snvs(transcript):
        assert variant.hgvs_c not in seen, f"duplicate {variant.hgvs_c}"
        seen.add(variant.hgvs_c)
        ref, alt = variant.hgvs_c.split(">")[0][-1], variant.hgvs_c.split(">")[1]
        assert ref != alt, f"no-op variant {variant.hgvs_c}"


@pytest.mark.parametrize("strand", STRANDS)
def test_intronic_window_covers_both_sides_of_every_junction(both_strands, strand):
    transcript = both_strands[strand].transcript
    positions = intronic_snv_positions(transcript, flank_bp=50)
    expected = 0
    for intron in transcript.introns:
        # Both windows, deduplicated where a short intron makes them overlap.
        donor = set(range(1, min(50, intron.length) + 1))
        acceptor = set(range(max(1, intron.length - 50 + 1), intron.length + 1))
        expected += len(donor | acceptor)
    assert len(positions) == expected
    assert len(set(positions)) == len(positions)


@pytest.mark.parametrize("strand", STRANDS)
def test_intronic_snvs_are_three_per_known_position(both_strands, strand):
    gene = both_strands[strand]
    positions = intronic_snv_positions(gene.transcript, flank_bp=50)
    variants = list(enumerate_intronic_snvs(gene.transcript, gene.flanks, flank_bp=50))
    assert len(variants) == 3 * len(positions)
    assert len({v.hgvs_c for v in variants}) == len(variants)


@pytest.mark.parametrize("strand", STRANDS)
def test_intronic_snvs_without_flanks_emit_four_alternatives_and_say_so(
    both_strands, strand
):
    """Positional coverage stays complete; the unknown reference is declared."""
    transcript = both_strands[strand].transcript
    positions = intronic_snv_positions(transcript, flank_bp=50)
    variants = list(enumerate_intronic_snvs(transcript, None, flank_bp=50))
    assert len(variants) == 4 * len(positions)
    assert all(v.attributes.get("reference_base_unknown") == "true" for v in variants)


def test_alternative_codons_are_exactly_the_54_multi_nucleotide_ones():
    for ref in ("GCC", "ATG", "TAA", "AAA"):
        alternatives = alternative_codons(ref)
        assert len(alternatives) == MNV_PER_CODON == 54
        assert all(len(differing) >= 2 for _, differing in alternatives)
        assert len({alt for alt, _ in alternatives}) == 54
        # 63 alternatives in total, 9 of which are single-nucleotide.
        single = [
            "".join(ref[:i] + b + ref[i + 1 :])
            for i in range(3)
            for b in BASES
            if b != ref[i]
        ]
        assert len(single) == 9
        assert not {alt for alt, _ in alternatives} & set(single)


@pytest.mark.parametrize("strand", STRANDS)
def test_mnv_cardinality_and_distance(both_strands, strand):
    transcript = both_strands[strand].transcript
    variants = list(enumerate_intracodon_mnvs(transcript))
    assert len(variants) == mnv_count(transcript) == 54 * transcript.n_codons
    assert all(v.mutational_distance in (2, 3) for v in variants)
    assert not any(v.mutational_distance == 1 for v in variants)
    assert all(v.kind is VariantKind.MNV for v in variants)
    assert len({v.hgvs_c for v in variants}) == len(variants)


@pytest.mark.parametrize("strand", STRANDS)
def test_frameshift_classes_account_for_every_codon_position(both_strands, strand):
    transcript = both_strands[strand].transcript
    enumeration = build_frameshift_classes(transcript)
    assert enumeration.total_codon_positions == transcript.protein_length
    reachable = {c.ptc_codon for c in enumeration.classes}
    assert reachable.isdisjoint(set(enumeration.unreachable))
    assert all(1 <= c.ptc_codon <= transcript.protein_length for c in enumeration.classes)


@pytest.mark.parametrize("strand", STRANDS)
def test_frameshift_class_members_really_produce_that_ptc(both_strands, strand):
    """Verified by translation, not by trusting the metadata that grouped them."""
    from vus_foresight.enumeration import verify_class_membership

    transcript = both_strands[strand].transcript
    enumeration = build_frameshift_classes(transcript)
    assert enumeration.classes, "the fixture should reach at least one PTC class"
    for cls in enumeration.classes:
        assert cls.representatives
        for indel in cls.representatives:
            assert verify_class_membership(transcript, indel, cls.ptc_codon), (
                f"{indel} does not actually produce a PTC at codon {cls.ptc_codon}"
            )


@pytest.mark.parametrize("strand", STRANDS)
def test_exon_cnv_interval_count(both_strands, strand):
    transcript = both_strands[strand].transcript
    n = len(transcript.exons)
    assert exon_interval_count(transcript) == n * (n + 1) // 2
    variants = list(enumerate_exon_cnvs(transcript))
    # One deletion and one duplication per interval.
    assert len(variants) == 2 * n * (n + 1) // 2
    assert len({v.hgvs_c for v in variants}) == len(variants)


@pytest.mark.parametrize("strand", STRANDS)
def test_inframe_deletions_are_distinct_and_in_frame(both_strands, strand):
    transcript = both_strands[strand].transcript
    variants = list(enumerate_inframe_deletions(transcript))
    assert len({v.hgvs_c for v in variants}) == len(variants)
    assert all(v.mutational_distance in (3, 6) for v in variants)
    assert all(v.kind is VariantKind.INFRAME_DELETION for v in variants)


# --------------------------------------------------------------------------
# The declared BRCA cardinalities of spec section 11. These run without any
# reference resource because the counts follow from the declared CDS length --
# which is exactly the number a config typo would get wrong.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("gene", "expected_snvs", "expected_mnvs", "expected_classes", "expected_intervals"),
    [
        ("BRCA1", 16_776, 100_656, 1_863, 276),
        ("BRCA2", 30_771, 184_626, 3_418, 378),
    ],
)
def test_declared_gene_configs_imply_the_specified_cardinalities(
    gene, expected_snvs, expected_mnvs, expected_classes, expected_intervals
):
    config = load_gene_config(CONFIG_DIR / "genes" / f"{gene}.yaml")
    tc = config.transcript
    assert 3 * tc.cds_length == expected_snvs
    assert 54 * (tc.cds_length // 3) == expected_mnvs
    assert tc.protein_length == expected_classes
    assert tc.exon_count * (tc.exon_count + 1) // 2 == expected_intervals


@pytest.mark.requires_reference
def test_brca1_coding_snv_count(brca1):
    _config, transcript = brca1
    assert len(list(enumerate_coding_snvs(transcript))) == 16_776


@pytest.mark.requires_reference
def test_brca2_coding_snv_count(brca2):
    _config, transcript = brca2
    assert len(list(enumerate_coding_snvs(transcript))) == 30_771
