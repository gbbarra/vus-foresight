"""Layer 2: coordinates and HGVS.

The riskiest layer in the system (spec section 11): a bug here produces HGVS
that looks right and is wrong. The oracles used are the ones that do not share
an implementation with the code under test -- round trips, independent
translation, and the strand symmetry between a plus- and a minus-strand gene.
"""

from __future__ import annotations

import pytest

from vus_foresight.annotate.annotator import CodingEdit, annotate_coding_edits
from vus_foresight.genome.sequence import reverse_complement, translate, translate_codon
from vus_foresight.genome.transcript import CPosition, format_c_position

STRANDS = ("+", "-")


@pytest.mark.parametrize("strand", STRANDS)
def test_genomic_to_c_to_genomic_is_the_identity(both_strands, strand):
    transcript = both_strands[strand].transcript
    for tx in range(1, transcript.length + 1):
        genomic = transcript.genomic_at(tx)
        assert transcript.tx_at_genomic(genomic) == tx
        c = transcript.c_at_tx(tx)
        assert transcript.tx_at_c(c) == tx


@pytest.mark.parametrize("strand", STRANDS)
def test_c_coordinates_increase_monotonically_within_each_exon(both_strands, strand):
    """Inside an exon, ascending ``c.`` must track ascending or descending
    genomic position without ever reversing mid-exon."""
    transcript = both_strands[strand].transcript
    step = transcript.step
    for index, (tx_start, tx_end) in enumerate(transcript._exon_tx_bounds):
        previous = transcript.genomic_at(tx_start)
        for tx in range(tx_start + 1, tx_end + 1):
            current = transcript.genomic_at(tx)
            assert current - previous == step, (
                f"exon {transcript.exons[index].label} is not genomically contiguous"
            )
            previous = current


@pytest.mark.parametrize("strand", STRANDS)
def test_translation_of_the_mutant_cds_matches_the_predicted_protein(both_strands, strand):
    """Independent path: apply the edit to the CDS, translate, compare.

    This never calls the annotator's consequence logic -- it rebuilds the coding
    sequence from scratch and reads the amino acid off a fresh translation.
    """
    transcript = both_strands[strand].transcript
    for cds_position in range(1, transcript.cds_length + 1, 7):
        ref = transcript.cds[cds_position - 1]
        for alt in "ACGT":
            if alt == ref:
                continue
            variant = annotate_coding_edits(
                transcript, (CodingEdit(cds_position=cds_position, ref=ref, alt=alt),)
            )
            mutant = transcript.cds[: cds_position - 1] + alt + transcript.cds[cds_position:]
            codon_index = (cds_position - 1) // 3 + 1
            expected_aa = translate(mutant)[codon_index - 1]
            assert variant.alt_aa == expected_aa
            assert variant.ref_aa == translate(transcript.cds)[codon_index - 1]


@pytest.mark.parametrize("strand", STRANDS)
def test_genomic_representation_uses_the_plus_strand(both_strands, strand):
    transcript = both_strands[strand].transcript
    gene = both_strands[strand]
    for cds_position in range(1, transcript.cds_length + 1, 11):
        ref = transcript.cds[cds_position - 1]
        alt = "A" if ref != "A" else "C"
        variant = annotate_coding_edits(
            transcript, (CodingEdit(cds_position=cds_position, ref=ref, alt=alt),)
        )
        chrom, position, g_ref, g_alt = variant.grch38_pos.split("-")
        assert chrom == transcript.chrom
        assert gene.locus[int(position)] == g_ref
        if strand == "+":
            assert (g_ref, g_alt) == (ref, alt)
        else:
            assert (g_ref, g_alt) == (reverse_complement(ref), reverse_complement(alt))


def test_strand_symmetry_of_the_transcript_level_annotation(plus_gene, minus_gene):
    """The same assertions on both genes, demanding identical results.

    BRCA1 is on the minus strand and BRCA2 on the plus, which makes this a free
    internal control: a reverse-complement bug shows up in one and not the
    other. The two fixtures share a transcript sequence, so every ``c.`` and
    ``p.`` description must match exactly while the genomic coordinates differ.
    """
    assert plus_gene.transcript.sequence == minus_gene.transcript.sequence

    plus = [
        annotate_coding_edits(
            plus_gene.transcript,
            (CodingEdit(cds_position=p, ref=plus_gene.transcript.cds[p - 1], alt=a),),
        )
        for p in range(1, plus_gene.transcript.cds_length + 1, 5)
        for a in "ACGT"
        if a != plus_gene.transcript.cds[p - 1]
    ]
    minus = [
        annotate_coding_edits(
            minus_gene.transcript,
            (CodingEdit(cds_position=p, ref=minus_gene.transcript.cds[p - 1], alt=a),),
        )
        for p in range(1, minus_gene.transcript.cds_length + 1, 5)
        for a in "ACGT"
        if a != minus_gene.transcript.cds[p - 1]
    ]

    assert [v.hgvs_c for v in plus] == [v.hgvs_c for v in minus]
    assert [v.hgvs_p for v in plus] == [v.hgvs_p for v in minus]
    assert [v.consequence for v in plus] == [v.consequence for v in minus]
    assert [v.grch38_pos for v in plus] != [v.grch38_pos for v in minus]


@pytest.mark.parametrize("strand", STRANDS)
def test_mnv_translates_the_whole_codon_not_the_component_substitutions(both_strands, strand):
    """The failure mode spec section 11 calls the likeliest in the whole layer.

    A codon where two separate substitutions give amino acids X and Y, but the
    two together give a third amino acid Z. Composing per-position effects gives
    X or Y; translating the assembled codon gives Z.
    """
    transcript = both_strands[strand].transcript
    checked = 0
    for codon_index in range(2, transcript.n_codons):
        ref_codon = transcript.codon_sequence(codon_index)
        first_cds, _, _ = transcript.codon_bounds_cds(codon_index)
        for i in range(3):
            for j in range(i + 1, 3):
                for alt_i in "ACGT":
                    if alt_i == ref_codon[i]:
                        continue
                    for alt_j in "ACGT":
                        if alt_j == ref_codon[j]:
                            continue
                        single_i = list(ref_codon)
                        single_i[i] = alt_i
                        single_j = list(ref_codon)
                        single_j[j] = alt_j
                        both = list(ref_codon)
                        both[i], both[j] = alt_i, alt_j

                        aa_i = translate_codon("".join(single_i))
                        aa_j = translate_codon("".join(single_j))
                        aa_both = translate_codon("".join(both))
                        if aa_both in (aa_i, aa_j):
                            continue

                        variant = annotate_coding_edits(
                            transcript,
                            (
                                CodingEdit(
                                    cds_position=first_cds + i,
                                    ref=ref_codon[i],
                                    alt=alt_i,
                                ),
                                CodingEdit(
                                    cds_position=first_cds + j,
                                    ref=ref_codon[j],
                                    alt=alt_j,
                                ),
                            ),
                        )
                        assert variant.alt_aa == aa_both, (
                            f"codon {codon_index} {ref_codon}: expected {aa_both}, "
                            f"composing per position would give {aa_i}/{aa_j}"
                        )
                        checked += 1
    assert checked > 0, "the fixture contains no third-amino-acid codon to test"


@pytest.mark.parametrize("strand", STRANDS)
def test_mnv_hgvs_uses_a_delins_over_the_minimal_span(both_strands, strand):
    transcript = both_strands[strand].transcript
    first_cds, _, _ = transcript.codon_bounds_cds(5)
    ref_codon = transcript.codon_sequence(5)
    alt_first = "A" if ref_codon[0] != "A" else "C"
    alt_last = "A" if ref_codon[2] != "A" else "C"

    # Positions 1 and 3 of the codon: the span cannot be trimmed.
    wide = annotate_coding_edits(
        transcript,
        (
            CodingEdit(cds_position=first_cds, ref=ref_codon[0], alt=alt_first),
            CodingEdit(cds_position=first_cds + 2, ref=ref_codon[2], alt=alt_last),
        ),
    )
    assert wide.hgvs_c == f"c.{first_cds}_{first_cds + 2}delins{alt_first}{ref_codon[1]}{alt_last}"

    # Positions 1 and 2: the span stops before the unchanged third base.
    alt_middle = "A" if ref_codon[1] != "A" else "C"
    narrow = annotate_coding_edits(
        transcript,
        (
            CodingEdit(cds_position=first_cds, ref=ref_codon[0], alt=alt_first),
            CodingEdit(cds_position=first_cds + 1, ref=ref_codon[1], alt=alt_middle),
        ),
    )
    assert narrow.hgvs_c == f"c.{first_cds}_{first_cds + 1}delins{alt_first}{alt_middle}"


def test_c_position_formatting():
    assert format_c_position(CPosition(base=100)) == "100"
    assert format_c_position(CPosition(base=100, offset=3)) == "100+3"
    assert format_c_position(CPosition(base=101, offset=-2)) == "101-2"
    assert format_c_position(CPosition(base=12, utr3=True)) == "*12"
    assert format_c_position(CPosition(base=12, offset=-1, utr3=True)) == "*12-1"
    assert format_c_position(CPosition(base=-20)) == "-20"


@pytest.mark.parametrize("strand", STRANDS)
def test_intron_numbering_splits_at_the_midpoint(both_strands, strand):
    """HGVS numbers the 5' half from the donor and the 3' half from the acceptor.

    For an odd-length intron the central base belongs to the 5' half. The
    fixture includes a 31 nt intron precisely to pin that down.
    """
    transcript = both_strands[strand].transcript
    odd = [i for i in transcript.introns if i.length % 2 == 1]
    assert odd, "the fixture should contain an odd-length intron"
    intron = odd[0]
    midpoint = (intron.length + 1) // 2

    assert transcript.intron_c_position(intron, 1).offset == 1
    assert transcript.intron_c_position(intron, midpoint).offset == midpoint
    assert transcript.intron_c_position(intron, midpoint + 1).offset == -(intron.length - midpoint)
    assert transcript.intron_c_position(intron, intron.length).offset == -1

    offsets = [transcript.intron_c_position(intron, n).offset for n in range(1, intron.length + 1)]
    assert len(set(offsets)) == intron.length
    assert 0 not in offsets
