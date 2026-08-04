"""Layer 2: the mandatory boundary fixtures of spec section 11.

First and last base of the CDS, the initiation codon, the terminator and
stop-loss, the four canonical splice positions of every junction, the NMD
boundary, and an exon-label discontinuity.
"""

from __future__ import annotations

import pytest

from vus_foresight.annotate.annotator import CodingEdit, annotate_coding_edits
from vus_foresight.annotate.indel import annotate_inframe_deletion
from vus_foresight.enumeration.snv import enumerate_intronic_snvs
from vus_foresight.variant import Consequence

STRANDS = ("+", "-")


@pytest.mark.parametrize("strand", STRANDS)
def test_first_cds_base_is_the_initiation_codon(both_strands, strand):
    transcript = both_strands[strand].transcript
    assert transcript.cds[:3] == "ATG"
    ref = transcript.cds[0]
    variant = annotate_coding_edits(
        transcript, (CodingEdit(cds_position=1, ref=ref, alt="C" if ref != "C" else "G"),)
    )
    assert variant.hgvs_c.startswith("c.1")
    assert variant.consequence is Consequence.START_LOST
    assert variant.hgvs_p == "p.Met1?"


@pytest.mark.parametrize("strand", STRANDS)
def test_last_cds_base_belongs_to_the_terminator(both_strands, strand):
    transcript = both_strands[strand].transcript
    assert transcript.cds[-3:] in {"TAA", "TAG", "TGA"}
    last = transcript.cds_length
    ref = transcript.cds[-1]
    variant = annotate_coding_edits(
        transcript, (CodingEdit(cds_position=last, ref=ref, alt="C" if ref != "C" else "G"),)
    )
    assert variant.hgvs_c.startswith(f"c.{last}")
    assert variant.codon_index == transcript.n_codons
    assert variant.consequence in (Consequence.STOP_LOST, Consequence.STOP_RETAINED)


@pytest.mark.parametrize("strand", STRANDS)
def test_stop_loss_is_named_as_an_extension(both_strands, strand):
    transcript = both_strands[strand].transcript
    first_of_stop = transcript.cds_length - 2
    ref = transcript.cds[first_of_stop - 1]
    alt = next(b for b in "ACGT" if b != ref and b != "T")
    variant = annotate_coding_edits(
        transcript, (CodingEdit(cds_position=first_of_stop, ref=ref, alt=alt),)
    )
    assert variant.consequence is Consequence.STOP_LOST
    assert variant.hgvs_p.startswith(f"p.Ter{transcript.n_codons}")
    assert variant.hgvs_p.endswith("extTer?")


@pytest.mark.parametrize("strand", STRANDS)
def test_the_four_canonical_splice_positions_of_every_junction(both_strands, strand):
    """+1, +2 are donor; -2, -1 are acceptor; +3..+8 and -8..-3 are splice region."""
    gene = both_strands[strand]
    transcript = gene.transcript
    by_offset: dict[int, set[Consequence]] = {}
    for variant in enumerate_intronic_snvs(transcript, gene.flanks, flank_bp=50):
        offset = int(variant.attributes["intron_offset"])
        by_offset.setdefault(offset, set()).add(variant.consequence)

    for offset in (1, 2):
        assert by_offset[offset] == {Consequence.SPLICE_DONOR}
    for offset in (-1, -2):
        assert by_offset[offset] == {Consequence.SPLICE_ACCEPTOR}
    for offset in (3, 8, -3, -8):
        assert by_offset[offset] == {Consequence.SPLICE_REGION}
    for offset in (9, -9, 20, -20):
        assert by_offset[offset] == {Consequence.INTRONIC}


@pytest.mark.parametrize("strand", STRANDS)
def test_canonical_donor_and_acceptor_dinucleotides_are_gt_and_ag(both_strands, strand):
    """The fixture's junctions are GT..AG on the *coding* strand.

    That is the assertion a reverse-complement bug breaks on the minus-strand
    fixture and not on the plus-strand one.
    """
    gene = both_strands[strand]
    transcript = gene.transcript
    for intron in transcript.introns:
        donor = "".join(gene.flanks.base_at(transcript, intron, n) for n in (1, 2))
        acceptor = "".join(
            gene.flanks.base_at(transcript, intron, n) for n in (intron.length - 1, intron.length)
        )
        assert donor == "GT", f"intron {intron.upstream_index} donor is {donor}"
        assert acceptor == "AG", f"intron {intron.upstream_index} acceptor is {acceptor}"


@pytest.mark.parametrize("strand", STRANDS)
def test_exonic_splice_region_covers_three_bases_at_each_abutting_edge(both_strands, strand):
    transcript = both_strands[strand].transcript
    region = transcript.exonic_splice_region_tx
    bounds = transcript._exon_tx_bounds
    # The 5' end of the first exon and the 3' end of the last abut no intron.
    assert bounds[0][0] not in region
    assert bounds[-1][1] not in region
    assert bounds[0][1] in region and bounds[0][1] - 2 in region
    assert bounds[-1][0] in region and bounds[-1][0] + 2 in region


@pytest.mark.parametrize("strand", STRANDS)
def test_nmd_boundary_sits_50_nt_into_the_penultimate_exon(both_strands, strand):
    transcript = both_strands[strand].transcript
    penultimate_end_tx = transcript._exon_tx_bounds[-2][1]
    expected_cds = penultimate_end_tx - 50 + 1 - transcript.cds_start_tx + 1
    assert transcript.nmd_escape_cds_start == max(1, expected_cds)

    boundary_codon = (transcript.nmd_escape_cds_start - 1) // 3 + 1
    assert transcript.ptc_escapes_nmd(transcript.protein_length)
    assert not transcript.ptc_escapes_nmd(1)
    assert transcript.ptc_escapes_nmd(boundary_codon + 1)


def test_exon_labels_may_skip_a_number(plus_gene):
    """BRCA1 has no exon 4. Labels are names, and nothing may do arithmetic on them."""
    labels = [e.label for e in plus_gene.transcript.exons]
    assert labels == ["1", "2", "3", "5", "6"]
    assert "4" not in labels
    # The transcript-order index and the clinical label disagree, on purpose.
    assert plus_gene.transcript.exons[3].label == "5"


@pytest.mark.parametrize("strand", STRANDS)
def test_inframe_deletion_is_shifted_to_the_most_3_prime_position(both_strands, strand):
    """HGVS requires the 3'-most representation of a deletion in a repeat."""
    transcript = both_strands[strand].transcript
    cds = transcript.cds
    # Any two positions three apart with the same base give a shiftable deletion.
    shiftable = [p for p in range(4, transcript.cds_length - 9) if cds[p - 1] == cds[p + 2]]
    assert shiftable, "the fixture should contain a shiftable 3 nt deletion"
    start = shiftable[0]
    variant = annotate_inframe_deletion(transcript, start, 3)
    reported_start = int(variant.hgvs_c[2:].split("_")[0])
    assert reported_start > start
