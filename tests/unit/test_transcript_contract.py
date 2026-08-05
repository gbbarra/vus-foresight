"""What the transcript model refuses to be built from.

Module-local, so it lives here rather than in a layer: layer 2 asks whether a
coordinate is *right*, this asks whether an impossible transcript can exist at
all. Both matter, and the second is cheaper to get wrong.

The refusals are the point. A transcript whose exons overlap, or whose CDS is
not a whole number of codons, or whose sequence is a different length from its
exons, would still answer every coordinate question asked of it -- with
answers that are internally consistent and biologically wrong. That is the
failure mode this project names as its worst, and the only place to stop it is
before the object exists.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from vus_foresight.genome.transcript import Exon, Transcript


def _fields(gene, **overrides):
    transcript = gene.transcript
    base = {
        "transcript_id": transcript.transcript_id,
        "gene": transcript.gene,
        "chrom": transcript.chrom,
        "strand": transcript.strand,
        "exons": transcript.exons,
        "cds_start_tx": transcript.cds_start_tx,
        "cds_end_tx": transcript.cds_end_tx,
        "sequence": transcript.sequence,
    }
    base.update(overrides)
    return base


# --------------------------------------------------------------------------
# The exon itself.
# --------------------------------------------------------------------------


def test_an_exon_cannot_end_before_it_starts():
    with pytest.raises(ValidationError, match="end 100 < start 200"):
        Exon(label="1", start=200, end=100)


def test_a_single_base_exon_is_legal():
    """The boundary of the previous rule: equal start and end is length one."""
    assert Exon(label="1", start=200, end=200).length == 1


# --------------------------------------------------------------------------
# The exon list.
# --------------------------------------------------------------------------


def test_a_transcript_needs_at_least_one_exon(minus_gene):
    with pytest.raises(ValidationError, match="at least one exon"):
        Transcript(**_fields(minus_gene, exons=(), sequence=""))


def test_the_sequence_must_be_exactly_as_long_as_the_exons(minus_gene):
    """A mismatch here shifts every coordinate downstream of it.

    It is also the exact symptom of a GTF and a FASTA from different releases,
    which is why the MANE extractor checks the same equality.
    """
    truncated = minus_gene.transcript.sequence[:-1]
    with pytest.raises(ValidationError, match="does not match summed exon length"):
        Transcript(**_fields(minus_gene, sequence=truncated))


def test_duplicate_exon_labels_are_refused(minus_gene):
    """Labels are clinical names and are used to report a position back."""
    exons = minus_gene.transcript.exons
    clashing = (*exons[:-1], exons[-1].model_copy(update={"label": exons[0].label}))
    with pytest.raises(ValidationError, match="duplicate exon labels"):
        Transcript(**_fields(minus_gene, exons=clashing))


def test_minus_strand_exons_must_descend(minus_gene):
    """Sorting a minus-strand gene by position reverses it, silently.

    The resulting transcript is internally consistent: every coordinate
    resolves, every HGVS string is well formed, and every one of them points at
    the wrong base. Nothing downstream can detect it.
    """
    ascending = tuple(sorted(minus_gene.transcript.exons, key=lambda e: e.start))
    with pytest.raises(ValidationError, match="out of transcript order on the minus strand"):
        Transcript(**_fields(minus_gene, exons=ascending))


def test_plus_strand_exons_must_ascend(plus_gene):
    descending = tuple(sorted(plus_gene.transcript.exons, key=lambda e: e.start, reverse=True))
    with pytest.raises(ValidationError, match="out of transcript order on the plus strand"):
        Transcript(**_fields(plus_gene, exons=descending))


def test_overlapping_exons_are_refused(plus_gene):
    """Adjacency is legal, overlap is not -- the boundary is one base.

    The second exon is slid back rather than stretched, so the summed length
    stays put: otherwise the sequence-length check fires first and this would
    pass for the wrong reason.
    """
    first, second, *rest = plus_gene.transcript.exons
    slid = second.model_copy(update={"start": first.end, "end": first.end + second.length - 1})
    with pytest.raises(ValidationError, match="overlap"):
        Transcript(**_fields(plus_gene, exons=(first, slid, *rest)))


# --------------------------------------------------------------------------
# The coding sequence.
# --------------------------------------------------------------------------


def test_the_cds_must_end_after_it_starts(minus_gene):
    start = minus_gene.transcript.cds_start_tx
    with pytest.raises(ValidationError, match="CDS end is not after CDS start"):
        Transcript(**_fields(minus_gene, cds_end_tx=start))


def test_the_cds_cannot_run_past_the_transcript(minus_gene):
    transcript = minus_gene.transcript
    over = len(transcript.sequence) + 3
    with pytest.raises(ValidationError, match="CDS end runs past the transcript"):
        Transcript(**_fields(minus_gene, cds_end_tx=over))


def test_the_cds_must_be_a_whole_number_of_codons(minus_gene):
    """Off by one here reframes the entire protein."""
    transcript = minus_gene.transcript
    with pytest.raises(ValidationError, match="is not a multiple of 3"):
        Transcript(**_fields(minus_gene, cds_end_tx=transcript.cds_end_tx - 1))


@pytest.mark.parametrize("bad", [0, -1])
def test_a_cds_bound_must_be_positive(minus_gene, bad):
    """1-based coordinates: there is no position zero to be off by one into."""
    with pytest.raises(ValidationError):
        Transcript(**_fields(minus_gene, cds_start_tx=bad))


# --------------------------------------------------------------------------
# Asking a built transcript for a position it does not have.
# --------------------------------------------------------------------------


def test_a_transcript_position_outside_the_transcript_is_refused(minus_gene):
    transcript = minus_gene.transcript
    for tx_pos in (0, len(transcript.sequence) + 1):
        with pytest.raises(ValueError, match=r"outside 1\.\.|out of range"):
            transcript.genomic_at(tx_pos)


def test_the_first_and_last_transcript_positions_do_resolve(minus_gene):
    """The boundary of the previous test, from the legal side.

    On the minus strand transcript position 1 is the exon's *end*, which is the
    one place where strand becomes arithmetic.
    """
    transcript = minus_gene.transcript
    assert transcript.genomic_at(1) == transcript.exons[0].end
    assert transcript.genomic_at(transcript.length) == transcript.exons[-1].start


def test_c_zero_does_not_exist(minus_gene):
    """HGVS numbering skips zero; accepting it would shift the frame by one."""
    from vus_foresight.genome.transcript import CPosition

    with pytest.raises(ValueError, match=r"c\.0 does not exist"):
        minus_gene.transcript.tx_at_c(CPosition(base=0))


@pytest.mark.parametrize("offset", [-1, 1])
def test_a_cds_position_outside_the_coding_sequence_is_refused(minus_gene, offset):
    transcript = minus_gene.transcript
    bad = transcript.cds_length + 1 if offset > 0 else 0
    with pytest.raises(ValueError, match=r"outside 1\.\.|does not exist"):
        transcript.cds_position_to_tx(bad)


def test_the_first_and_last_cds_positions_do_resolve(minus_gene):
    transcript = minus_gene.transcript
    assert transcript.cds_position_to_tx(1) == transcript.cds_start_tx
    assert transcript.cds_position_to_tx(transcript.cds_length) == transcript.cds_end_tx


@pytest.mark.parametrize("codon", [0, "past the end"])
def test_a_codon_index_outside_the_protein_is_refused(minus_gene, codon):
    transcript = minus_gene.transcript
    index = transcript.n_codons + 1 if codon == "past the end" else codon
    with pytest.raises(ValueError, match="outside 1"):
        transcript.codon_sequence(index)


def test_the_first_and_last_codons_do_resolve(minus_gene):
    """Codon 1 is the start codon and the last one is the terminator."""
    from vus_foresight.genome.sequence import translate_codon

    transcript = minus_gene.transcript
    assert transcript.codon_sequence(1) == "ATG"
    assert translate_codon(transcript.codon_sequence(transcript.n_codons)) == "*"
