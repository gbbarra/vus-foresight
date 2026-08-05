"""Edge cases in the frameshift / in-frame-deletion machinery.

Oracle: hand-translation of a sequence short enough to read. Every expectation
below is derived by writing the mutant sequence out, splitting it into codons
with the standard genetic code, and counting -- never by running the code and
recording what came back. That is what makes these falsifiable: the way they can
be *wrong* is that the hand translation was done wrong, which a second reader can
check in a minute, not that they were bent to agree with the implementation.

Scope note. Three things in this module are "plausible and wrong" hazards, and
they decide what is pinned here:

* :class:`StopIndex` is a *tuple* behind ``__getitem__``. A query outside the
  table must answer "no terminator", but a tuple answers a negative index by
  wrapping to the 3' end -- which would silently report a real stop codon from
  the wrong end of the transcript. The guard is pinned in both directions.
* :func:`ptc_after_indel` returns a 1-based *codon* number computed from a
  0-based *nucleotide* offset. An off-by-one here is invisible: the number stays
  in range and still looks like a residue.
* The reference-sequence guards (``cds_position >= 1``, positive deletion
  length, deletion inside the CDS) are what stop a coordinate bug upstream from
  being absorbed into a well-formed but wrong HGVS string.

The end-to-end 3'-normalisation of a deletion, and the strand symmetry of the
resulting genomic coordinate, are covered in ``tests/layer2_annotation`` and are
not repeated.
"""

from __future__ import annotations

import pytest

from vus_foresight.annotate.indel import (
    ProteinChange,
    annotate_inframe_deletion,
    build_stop_index,
    first_stop_index,
    normalize_deletion_3prime,
    protein_change,
    ptc_after_indel,
)
from vus_foresight.genome.transcript import Exon, Transcript
from vus_foresight.variant import Consequence

#: ``ATG AAA TAA GGG``. Chosen because each of the three reading frames gives a
#: different answer: frame 0 stops at the ``TAA`` (offset 6), frame 1 stops
#: immediately on the ``TGA`` spanning the first two codons (offset 1), and
#: frame 2 never meets a terminator at all.
THREE_FRAMES = "ATGAAATAAGGG"

#: Exon start used by the hand-built transcripts, so a genomic assertion can be
#: read as ``EXON_START + (transcript position - 1)``.
EXON_START = 1000


def _transcript_from_cds(cds: str) -> Transcript:
    """A single-exon, plus-strand transcript whose whole sequence is the CDS.

    Single-exon on purpose: with no intron there is no splice region and no
    genomic discontinuity, so anything the annotator reports about a deletion is
    attributable to the deletion alone.
    """
    return Transcript(
        transcript_id="NM_000000.1",
        gene="TOYX",
        chrom="chr1",
        strand="+",
        exons=(Exon(label="1", start=EXON_START, end=EXON_START + len(cds) - 1),),
        cds_start_tx=1,
        cds_end_tx=len(cds),
        sequence=cds,
    )


# ----------------------------------------------------------------------
# StopIndex lookup
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "offset",
    [-1, -3, -12, 12, 13, 500],
    ids=["just-before", "negative-codon", "negative-whole-length", "one-past", "past", "far-past"],
)
def test_a_query_outside_the_indexed_sequence_reports_no_terminator(offset):
    """A negative offset must not wrap to the 3' end of the table.

    ``table`` is a tuple, so without the guard ``index[-3]`` would answer with
    the last entry -- a terminator that exists, at a position nowhere near the
    query. That is the worst shape of bug this project has: a real number,
    attached to the wrong coordinate.
    """
    index = build_stop_index(THREE_FRAMES)

    result = index[offset]

    assert result is None


@pytest.mark.parametrize(
    ("offset", "expected"),
    [(0, 6), (3, 6), (1, 1), (4, None), (2, None), (11, None)],
    ids=[
        "frame0-from-the-start",
        "frame0-second-codon",
        "frame1-immediate",
        "frame1-after-the-stop",
        "frame2-never",
        "last-base",
    ],
)
def test_one_table_answers_all_three_reading_frames(offset, expected):
    """``ATG AAA TAA GGG`` read in frame 1 is ``TGA``, a stop at offset 1.

    The recurrence steps by three, so the entry for an offset only ever sees
    offsets congruent to it mod 3. If it ever consulted ``table[j + 1]``, frame 2
    of this sequence would inherit frame 0's terminator and every frameshift PTC
    in a real gene would be off by a fraction of a codon.
    """
    index = build_stop_index(THREE_FRAMES)

    result = index[offset]

    assert result == expected


@pytest.mark.parametrize(
    ("offset", "expected"),
    [(0, 6), (1, 1), (2, None), (-1, None), (99, None)],
    ids=["frame0", "frame1", "frame2", "negative", "past-the-end"],
)
def test_a_one_off_query_agrees_with_the_precomputed_table(offset, expected):
    """The convenience wrapper must not be a second, divergent implementation.

    It exists only so a caller with a single question does not have to build a
    table; the day it stops agreeing with the table, two parts of the pipeline
    disagree about where a gene's terminators are.
    """
    result = first_stop_index(THREE_FRAMES, offset)

    assert result == expected


# ----------------------------------------------------------------------
# ptc_after_indel
# ----------------------------------------------------------------------


@pytest.mark.parametrize("cds_position", [0, -1, -100], ids=["zero", "minus-one", "far-negative"])
def test_a_cds_position_below_one_is_rejected_rather_than_read_as_an_offset(cds_position):
    """``c.0`` does not exist, and a 0-based caller must fail loudly.

    Without the guard, ``cds_position=0`` would index ``sequence[-1]`` and build
    a head codon out of the last base of the 3' UTR -- producing a PTC codon
    number that is in range and meaningless.
    """
    index = build_stop_index(THREE_FRAMES)

    with pytest.raises(ValueError, match=r"cds_position is 1-based and must be >= 1"):
        ptc_after_indel(index, cds_position, 1)


@pytest.mark.parametrize(
    ("cds_position", "deleted", "inserted", "expected"),
    [
        (4, 1, "", 2),
        (4, 0, "G", None),
    ],
    ids=["deletion-creates-a-stop-in-the-new-frame", "insertion-runs-out-of-sequence"],
)
def test_the_ptc_codon_is_counted_from_the_mutant_sequence_not_the_reference(
    cds_position, deleted, inserted, expected
):
    """Hand translation of ``ATGATAGCCC`` under each edit.

    Deleting c.4 (``A``) gives ``ATG TAG CCC``: the terminator is the second
    codon of the mutant protein, so the answer is 2 -- note it is *not* 4/3
    rounded, and not the reference codon number, which is where an off-by-one
    hides. Inserting a ``G`` before c.4 gives ``ATG GAT AGC CC``, which reaches
    the end of the available sequence without a terminator, and ``None`` is the
    honest answer rather than the last codon.
    """
    index = build_stop_index("ATGATAGCCC")

    result = ptc_after_indel(index, cds_position, deleted, inserted)

    assert result == expected


@pytest.mark.parametrize(
    ("cds_position", "deleted"),
    [(5, 2), (6, 1)],
    ids=["two-base-deletion", "one-base-deletion"],
)
def test_a_hybrid_codon_that_runs_past_the_end_of_the_sequence_reports_no_terminator(
    cds_position, deleted
):
    """Deleting at the 3' edge leaves a partial codon, which is not a stop.

    ``ATGAAA`` minus c.5_6 is ``ATGA``: one complete codon and one orphan base.
    The frame has genuinely run off the end of everything the transcript
    provides, so the only truthful answer is "no terminator found". Translating
    the partial codon instead -- or padding it -- would invent a PTC position,
    and PVS1 reads exactly that position.
    """
    index = build_stop_index("ATGAAA")

    result = ptc_after_indel(index, cds_position, deleted)

    assert result is None


# ----------------------------------------------------------------------
# 3'-most normalisation
# ----------------------------------------------------------------------


@pytest.mark.parametrize("length", [0, -1, -3], ids=["zero", "minus-one", "minus-a-codon"])
def test_a_deletion_of_no_bases_is_rejected_instead_of_looping(length):
    """A non-positive length is a caller bug, and silence would be worse than an error.

    With ``length = 0`` the shift condition compares a base with itself, so the
    loop would run until it hit ``limit`` and return a position unrelated to the
    input.
    """
    with pytest.raises(ValueError, match="deletion length must be positive"):
        normalize_deletion_3prime("AAAAAA", 1, length, limit=6)


# ----------------------------------------------------------------------
# ProteinChange
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("change", "expected_last"),
    [
        (ProteinChange(first_residue=10, ref_segment="AAA", alt_segment=""), 12),
        (ProteinChange(first_residue=10, ref_segment="A", alt_segment=""), 10),
        (ProteinChange(first_residue=10, ref_segment="AA", alt_segment="Q"), 11),
        (ProteinChange(first_residue=10, ref_segment="", alt_segment="QQ"), 10),
    ],
    ids=["three-residues", "one-residue", "delins", "pure-insertion"],
)
def test_the_last_changed_residue_spans_the_reference_segment_not_the_replacement(
    change, expected_last
):
    """``p.Ala10_Ala12del`` names reference residues; the alternate never sets the span.

    The pure-insertion row is the degenerate one: with no reference residues
    removed there is no last residue to name, and the span collapses onto the
    first rather than running backwards to residue 9. A ``last < first`` here
    would render an inverted range in HGVS.
    """
    result = change.last_residue

    assert result == expected_last


@pytest.mark.parametrize(
    ("ref_protein", "alt_protein", "expected"),
    [
        ("MAAA*", "MAA*", ProteinChange(first_residue=4, ref_segment="A", alt_segment="")),
        ("MLPQ*", "MLQ*", ProteinChange(first_residue=3, ref_segment="P", alt_segment="")),
        ("MLPQ*", "MLRQ*", ProteinChange(first_residue=3, ref_segment="P", alt_segment="R")),
    ],
    ids=["repeat-run", "unique-residue", "delins"],
)
def test_a_deletion_inside_a_repeat_run_names_the_last_copy_not_the_first(
    ref_protein, alt_protein, expected
):
    """``MAAA*`` -> ``MAA*`` is ``p.Ala4del``, the third alanine.

    Any of the three alanines could be called the deleted one; HGVS fixes the
    choice as the 3'-most, and trimming the common *prefix* first is what
    implements it. Trimming the suffix first would name residue 2 -- an
    arithmetically defensible description that matches no curated record.
    """
    result = protein_change(ref_protein, alt_protein)

    assert result == expected


def test_two_identical_peptides_leave_the_residue_number_past_the_end_of_the_protein():
    """Expectation corrected: I predicted residue 5, the terminator; it is 6.

    With nothing to trim, the common prefix consumes the whole peptide and
    ``first_residue`` lands one past the last position -- ``len("MLPQ*") + 1``.
    That is a sentinel, not a residue, and it is only safe because the caller
    tests ``is_silent`` before reading it. Pinned so the sentinel cannot quietly
    become residue 5, which *is* a real position and would be rendered.
    """
    result = protein_change("MLPQ*", "MLPQ*")

    assert result == ProteinChange(first_residue=6, ref_segment="", alt_segment="")
    assert result.is_silent is True


# ----------------------------------------------------------------------
# annotate_inframe_deletion guards
# ----------------------------------------------------------------------


@pytest.mark.parametrize("length", [1, 2, 4, 5], ids=["one", "two", "four", "five"])
def test_a_deletion_that_is_not_a_whole_number_of_codons_is_not_an_inframe_deletion(length):
    """A 4 nt deletion shifts the frame, and would be scored under the wrong criteria.

    It is a frameshift, enumerated as a PTC class elsewhere. Annotating it here
    would emit ``p.Xaa10_Xaa11del`` for a variant whose real effect is a
    premature terminator somewhere downstream.
    """
    transcript = _transcript_from_cds("ATGGCTGCTGCTGCTTAA")

    with pytest.raises(
        ValueError, match=rf"in-frame deletion length must be a multiple of 3, got {length}"
    ):
        annotate_inframe_deletion(transcript, 4, length)


@pytest.mark.parametrize(
    ("cds_start", "length", "span"),
    [
        (0, 3, "0..2"),
        (-2, 3, "-2..0"),
        (14, 3, "14..16"),
        (16, 3, "16..18"),
        (13, 6, "13..18"),
    ],
    ids=[
        "before-the-cds",
        "negative",
        "straddling-the-terminator",
        "the-terminator-itself",
        "six-nt-overrunning",
    ],
)
def test_a_deletion_reaching_the_terminator_codon_is_refused_as_a_different_variant_class(
    cds_start, length, span
):
    """Deleting the stop codon is a stop-loss event, scored on a different scale.

    The permitted window ends at ``cds_length - 3`` for exactly that reason: a
    deletion that consumes the terminator produces a protein extension, not a
    shortened one, and describing it as ``p.Xaa5_Ter6del`` would hand a curator
    an in-frame deletion where there is none.
    """
    transcript = _transcript_from_cds("ATGGCTGCTGCTGCTTAA")

    with pytest.raises(
        ValueError, match=rf"in-frame deletion {span} falls outside CDS positions 1\.\.15"
    ):
        annotate_inframe_deletion(transcript, cds_start, length)


def test_the_codon_immediately_before_the_terminator_is_still_deletable():
    """The boundary of the window is inclusive, and ``c.13_15del`` is a real variant.

    ``ATG GCT GCT GCT GCT TAA`` loses its fifth codon, so the protein goes from
    ``MAAAA`` to ``MAAA``. The 3'-most rule names the last alanine, residue 5;
    the deletion does not slide further because position 16 is the terminator
    and the window forbids it.
    """
    transcript = _transcript_from_cds("ATGGCTGCTGCTGCTTAA")

    variant = annotate_inframe_deletion(transcript, 13, 3)

    assert (variant.hgvs_c, variant.hgvs_p) == ("c.13_15del", "p.Ala5del")
    assert variant.cds_position == 13
    assert variant.grch38_pos == "chr1-1012-GCT-del"


def test_a_deletion_downstream_of_a_reference_terminator_gets_no_protein_description():
    """A CDS that already stops early makes the deletion invisible to the protein.

    ``ATG TAA GCT GCT CCC TAA`` terminates at codon 2, so ``c.10_12del`` -- which
    normalises 3'-ward out of the ``GCTGCT`` repeat -- removes a codon that is
    never translated. Reference and alternate peptides are both ``M``, and the
    annotator reports no ``p.`` at all rather than inventing ``p.Ala4del`` for a
    residue the ribosome never reaches.

    This is the only route to that arm: for a well-formed CDS the alternate
    peptide is always one residue shorter than the reference, so an in-frame
    deletion cannot be silent. Reaching it means the reference is broken, which
    is worth reporting -- see the note about the missing marker attribute.
    """
    transcript = _transcript_from_cds("ATGTAAGCTGCTCCCTAA")

    variant = annotate_inframe_deletion(transcript, 7, 3)

    assert variant.hgvs_c == "c.10_12del"
    assert variant.hgvs_p is None
    assert variant.consequence_terms == (Consequence.INFRAME_DELETION,)
    assert variant.codon_index == 4
