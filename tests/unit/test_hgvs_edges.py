"""Edge cases in HGVS string construction.

Oracle: the HGVS Sequence Variant Nomenclature standard, plus the conventions
this module's own docstring fixes (three-letter amino acids, no parentheses
around a predicted effect, ``fsTer`` counting the new codon through the
terminator). Every expectation below is a literal string a VCEP curator could
read off a ClinVar record. That is what makes these tests falsifiable: the way
they can be *wrong* is that the standard was transcribed wrongly here, not that
they were bent to agree with whatever the code already emitted.

Scope note, because it decides what is even testable in this module: nothing
here *computes* a ``c.`` coordinate. Intronic offsets, 5' UTR negatives and
3' UTR ``*n`` anchors are rendered upstream by
:func:`vus_foresight.genome.transcript.format_c_position` and reach these
functions as opaque tokens. What this module can still get wrong -- and what is
pinned below -- is handling a token that is not a bare integer: the
single-position collapse in ``del``/``delins`` compares whole tokens, so a
regression to a numeric or prefix comparison would emit ``c.100+3_100+3del``
or split ``c.*15`` at the asterisk. Those are exactly the "looks right, is
wrong" failures this layer exists to catch.
"""

from __future__ import annotations

import pytest

from vus_foresight.annotate.hgvs import (
    format_deletion,
    format_delins,
    format_exon_cnv,
    format_protein_delins,
    format_protein_frameshift,
    format_protein_stop_lost,
    format_substitution,
)
from vus_foresight.genome.transcript import CPosition, format_c_position


@pytest.mark.parametrize(
    ("position", "expected"),
    [
        (CPosition(base=100, offset=3), "c.100+3G>T"),
        (CPosition(base=101, offset=-2), "c.101-2G>T"),
        (CPosition(base=-20), "c.-20G>T"),
        (CPosition(base=-20, offset=1), "c.-20+1G>T"),
        (CPosition(base=15, utr3=True), "c.*15G>T"),
        (CPosition(base=15, utr3=True, offset=-2), "c.*15-2G>T"),
    ],
    ids=["donor", "acceptor", "utr5", "utr5-intronic", "utr3", "utr3-intronic"],
)
def test_a_substitution_carries_a_noncoding_anchor_through_untouched(position, expected):
    """A splice or UTR variant is described by its anchor; mangling it moves the variant."""
    token = format_c_position(position)

    result = format_substitution(token, "G", "T")

    assert result == expected


@pytest.mark.parametrize(
    ("start", "end", "expected"),
    [
        (CPosition(base=100, offset=3), CPosition(base=100, offset=3), "c.100+3del"),
        (CPosition(base=100, offset=3), CPosition(base=100, offset=5), "c.100+3_100+5del"),
        (CPosition(base=15, utr3=True), CPosition(base=15, utr3=True), "c.*15del"),
        (CPosition(base=15, utr3=True), CPosition(base=17, utr3=True), "c.*15_*17del"),
        (CPosition(base=-20), CPosition(base=-20), "c.-20del"),
        (CPosition(base=-20), CPosition(base=-18), "c.-20_-18del"),
        (CPosition(base=-3), CPosition(base=2), "c.-3_2del"),
    ],
    ids=[
        "intronic-single",
        "intronic-span",
        "utr3-single",
        "utr3-span",
        "utr5-single",
        "utr5-span",
        "across-the-start-codon",
    ],
)
def test_a_deletion_uses_a_range_only_when_its_two_ends_are_different_positions(
    start, end, expected
):
    """``c.*15_*15del`` is not valid HGVS, and would not match any curated record."""
    first, last = format_c_position(start), format_c_position(end)

    result = format_deletion(first, last)

    assert result == expected


@pytest.mark.parametrize(
    ("start", "end", "inserted", "expected"),
    [
        (CPosition(base=100, offset=3), CPosition(base=100, offset=3), "A", "c.100+3delinsA"),
        (
            CPosition(base=100, offset=3),
            CPosition(base=100, offset=4),
            "AT",
            "c.100+3_100+4delinsAT",
        ),
        (CPosition(base=15, utr3=True), CPosition(base=15, utr3=True), "AT", "c.*15delinsAT"),
        (CPosition(base=-20), CPosition(base=-18), "G", "c.-20_-18delinsG"),
    ],
    ids=["intronic-single", "intronic-span", "utr3-single", "utr5-span"],
)
def test_a_delins_uses_a_range_only_when_its_two_ends_are_different_positions(
    start, end, inserted, expected
):
    """Same collapse rule as ``del``; the two must not drift apart."""
    first, last = format_c_position(start), format_c_position(end)

    result = format_delins(first, last, inserted)

    assert result == expected


@pytest.mark.parametrize(
    ("first_aa", "first_residue", "last_aa", "last_residue", "inserted", "expected"),
    [
        ("L", 100, "L", 100, "Q", "p.Leu100delinsGln"),
        ("L", 100, "L", 100, "QH", "p.Leu100delinsGlnHis"),
        ("L", 100, "P", 101, "Q", "p.Leu100_Pro101delinsGln"),
        ("L", 100, "P", 102, "QH", "p.Leu100_Pro102delinsGlnHis"),
        ("M", 1, "M", 1, "*", "p.Met1delinsTer"),
    ],
    ids=[
        "one-residue-one-aa",
        "one-residue-two-aa",
        "two-residues",
        "three-residues",
        "terminator",
    ],
)
def test_a_protein_delins_uses_a_residue_range_only_when_first_and_last_differ(
    first_aa, first_residue, last_aa, last_residue, inserted, expected
):
    """``p.Leu100_Leu100delinsGln`` would be a self-referential range, not a description.

    The single-residue arm is the one an in-frame deletion of exactly one codon
    that also changes its neighbour reaches, which is a real annotator path.
    """
    result = format_protein_delins(first_aa, first_residue, last_aa, last_residue, inserted)

    assert result == expected


@pytest.mark.parametrize(
    ("codons_to_stop", "expected"),
    [
        (12, "p.Arg100SerfsTer12"),
        (1, "p.Arg100SerfsTer1"),
        (None, "p.Arg100SerfsTer?"),
    ],
    ids=["known", "immediate-terminator", "unknown"],
)
def test_a_frameshift_reports_an_unknown_terminator_distance_as_a_question_mark(
    codons_to_stop, expected
):
    """Guessing a number here would put a fabricated PTC distance in front of a curator."""
    result = format_protein_frameshift("R", 100, "S", codons_to_stop)

    assert result == expected


def test_a_frameshift_whose_new_residue_is_unknown_is_named_without_a_terminator_distance():
    """The short ``p.Arg100fs`` form, per the module docstring.

    Reaching this arm requires passing ``None`` for ``alt_aa``, which the
    parameter's own annotation (``alt_aa: str``) forbids -- see the report. The
    behaviour is documented and is pinned here as documented; the annotation is
    reported as the defect rather than repaired.

    Note the terminator distance is discarded even when the caller supplies one:
    without the new residue there is no ``fsTer`` form to hang it on.
    """
    result = format_protein_frameshift("R", 100, None, 12)

    assert result == "p.Arg100fs"


def test_a_stop_loss_names_the_new_residue_and_leaves_the_extension_length_unknown():
    """The extension needs the 3' UTR reading frame, which this pipeline does not read."""
    result = format_protein_stop_lost(1864, "L")

    assert result == "p.Ter1864LeuextTer?"


@pytest.mark.parametrize(
    ("first_label", "last_label", "kind", "expected"),
    [
        ("3", "3", "deletion", "BRCA1 exon 3 deletion"),
        ("3", "5", "deletion", "BRCA1 exon 3-5 deletion"),
        ("3", "5", "duplication", "BRCA1 exon 3-5 duplication"),
        ("11A", "11A", "duplication", "BRCA1 exon 11A duplication"),
    ],
    ids=["single", "range", "range-duplication", "non-numeric-label"],
)
def test_an_exon_cnv_names_a_range_only_when_it_spans_more_than_one_exon(
    first_label, last_label, kind, expected
):
    """Clinical exon labels are not always integers, so the span rule compares labels."""
    result = format_exon_cnv("BRCA1", first_label, last_label, kind)

    assert result == expected
