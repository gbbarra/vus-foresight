"""What a snapshot file may look like, and what happens when it does not.

Module-local. Every external dataset enters through this reader, and its
tolerances decide what "no data" means downstream: a cell it drops becomes a
MISSING field, which becomes a criterion that cannot be evaluated, which
becomes a blocking reason on the map. So the parsing rules are not an
implementation detail -- they are the definition of absence.
"""

from __future__ import annotations

import pytest

from vus_foresight.adapters.tabular import FrequencyAdapter, PredictorAdapter


def _write(tmp_path, text: str, name: str = "snapshot.tsv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _lookup(adapter, variant, transcript):
    return adapter.lookup(variant, transcript)


@pytest.fixture()
def variant(minus_gene):
    from vus_foresight.enumeration import enumerate_coding_snvs

    return next(iter(enumerate_coding_snvs(minus_gene.transcript)))


# --------------------------------------------------------------------------
# Refusals.
# --------------------------------------------------------------------------


def test_a_missing_snapshot_is_an_error_rather_than_an_empty_table(tmp_path):
    """Adapters never fetch. An absent file is a setup mistake, not "no data".

    Treating it as an empty table would turn "you forgot to materialise gnomAD"
    into "no variant in this gene has a frequency", which is a map that looks
    complete and is not.
    """
    with pytest.raises(FileNotFoundError, match="materialise the snapshot first"):
        FrequencyAdapter(tmp_path / "never_written.tsv", "test")


def test_a_snapshot_with_no_header_row_is_refused(tmp_path):
    path = _write(tmp_path, "")
    with pytest.raises(ValueError, match="has no header row"):
        FrequencyAdapter(path, "test")


# --------------------------------------------------------------------------
# What counts as absent.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("cell", ["", "NA", "null", ".", "  "])
def test_the_conventional_spellings_of_absent_all_become_absent(
    tmp_path, minus_gene, variant, cell
):
    """A key that is not published is what the engine reads as MISSING."""
    path = _write(tmp_path, f"grch38_pos\tbayesdel\n{variant.grch38_pos}\t{cell}\n")
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert "bayesdel" not in payload


def test_zero_is_a_value_and_not_an_absence(tmp_path, minus_gene, variant):
    """The distinction the whole MISSING sentinel exists for.

    A frequency of exactly zero is evidence -- it is what BA1 and BS1 are
    weighed against. Losing it to a falsy check would silently convert "never
    observed in gnomAD" into "gnomAD was not consulted".
    """
    path = _write(tmp_path, f"grch38_pos\tgnomad.faf95_popmax\n{variant.grch38_pos}\t0.0\n")
    payload = _lookup(FrequencyAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload["gnomad"]["faf95_popmax"] == 0.0


@pytest.mark.parametrize(
    ("cell", "expected"),
    [("true", True), ("FALSE", False), ("7", 7), ("0.001", 0.001), ("abnormal", "abnormal")],
)
def test_each_cell_is_parsed_to_its_narrowest_type(tmp_path, minus_gene, variant, cell, expected):
    path = _write(tmp_path, f"grch38_pos\tvalue\n{variant.grch38_pos}\t{cell}\n")
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload["value"] == expected
    assert type(payload["value"]) is type(expected)


def test_a_dotted_column_becomes_a_nested_key(tmp_path, minus_gene, variant):
    """Which is how a specification addresses 'frequency.gnomad.faf95_popmax'."""
    path = _write(tmp_path, f"grch38_pos\ta.b.c\n{variant.grch38_pos}\t1\n")
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload == {"a": {"b": {"c": 1}}}


# --------------------------------------------------------------------------
# Ragged rows, which is the case flagged when zip's strictness was declared.
# --------------------------------------------------------------------------


def test_a_row_shorter_than_the_header_yields_only_the_cells_it_has(tmp_path, minus_gene, variant):
    """Truncation is deliberate here: the absent columns read as MISSING."""
    path = _write(tmp_path, f"grch38_pos\tone\ttwo\tthree\n{variant.grch38_pos}\t1\n")
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload == {"one": 1}


def test_a_row_longer_than_the_header_silently_drops_its_extra_cells(tmp_path, minus_gene, variant):
    """Characterisation of a real gap, recorded rather than fixed here.

    Extra cells are discarded without a word. A snapshot with a column the
    header does not name is a malformed file, and the reader treats it as a
    well-formed one -- so a column added upstream and not reflected in the
    header vanishes, and the map reports that evidence as missing.

    This is not the phase that changes reader behaviour: doing so would alter
    what every existing snapshot parses to. The test exists so the behaviour is
    a decision on record, and so that changing it has to change this test.
    """
    path = _write(tmp_path, f"grch38_pos\tone\n{variant.grch38_pos}\t1\t2\t3\n")
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload == {"one": 1}


def test_a_comment_row_is_skipped(tmp_path, minus_gene, variant):
    path = _write(
        tmp_path,
        f"grch38_pos\tone\n# a note about provenance\t\n{variant.grch38_pos}\t1\n",
    )
    payload = _lookup(PredictorAdapter(path, "test"), variant, minus_gene.transcript)

    assert payload == {"one": 1}


def test_a_variant_absent_from_the_snapshot_gets_the_default_payload(tmp_path, minus_gene, variant):
    """ "Not in the file" and "in the file with no value" are the same answer."""
    path = _write(tmp_path, "grch38_pos\tone\nchr99-1-A-T\t1\n")
    adapter = PredictorAdapter(path, "test")

    assert _lookup(adapter, variant, minus_gene.transcript) == {}


def test_a_declared_default_is_returned_for_an_unlisted_variant(tmp_path, minus_gene, variant):
    """gnomAD's absence means "not observed", which is a value, not a gap."""
    path = _write(tmp_path, "grch38_pos\tgnomad.faf95_popmax\nchr99-1-A-T\t0.5\n")
    adapter = FrequencyAdapter(
        path, "test", default_payload={"gnomad": {"observed": False, "faf95_popmax": 0.0}}
    )

    payload = _lookup(adapter, variant, minus_gene.transcript)
    assert payload == {"gnomad": {"observed": False, "faf95_popmax": 0.0}}
