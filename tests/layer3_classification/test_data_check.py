"""Layer 3: the pre-flight check on materialised data.

The point is coverage, not presence. A snapshot built against the wrong
transcript version is empty; a region slice with the wrong coordinates covers
nothing. Both reach the engine looking exactly like "this variant has no data",
so both have to be caught before a map is computed rather than after.
"""

from __future__ import annotations

import pytest

from vus_foresight.adapters import FrequencyAdapter, PredictorAdapter
from vus_foresight.datacheck import DataReport, check_reference, check_sources
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.genome.reference import load_gene_config

from ..conftest import CONFIG_DIR, DATA_ROOT


def _sample(gene, size=40):
    return list(enumerate_coding_snvs(gene.transcript))[:size]


def test_a_gene_without_imported_coordinates_reports_the_next_command():
    config = load_gene_config(CONFIG_DIR / "genes" / "BRCA1.yaml")
    status = check_reference(config, data_root=DATA_ROOT)
    assert not status.present
    assert "reference import" in (status.detail or "")


def test_an_imported_gene_reports_its_geometry(minus_config, minus_gene, tmp_path):
    """Built from the synthetic config, whose coordinates are inline."""
    status = check_reference(minus_config, data_root=tmp_path)
    # The synthetic gene declares no sequence resource, so the reference is
    # reported absent rather than half-built.
    assert not status.present


def test_a_populated_source_reports_high_coverage(tmp_path, minus_gene):
    sample = _sample(minus_gene)
    path = tmp_path / "frequency.tsv"
    path.write_text(
        "grch38_pos\tgnomad.faf95_popmax\n"
        + "".join(f"{v.grch38_pos}\t0.0\n" for v in sample),
        encoding="utf-8",
    )
    statuses = check_sources(
        [FrequencyAdapter(path, "test")], minus_gene.transcript, sample
    )
    assert statuses[0].coverage == 1.0
    assert statuses[0].symbol == "ok"
    assert statuses[0].records == len(sample)


def test_a_slice_with_the_wrong_coordinates_is_loud_not_silent(tmp_path, minus_gene):
    """The failure this command exists for.

    Thousands of rows, none of which touch the transcript. Reaching the engine,
    every criterion would come out NOT_EVALUABLE and the map would look like a
    gene nobody has ever measured.
    """
    sample = _sample(minus_gene)
    path = tmp_path / "wrong_region.tsv"
    path.write_text(
        "grch38_pos\tgnomad.faf95_popmax\n"
        + "".join(f"chr99-{9_000_000 + i}-A-G\t0.0\n" for i in range(2000)),
        encoding="utf-8",
    )
    statuses = check_sources(
        [FrequencyAdapter(path, "test")], minus_gene.transcript, sample
    )
    assert statuses[0].present
    assert statuses[0].records == 2000
    assert statuses[0].coverage == 0.0
    assert statuses[0].symbol == "!!"

    report = DataReport(gene="TOYM", reference=statuses[0], sources=statuses, sampled=len(sample))
    assert report.empty_but_present == statuses
    assert "matches none of the sampled variants" in report.summary()


def test_partial_coverage_is_flagged_as_neither_fine_nor_broken(tmp_path, minus_gene):
    sample = _sample(minus_gene, size=40)
    path = tmp_path / "partial.tsv"
    path.write_text(
        "grch38_pos\tbayesdel\n"
        + "".join(f"{v.grch38_pos}\t0.4\n" for v in sample[:8]),
        encoding="utf-8",
    )
    statuses = check_sources(
        [PredictorAdapter(path, "test")], minus_gene.transcript, sample
    )
    assert statuses[0].coverage == pytest.approx(0.2)
    assert statuses[0].symbol == " ~"


def test_a_default_payload_counts_as_coverage_because_it_is_an_assertion(
    tmp_path, minus_gene
):
    """gnomAD not listing a variant is a positive statement, not a gap.

    An adapter carrying a default payload answers for every variant, and that is
    correct: "looked for, not seen" is exactly the evidence PM2 rests on.
    """
    sample = _sample(minus_gene)
    path = tmp_path / "empty.tsv"
    path.write_text("grch38_pos\tgnomad.faf95_popmax\n", encoding="utf-8")

    without = check_sources(
        [FrequencyAdapter(path, "test")], minus_gene.transcript, sample
    )
    assert without[0].coverage == 0.0

    with_default = check_sources(
        [
            FrequencyAdapter(
                path, "test", default_payload={"gnomad": {"faf95_popmax": 0.0}}
            )
        ],
        minus_gene.transcript,
        sample,
    )
    assert with_default[0].coverage == 1.0


def test_the_report_is_usable_only_when_the_reference_is_present(minus_gene, tmp_path):
    absent = check_reference(
        load_gene_config(CONFIG_DIR / "genes" / "BRCA2.yaml"), data_root=DATA_ROOT
    )
    report = DataReport(gene="BRCA2", reference=absent)
    assert not report.usable
    assert "BRCA2" in report.summary()
