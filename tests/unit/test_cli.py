"""The command line, which had no test at all.

Local to one module, so it lives here rather than in a layer: the question is
not "is this classification right" -- the layers answer that -- but "does the
interface do what it says". Which is its own risk. The CLI is where a curator
meets the system, and its refusals are load-bearing: refusing to write a map
from uncurated thresholds, refusing a snapshot built against the wrong
transcript version, refusing to guess when a required file is missing. A
refusal that silently stopped working would not fail any other test here.

Characterisation, since the code already exists (CLAUDE.md section 2): each
test is written against the behaviour expected of the command, and a
disagreement is treated as a finding to adjudicate rather than an expectation
to adjust.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from vus_foresight.adapters.clinvar import ClinVarRecord
from vus_foresight.adapters.clinvar_import import write_snapshot
from vus_foresight.cli import app

REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC_DIR = REPO_ROOT / "config" / "specs"

runner = CliRunner()


def _run(*args: str):
    return runner.invoke(app, list(args))


@pytest.fixture()
def gene_on_disk(tmp_path, minus_gene):
    """A synthetic gene config plus its FASTA, ready for the commands to read.

    Written out rather than taken from a fixture object because the commands
    take paths: what is under test includes the loading, not only the maths.
    """
    transcript = minus_gene.transcript
    fasta = tmp_path / "TOY1.fa"
    fasta.write_text(f">{transcript.transcript_id}\n{transcript.sequence}\n", encoding="ascii")
    config = {
        "gene": "TOY1",
        "assembly": "TOY",
        "spec": "toy_v0.1.0",
        "lof_mechanism": "established",
        "transcript": {
            "id": transcript.transcript_id,
            "chrom": transcript.chrom,
            "strand": transcript.strand,
            "cds_length": transcript.cds_length,
            "protein_length": transcript.protein_length,
            "exon_count": len(transcript.exons),
            "exon_labels": [e.label for e in transcript.exons],
            "cds_start_tx": transcript.cds_start_tx,
            "cds_end_tx": transcript.cds_end_tx,
            "exons": [{"label": e.label, "start": e.start, "end": e.end} for e in transcript.exons],
        },
        "sequence": {"path": "TOY1.fa"},
    }
    path = tmp_path / "TOY1.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


@pytest.fixture()
def unpinned_gene(tmp_path, minus_gene):
    """The state a fresh clone is in: a config with no coordinates imported."""
    transcript = minus_gene.transcript
    config = {
        "gene": "TOY2",
        "assembly": "TOY",
        "spec": "toy_v0.1.0",
        "lof_mechanism": "established",
        "transcript": {
            "id": transcript.transcript_id,
            "chrom": transcript.chrom,
            "strand": transcript.strand,
            "cds_length": transcript.cds_length,
            "protein_length": transcript.protein_length,
            "exon_count": len(transcript.exons),
            "exon_labels": [e.label for e in transcript.exons],
        },
        "sequence": {"path": "missing.fa"},
    }
    path = tmp_path / "TOY2.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# selftest: the whole pipeline with no reference data at all.
# --------------------------------------------------------------------------


def test_selftest_runs_the_pipeline_on_a_gene_that_does_not_exist(tmp_path):
    out = tmp_path / "demo.parquet"
    result = _run("selftest", "--spec-dir", str(SPEC_DIR), "--out", str(out))

    assert result.exit_code == 0, result.output
    assert "RESOLVED_NOT_BLOCKED" in result.output
    assert out.exists()


def test_selftest_needs_no_network_and_no_reference(tmp_path):
    """The domain-isolation claim, executable: it works on a fictional gene."""
    result = _run("selftest", "--spec-dir", str(SPEC_DIR))

    assert result.exit_code == 0, result.output
    for token in ("BRCA1", "BRCA2", "ENIGMA"):
        assert token not in result.output


# --------------------------------------------------------------------------
# enumerate
# --------------------------------------------------------------------------


def test_enumerate_reports_the_closed_form_counts(gene_on_disk, tmp_path, minus_gene):
    result = _run("enumerate", "--gene", str(gene_on_disk), "--data-root", str(tmp_path))

    assert result.exit_code == 0, result.output
    assert f"{3 * minus_gene.transcript.cds_length:,}" in result.output
    # 54 per codon, terminator codon included -- so n_codons, not the
    # residue count. The first version of this assertion used
    # protein_length and was wrong by exactly one codon.
    assert f"{54 * minus_gene.transcript.n_codons:,}" in result.output


def test_enumerate_writes_one_row_per_variant(gene_on_disk, tmp_path):
    out = tmp_path / "variants.tsv"
    result = _run(
        "enumerate", "--gene", str(gene_on_disk), "--data-root", str(tmp_path), "--out", str(out)
    )

    assert result.exit_code == 0, result.output
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "variant_id",
        "hgvs_c",
        "hgvs_p",
        "consequence",
        "kind",
        "grch38_pos",
    ]
    assert f"wrote {len(lines) - 1:,} rows" in result.output


def test_a_gene_without_imported_coordinates_exits_two_and_says_what_to_run(
    unpinned_gene, tmp_path
):
    """Not a crash: an expected state with an actionable message."""
    result = _run("enumerate", "--gene", str(unpinned_gene), "--data-root", str(tmp_path))

    assert result.exit_code == 2
    assert "reference import" in result.output


# --------------------------------------------------------------------------
# map: the refusals matter more than the happy path.
# --------------------------------------------------------------------------


def test_map_refuses_an_unverified_specification(tmp_path, minus_gene, gene_on_disk):
    """An uncurated threshold must not silently become a published number."""
    config = yaml.safe_load(gene_on_disk.read_text(encoding="utf-8"))
    config["spec"] = "enigma_brca_v1.1.0"
    gene_on_disk.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    result = _run(
        "map",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(tmp_path / "out"),
    )

    assert result.exit_code != 0
    assert "verified: false" in result.output
    assert "--allow-unverified" in result.output
    assert not (tmp_path / "out").exists()


def test_map_with_no_enumeration_tier_selected_is_refused(tmp_path, gene_on_disk):
    result = _run(
        "map",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(tmp_path / "out"),
        "--tiers",
        "",
    )

    assert result.exit_code != 0
    assert "no enumeration tiers selected" in result.output


def test_map_writes_a_partition_and_reuses_by_signature(tmp_path, gene_on_disk):
    out_dir = tmp_path / "out"
    result = _run(
        "map",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(out_dir),
        "--tiers",
        "1",
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "gene=TOY1" / "gap_map.parquet").exists()
    assert "distinct evidence profiles" in result.output


def test_the_copy_number_map_is_written_to_its_own_file(tmp_path, gene_on_disk):
    """A different scale. Mixing it into the sequence rows would be silent."""
    out_dir = tmp_path / "out"
    result = _run(
        "map",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(out_dir),
        "--tiers",
        "2",
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "gene=TOY1" / "gap_map_cnv.parquet").exists()
    assert "a different scale" in result.output


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------


@pytest.fixture()
def written_map(tmp_path, gene_on_disk):
    out_dir = tmp_path / "out"
    result = _run(
        "map",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(out_dir),
        "--tiers",
        "1",
    )
    assert result.exit_code == 0, result.output
    return out_dir / "gene=TOY1" / "gap_map.parquet"


def test_report_writes_the_four_aggregations(tmp_path, written_map):
    out_dir = tmp_path / "results"
    result = _run("report", str(written_map), "--out-dir", str(out_dir), "--top", "5")

    assert result.exit_code == 0, result.output
    assert sorted(p.name for p in out_dir.iterdir()) == [
        "available_uningested.tsv",
        "blocking_summary.tsv",
        "consequence_class.tsv",
        "equivalence_top5.tsv",
    ]
    header = (out_dir / "blocking_summary.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert header.split("\t") == ["gene", "blocking_reason", "variants"]


def test_the_blocking_summary_accounts_for_every_row(tmp_path, written_map):
    """The aggregation is the product; a row it drops is a row nobody sees."""
    out_dir = tmp_path / "results"
    _run("report", str(written_map), "--out-dir", str(out_dir))

    lines = (out_dir / "blocking_summary.tsv").read_text(encoding="utf-8").splitlines()[1:]
    total = sum(int(line.split("\t")[2]) for line in lines)
    from vus_foresight.output import read_parquet

    assert total == read_parquet(written_map).height


# --------------------------------------------------------------------------
# timeline
# --------------------------------------------------------------------------


def test_timeline_needs_at_least_two_snapshots(written_map):
    result = _run("timeline", f"2020-01-01={written_map}")

    assert result.exit_code != 0
    assert "at least two snapshots" in result.output


def test_timeline_refuses_an_argument_that_is_not_label_equals_path(written_map):
    result = _run("timeline", str(written_map), f"2024-01-01={written_map}")

    assert result.exit_code != 0
    assert "expected 'label=path'" in result.output


def test_timeline_over_an_unchanged_pair_reports_no_transitions(tmp_path, written_map):
    out = tmp_path / "transitions.tsv"
    result = _run(
        "timeline",
        f"2020-01-01={written_map}",
        f"2024-01-01={written_map}",
        "--out",
        str(out),
    )

    assert result.exit_code == 0, result.output
    assert "wrote 0 transitions" in result.output
    assert out.read_text(encoding="utf-8").splitlines() == [
        "\t".join(  # noqa: FLY002 -- a column list diffs one line per column
            [
                "label_before",
                "label_after",
                "variant_id",
                "hgvs_p",
                "consequence",
                "class_before",
                "class_after",
                "points_before",
                "points_after",
                "blocking_before",
                "blocking_after",
                "criteria_gained",
                "criteria_lost",
                "cause",
            ]
        )
    ]


# --------------------------------------------------------------------------
# validate
# --------------------------------------------------------------------------


def test_validate_without_outcomes_or_a_second_map_is_refused(written_map):
    result = _run("validate", "--map-at-t", str(written_map))

    assert result.exit_code != 0
    assert "--outcomes" in result.output


def test_validate_in_time_series_mode_needs_the_second_snapshot(tmp_path, written_map):
    snapshot = tmp_path / "cv.tsv"
    write_snapshot([], snapshot)
    result = _run(
        "validate",
        "--map-at-t",
        str(written_map),
        "--map-at-t-plus-n",
        str(written_map),
        "--clinvar-at-t",
        str(snapshot),
    )

    assert result.exit_code != 0
    assert "--clinvar-at-t-plus-n is required" in result.output


def test_validate_writes_the_study_as_tsv(tmp_path, written_map, minus_gene):
    before = tmp_path / "cv_before.tsv"
    after = tmp_path / "cv_after.tsv"
    write_snapshot([ClinVarRecord("c.10A>G", "p.Thr4Ala", "uncertain", 2, 4)], before)
    write_snapshot([ClinVarRecord("c.10A>G", "p.Thr4Ala", "pathogenic", 3, 4)], after)

    out_dir = tmp_path / "study"
    result = _run(
        "validate",
        "--map-at-t",
        str(written_map),
        "--map-at-t-plus-n",
        str(written_map),
        "--clinvar-at-t",
        str(before),
        "--clinvar-at-t-plus-n",
        str(after),
        "--spec-dir",
        str(SPEC_DIR),
        "--out-dir",
        str(out_dir),
    )

    assert result.exit_code == 0, result.output
    metrics = dict(
        line.split("\t")
        for line in (out_dir / "metrics.tsv").read_text(encoding="utf-8").splitlines()[1:]
    )
    assert metrics["resolved"] == "1"
    assert metrics["considered"] == "1"


# --------------------------------------------------------------------------
# clinvar build
# --------------------------------------------------------------------------


def _variant_summary(tmp_path, *, transcript: str, rows: int = 1) -> Path:
    header = [
        "#AlleleID",
        "Type",
        "Name",
        "GeneID",
        "GeneSymbol",
        "HGNC_ID",
        "ClinicalSignificance",
        "ClinSigSimple",
        "LastEvaluated",
        "RS# (dbSNP)",
        "nsv/esv (dbVar)",
        "RCVaccession",
        "PhenotypeIDS",
        "PhenotypeList",
        "Origin",
        "OriginSimple",
        "Assembly",
        "ChromosomeAccession",
        "Chromosome",
        "Start",
        "Stop",
        "ReferenceAllele",
        "AlternateAllele",
        "Cytogenetic",
        "ReviewStatus",
        "NumberSubmitters",
        "Guidelines",
        "TestedInGTR",
        "OtherIDs",
        "SubmitterCategories",
        "VariationID",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    lines = ["\t".join(header)]
    for i in range(rows):
        record = dict.fromkeys(header, "")
        record["Name"] = f"{transcript}(TOY1):c.{10 + i}A>G (p.Thr4Ala)"
        record["ClinicalSignificance"] = "Pathogenic"
        record["Assembly"] = "GRCh38"
        record["ReviewStatus"] = "criteria provided, multiple submitters, no conflicts"
        lines.append("\t".join(record[k] for k in header))
    path = tmp_path / "variant_summary.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_clinvar_build_needs_exactly_one_of_transcript_or_gene(tmp_path, gene_on_disk):
    source = _variant_summary(tmp_path, transcript="NM_999002.1")
    out = tmp_path / "snap.tsv"

    neither = _run("clinvar", "build", "--source", str(source), "--out", str(out))
    assert neither.exit_code != 0
    assert "exactly one" in neither.output

    both = _run(
        "clinvar",
        "build",
        "--source",
        str(source),
        "--out",
        str(out),
        "--transcript",
        "NM_999002.1",
        "--gene",
        str(gene_on_disk),
    )
    assert both.exit_code != 0
    assert "exactly one" in both.output


def test_clinvar_build_keeps_records_for_the_declared_transcript(tmp_path):
    source = _variant_summary(tmp_path, transcript="NM_999002.1", rows=3)
    out = tmp_path / "snap.tsv"
    result = _run(
        "clinvar",
        "build",
        "--source",
        str(source),
        "--out",
        str(out),
        "--transcript",
        "NM_999002.1",
    )

    assert result.exit_code == 0, result.output
    assert len(out.read_text(encoding="utf-8").splitlines()) == 4  # header + 3


def test_a_snapshot_that_keeps_nothing_exits_two_rather_than_writing_an_empty_file(tmp_path):
    """The usual symptom of a transcript-version mismatch, and it must be loud.

    An empty snapshot is indistinguishable downstream from "this variant has no
    neighbours", so a silent one would quietly disable PS1 and PM5 for a whole
    gene and still produce a map that looks complete.
    """
    source = _variant_summary(tmp_path, transcript="NM_999002.3")
    out = tmp_path / "snap.tsv"
    result = _run(
        "clinvar",
        "build",
        "--source",
        str(source),
        "--out",
        str(out),
        "--transcript",
        "NM_999002.4",
    )

    assert result.exit_code == 2
    assert "check the transcript version" in result.output


# --------------------------------------------------------------------------
# data check
# --------------------------------------------------------------------------


def test_data_check_exits_two_when_the_reference_is_absent(unpinned_gene, tmp_path):
    result = _run("data", "check", "--gene", str(unpinned_gene), "--data-root", str(tmp_path))

    assert result.exit_code == 2


def test_data_check_passes_on_a_pinned_gene(gene_on_disk, tmp_path):
    result = _run("data", "check", "--gene", str(gene_on_disk), "--data-root", str(tmp_path))

    assert result.exit_code == 0, result.output


def test_data_check_exits_one_when_a_declared_source_covers_nothing(gene_on_disk, tmp_path):
    """Present and empty is its own state: it looks exactly like "no data" to
    the engine, which is why it is caught before a map rather than during one."""
    empty = tmp_path / "predictor.tsv"
    empty.write_text("grch38_pos\tbayesdel\n", encoding="utf-8")
    result = _run(
        "data",
        "check",
        "--gene",
        str(gene_on_disk),
        "--data-root",
        str(tmp_path),
        "--predictor",
        str(empty),
    )

    assert result.exit_code == 1, result.output
