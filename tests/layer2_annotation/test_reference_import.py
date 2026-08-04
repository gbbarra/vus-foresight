"""Layer 2: the reference importer.

The importer is the only path by which coordinates enter the system, so its
refusals matter as much as its successes: a config that disagrees with the FASTA
must fail loudly rather than build a transcript that is internally consistent
and wrong.

The round trip is exercised against a synthetic gene, whose true coordinates are
known, so the test has an exact oracle without needing MANE Select.
"""

from __future__ import annotations

import pytest
import yaml

from vus_foresight.genome.importer import import_reference, locate_cds, read_exon_table
from vus_foresight.genome.reference import (
    ReferenceUnavailable,
    build_transcript,
    load_gene_config,
    read_fasta,
)


def _write_resources(tmp_path, gene, name="TOY1"):
    transcript = gene.transcript
    fasta = tmp_path / f"{name}.fa"
    fasta.write_text(
        f">{transcript.transcript_id}\n"
        + "\n".join(transcript.sequence[i : i + 60] for i in range(0, len(transcript.sequence), 60))
        + "\n",
        encoding="ascii",
    )
    exons = tmp_path / f"{name}_exons.tsv"
    exons.write_text(
        "label\tstart\tend\n"
        + "".join(f"{e.label}\t{e.start}\t{e.end}\n" for e in transcript.exons),
        encoding="utf-8",
    )
    return fasta, exons


def _write_config(tmp_path, gene, *, name="TOY1"):
    transcript = gene.transcript
    config = {
        "gene": name,
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
        "sequence": {"path": f"{name}.fa"},
    }
    path = tmp_path / f"{name}.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


@pytest.mark.parametrize("strand", ["+", "-"])
def test_import_round_trips_to_the_original_transcript(tmp_path, both_strands, strand):
    gene = both_strands[strand]
    fasta, exons = _write_resources(tmp_path, gene)
    config_path = _write_config(tmp_path, gene)

    block = import_reference(
        load_gene_config(config_path),
        sequence_path=fasta,
        exon_table_path=exons,
        config_path=config_path,
    )
    assert block["cds_start_tx"] == gene.transcript.cds_start_tx
    assert block["cds_end_tx"] == gene.transcript.cds_end_tx

    rebuilt = build_transcript(load_gene_config(config_path), data_root=tmp_path)
    assert rebuilt.sequence == gene.transcript.sequence
    assert rebuilt.exons == gene.transcript.exons
    assert rebuilt.cds == gene.transcript.cds
    assert rebuilt.strand == gene.transcript.strand


def test_import_records_a_digest_that_later_runs_are_checked_against(tmp_path, minus_gene):
    fasta, exons = _write_resources(tmp_path, minus_gene)
    config_path = _write_config(tmp_path, minus_gene)
    import_reference(
        load_gene_config(config_path),
        sequence_path=fasta,
        exon_table_path=exons,
        config_path=config_path,
    )
    # A silent change to the reference must fail the run, not shift coordinates.
    fasta.write_text(
        ">tampered\n" + "A" * len(minus_gene.transcript.sequence) + "\n", encoding="ascii"
    )
    with pytest.raises(ReferenceUnavailable, match="sha256 mismatch"):
        build_transcript(load_gene_config(config_path), data_root=tmp_path)


def test_exon_table_order_is_taken_as_transcript_order(tmp_path, minus_gene):
    _, exons = _write_resources(tmp_path, minus_gene)
    read = read_exon_table(exons)
    assert read == minus_gene.transcript.exons
    # Minus strand: the coordinates descend, and nothing re-sorts them.
    assert read[0].start > read[1].start


def test_import_refuses_a_mismatched_exon_count(tmp_path, minus_gene):
    fasta, exons = _write_resources(tmp_path, minus_gene)
    config_path = _write_config(tmp_path, minus_gene)
    lines = exons.read_text(encoding="utf-8").splitlines()
    exons.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ReferenceUnavailable, match="exons"):
        import_reference(load_gene_config(config_path), sequence_path=fasta, exon_table_path=exons)


def test_locate_cds_refuses_an_ambiguous_reading_frame():
    # Two ATG..stop frames of the same length: refuse rather than pick.
    ambiguous = "ATGAAATAA" + "ATGAAATAA"
    with pytest.raises(ReferenceUnavailable, match="found 2"):
        locate_cds(ambiguous, 9)
    with pytest.raises(ReferenceUnavailable, match="found 0"):
        locate_cds("ACGTACGTACGT", 9)


def test_a_declared_but_missing_flank_file_degrades_rather_than_fails(tmp_path, minus_gene):
    """Flanks are optional by design, so the loader must not contradict that.

    The intronic enumeration's cardinality does not depend on the reference
    base; without flanks it emits every position and marks each row
    ``reference_base_unknown``. Refusing to enumerate at all would make a
    documented degraded mode unreachable.
    """
    from vus_foresight.genome.reference import SequenceResource, load_flanks
    from vus_foresight.testing import synthetic_gene_config

    config = synthetic_gene_config(minus_gene).model_copy(
        update={"flanks": SequenceResource(path="never_materialised.tsv")}
    )
    assert load_flanks(config, data_root=tmp_path).is_empty
    with pytest.raises(ReferenceUnavailable, match="is missing"):
        load_flanks(config, data_root=tmp_path, required=True)


def test_enumeration_is_complete_without_flanks(minus_gene):
    from vus_foresight.enumeration import enumerate_intronic_snvs, intronic_snv_positions

    positions = intronic_snv_positions(minus_gene.transcript, flank_bp=50)
    variants = list(enumerate_intronic_snvs(minus_gene.transcript, None, flank_bp=50))
    assert len(variants) == 4 * len(positions)
    assert all(v.attributes.get("reference_base_unknown") == "true" for v in variants)


CURATED = """\
# BRCA1 is on the minus strand: exon coordinates descend in transcript order.
gene: TOY1
assembly: TOY
spec: toy_v0.1.0
lof_mechanism: established
transcript:
  id: {transcript_id}
  chrom: {chrom}
  strand: '{strand}'          # not sorted; transcript order is the truth
  cds_length: {cds_length}
  protein_length: {protein_length}
  exon_count: {exon_count}
  # Exon 4 is absent from the reference numbering, hence labels as strings.
  exon_labels: {exon_labels}
sequence:
  path: TOY1.fa
  provenance: synthetic
notes: |
  Legacy numbering is never accepted implicitly.
"""


def _curated_config(tmp_path, gene, name="TOY1"):
    transcript = gene.transcript
    path = tmp_path / f"{name}.yaml"
    path.write_text(
        CURATED.format(
            transcript_id=transcript.transcript_id,
            chrom=transcript.chrom,
            strand=transcript.strand,
            cds_length=transcript.cds_length,
            protein_length=transcript.protein_length,
            exon_count=len(transcript.exons),
            exon_labels="[" + ", ".join(f"'{e.label}'" for e in transcript.exons) + "]",
        ),
        encoding="utf-8",
    )
    return path


def test_write_back_preserves_the_curation_comments(tmp_path, minus_gene):
    """The gene config is documentation as much as data.

    Reserialising it deletes every comment, and the notes it deletes are exactly
    the ones that stop the next reader from "fixing" a descending exon list.
    """
    fasta, exons = _write_resources(tmp_path, minus_gene)
    config_path = _curated_config(tmp_path, minus_gene)

    import_reference(
        load_gene_config(config_path),
        sequence_path=fasta,
        exon_table_path=exons,
        config_path=config_path,
    )
    written = config_path.read_text(encoding="utf-8")

    assert "# BRCA1 is on the minus strand" in written
    assert "# Exon 4 is absent from the reference numbering" in written
    assert "transcript order is the truth" in written
    assert "provenance: synthetic" in written
    assert "Legacy numbering is never accepted implicitly." in written

    # And the derived keys are actually there, read back by the real loader.
    rebuilt = build_transcript(load_gene_config(config_path), data_root=tmp_path)
    assert rebuilt.exons == minus_gene.transcript.exons
    assert rebuilt.cds == minus_gene.transcript.cds


def test_write_back_is_idempotent(tmp_path, minus_gene):
    fasta, exons = _write_resources(tmp_path, minus_gene)
    config_path = _curated_config(tmp_path, minus_gene)

    def run():
        import_reference(
            load_gene_config(config_path),
            sequence_path=fasta,
            exon_table_path=exons,
            config_path=config_path,
        )
        return config_path.read_text(encoding="utf-8")

    first = run()
    assert run() == first, "a second import must not churn the file"


def test_write_back_replaces_a_stale_exon_list_rather_than_appending(tmp_path, minus_gene):
    from vus_foresight.genome.importer import splice_derived

    fasta, exons = _write_resources(tmp_path, minus_gene)
    config_path = _curated_config(tmp_path, minus_gene)
    block = import_reference(
        load_gene_config(config_path),
        sequence_path=fasta,
        exon_table_path=exons,
        config_path=config_path,
    )
    stale = config_path.read_text(encoding="utf-8").replace("start: ", "start: 1", 1)
    respliced = yaml.safe_load(splice_derived(stale, block))
    assert respliced["transcript"]["exons"] == block["exons"]
    assert len(respliced["transcript"]["exons"]) == len(minus_gene.transcript.exons)


def test_splice_refuses_a_layout_it_cannot_edit_safely(tmp_path, minus_gene):
    """Silently reserialising would be worse than failing."""
    from vus_foresight.genome.importer import ImportedReference, splice_derived

    block = ImportedReference(
        {
            "cds_start_tx": 4,
            "cds_end_tx": 9,
            "exons": [{"label": "1", "start": 1, "end": 12}],
            "sequence_sha256": "abc",
            "sequence_length": 12,
        }
    )
    # A flow mapping is valid YAML and outside what a line splice can edit.
    with pytest.raises(ReferenceUnavailable, match="could not be spliced"):
        splice_derived("transcript: {id: NM_1.1, exon_count: 1}\n", block)


def test_read_fasta_refuses_a_multi_record_file(tmp_path):
    path = tmp_path / "two.fa"
    path.write_text(">a\nACGT\n>b\nTGCA\n", encoding="ascii")
    with pytest.raises(ReferenceUnavailable, match="exactly one FASTA record"):
        read_fasta(path)
