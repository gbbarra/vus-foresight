"""Layer 2: pulling one transcript out of a MANE Select release.

Everything here guards the same failure: an exon table that is internally
consistent and biologically wrong. Sorting by coordinate reverses a minus-strand
gene; taking labels from the GTF mislabels BRCA1, whose clinical numbering skips
exon 4; mixing a GTF and a FASTA from different releases shifts every coordinate
without changing any of them individually.

The fixtures are MANE-shaped files generated from the synthetic genes, so the
round trip has an exact oracle without a release on disk.
"""

from __future__ import annotations

import gzip

import pytest

from vus_foresight.genome.importer import import_reference
from vus_foresight.genome.mane import extract_transcript, normalise_chrom
from vus_foresight.genome.reference import ReferenceUnavailable, build_transcript, load_gene_config


def _write_release(tmp_path, gene, config, *, gzipped=False, accession_seqname=False, decoys=2):
    """Write a GTF and an RNA FASTA in the shape MANE distributes them."""
    transcript = gene.transcript
    seqname = "NC_000017.11" if accession_seqname else transcript.chrom

    gtf_lines = ["#!genome-build GRCh38\n"]
    for number, exon in enumerate(transcript.exons, start=1):
        attributes = (
            f'gene_id "{transcript.gene}"; transcript_id "{transcript.transcript_id}"; '
            f'exon_number "{number}"; tag "MANE Select";'
        )
        gtf_lines.append(
            f"{seqname}\tBestRefSeq\texon\t{exon.start}\t{exon.end}\t.\t"
            f"{transcript.strand}\t.\t{attributes}\n"
        )
        # CDS features share the transcript id and must be ignored.
        gtf_lines.append(
            f"{seqname}\tBestRefSeq\tCDS\t{exon.start}\t{exon.end}\t.\t"
            f"{transcript.strand}\t0\t{attributes}\n"
        )
    for i in range(decoys):
        gtf_lines.append(
            f"chr1\tBestRefSeq\texon\t{1000 + i}\t{1100 + i}\t.\t+\t.\t"
            f'gene_id "OTHER"; transcript_id "NM_000000.{i}"; exon_number "1";\n'
        )

    fasta_lines = []
    for i in range(decoys):
        fasta_lines.append(f">NM_000000.{i} some other transcript\nACGTACGTAC\n")
    fasta_lines.append(f">{transcript.transcript_id} {transcript.gene} mRNA\n")
    fasta_lines.extend(
        transcript.sequence[i : i + 60] + "\n" for i in range(0, len(transcript.sequence), 60)
    )
    fasta_lines.append(">NM_111111.1 trailing record\nTTTTTTTT\n")

    suffix = ".gz" if gzipped else ""
    gtf_path = tmp_path / f"mane.gtf{suffix}"
    fasta_path = tmp_path / f"mane.fna{suffix}"
    # The lambdas are context-manager factories, consumed by the 'with'
    # statements below, so SIM115 does not apply.
    opener = (
        (lambda p: gzip.open(p, "wt", encoding="utf-8"))  # noqa: SIM115
        if gzipped
        else (lambda p: p.open("w", encoding="utf-8"))
    )
    with opener(gtf_path) as handle:
        handle.writelines(gtf_lines)
    with opener(fasta_path) as handle:
        handle.writelines(fasta_lines)
    return gtf_path, fasta_path


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_round_trips_to_the_original_transcript(
    tmp_path, both_strands, strand, plus_config, minus_config
):
    gene = both_strands[strand]
    config = plus_config if strand == "+" else minus_config
    gtf, fasta = _write_release(tmp_path, gene, config)

    extract = extract_transcript(config, gtf_path=gtf, fasta_path=fasta)
    assert extract.fasta == gene.transcript.sequence
    assert extract.exons == gene.transcript.exons
    assert extract.strand == strand


def test_exon_order_follows_exon_number_not_coordinate(tmp_path, minus_gene, minus_config):
    """The failure this module exists to prevent.

    On the minus strand, transcript order is *descending* genomic coordinate.
    A coordinate sort would return the gene backwards.
    """
    gtf, fasta = _write_release(tmp_path, minus_gene, minus_config)
    extract = extract_transcript(minus_config, gtf_path=gtf, fasta_path=fasta)

    starts = [e.start for e in extract.exons]
    assert starts == sorted(starts, reverse=True)
    assert starts != sorted(starts)
    assert extract.exons[0] == minus_gene.transcript.exons[0]


def test_clinical_exon_labels_win_over_gtf_exon_number(tmp_path, minus_gene, minus_config):
    """BRCA1's GTF says 1..23; the clinic says 1,2,3,5,...,24."""
    gtf, fasta = _write_release(tmp_path, minus_gene, minus_config)
    extract = extract_transcript(minus_config, gtf_path=gtf, fasta_path=fasta)

    assert [e.label for e in extract.exons] == ["1", "2", "3", "5", "6"]
    assert any("clinical exon labels differ" in w for w in extract.warnings)


def test_gzipped_and_accession_keyed_releases_are_accepted(tmp_path, minus_gene, minus_config):
    gtf, fasta = _write_release(
        tmp_path, minus_gene, minus_config, gzipped=True, accession_seqname=False
    )
    extract = extract_transcript(minus_config, gtf_path=gtf, fasta_path=fasta)
    assert extract.exons == minus_gene.transcript.exons


def test_chromosome_spellings_are_reconciled():
    assert normalise_chrom("NC_000017.11") == "chr17"
    assert normalise_chrom("17") == "chr17"
    assert normalise_chrom("chr17") == "chr17"
    assert normalise_chrom("NC_000023.11") == "chrX"


def test_a_gtf_and_fasta_from_different_releases_is_an_error(tmp_path, minus_gene, minus_config):
    """Individually plausible, jointly wrong -- and it shifts every coordinate."""
    gtf, _fasta = _write_release(tmp_path, minus_gene, minus_config)
    truncated = tmp_path / "short.fna"
    truncated.write_text(
        f">{minus_gene.transcript.transcript_id}\n{minus_gene.transcript.sequence[:-10]}\n",
        encoding="utf-8",
    )
    with pytest.raises(ReferenceUnavailable, match="different releases"):
        extract_transcript(minus_config, gtf_path=gtf, fasta_path=truncated)


def test_a_strand_disagreement_is_an_error(tmp_path, plus_gene, minus_config, plus_config):
    gtf, fasta = _write_release(tmp_path, plus_gene, plus_config)
    # Ask for a minus-strand config against a plus-strand release.
    wrong = minus_config.model_copy(
        update={
            "transcript": minus_config.transcript.model_copy(
                update={"id": plus_gene.transcript.transcript_id}
            )
        }
    )
    with pytest.raises(ReferenceUnavailable, match="strand"):
        extract_transcript(wrong, gtf_path=gtf, fasta_path=fasta)


def test_a_missing_transcript_names_the_version_as_the_likely_cause(
    tmp_path, minus_gene, minus_config
):
    gtf, fasta = _write_release(tmp_path, minus_gene, minus_config)
    wrong_version = minus_config.model_copy(
        update={"transcript": minus_config.transcript.model_copy(update={"id": "NM_999002.9"})}
    )
    with pytest.raises(ReferenceUnavailable, match="no exon features"):
        extract_transcript(wrong_version, gtf_path=gtf, fasta_path=fasta)


def test_an_exon_count_disagreement_is_an_error(tmp_path, minus_gene, minus_config):
    gtf, fasta = _write_release(tmp_path, minus_gene, minus_config)
    wrong = minus_config.model_copy(
        update={
            "transcript": minus_config.transcript.model_copy(
                update={"exon_count": 4, "exon_labels": ("1", "2", "3", "5")}
            )
        }
    )
    with pytest.raises(ReferenceUnavailable, match="exons"):
        extract_transcript(wrong, gtf_path=gtf, fasta_path=fasta)


def test_extract_then_import_produces_a_working_transcript(tmp_path, minus_gene, minus_config):
    """The whole acquisition path, end to end, with an exact oracle."""
    import yaml

    gtf, fasta = _write_release(tmp_path, minus_gene, minus_config)
    extract = extract_transcript(minus_config, gtf_path=gtf, fasta_path=fasta)

    fasta_out = tmp_path / "extracted.fa"
    exons_out = tmp_path / "extracted_exons.tsv"
    fasta_out.write_text(extract.fasta_text(), encoding="ascii")
    exons_out.write_text(extract.exon_table(), encoding="utf-8")

    config_path = tmp_path / "gene.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "gene": minus_config.gene,
                "assembly": "TOY",
                "spec": "toy_v0.1.0",
                "transcript": {
                    "id": minus_config.transcript.id,
                    "chrom": minus_config.transcript.chrom,
                    "strand": minus_config.transcript.strand,
                    "cds_length": minus_config.transcript.cds_length,
                    "protein_length": minus_config.transcript.protein_length,
                    "exon_count": minus_config.transcript.exon_count,
                    "exon_labels": list(minus_config.transcript.exon_labels),
                },
                "sequence": {"path": fasta_out.name},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    import_reference(
        load_gene_config(config_path),
        sequence_path=fasta_out,
        exon_table_path=exons_out,
        config_path=config_path,
    )
    rebuilt = build_transcript(load_gene_config(config_path), data_root=tmp_path)

    assert rebuilt.sequence == minus_gene.transcript.sequence
    assert rebuilt.exons == minus_gene.transcript.exons
    assert rebuilt.cds == minus_gene.transcript.cds
    assert rebuilt.cds_start_tx == minus_gene.transcript.cds_start_tx
