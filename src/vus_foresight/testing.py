"""Synthetic genes and a ready-made demo run.

Shipped inside the package rather than in ``tests/`` for two reasons: the CLI's
``selftest`` command uses it, and the domain-isolation guarantee of spec section
11 is a property of the *library*, so the fixture that demonstrates it belongs
with the library.

The locus is built genomically and then spliced, rather than the other way
round, so the minus-strand fixture is a genuine reverse complement and not a
plus-strand transcript with a strand flag bolted on. That is what makes the
BRCA1(-)/BRCA2(+) symmetry test a real control.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

from .adapters.tabular import FrequencyAdapter, FunctionalAdapter, PredictorAdapter, SpliceAdapter
from .engine.pipeline import MapRunner, default_registry
from .engine.spec import load_spec
from .genome.reference import FunctionalRegion, GeneConfig, TranscriptConfig
from .genome.sequence import CODON_TABLE, complement, reverse_complement
from .genome.transcript import Exon, FlankSequences, Transcript
from .variant import Variant

__all__ = [
    "TOY_EXON_LENGTHS",
    "TOY_INTRON_LENGTHS",
    "TOY_PROTEIN_LENGTH",
    "TOY_UTR3",
    "TOY_UTR5",
    "SyntheticGene",
    "build_demo_runner",
    "build_synthetic_gene",
    "synthetic_gene_config",
]

#: Transcript-order exon lengths, summing to the transcript length.
TOY_EXON_LENGTHS: tuple[int, ...] = (60, 40, 30, 33, 35)
#: Deliberately mixed. Two introns are shorter than 2 x 50 so the flanking
#: windows overlap and the enumerator's deduplication is genuinely exercised;
#: one is odd-length, which exercises the HGVS mid-intron numbering split.
TOY_INTRON_LENGTHS: tuple[int, ...] = (120, 31, 200, 90)
TOY_UTR5 = 30
TOY_UTR3 = 45
TOY_PROTEIN_LENGTH = 40

_NON_STOP_CODONS = tuple(sorted(c for c, aa in CODON_TABLE.items() if aa != "*"))


@dataclass(frozen=True, slots=True)
class SyntheticGene:
    transcript: Transcript
    flanks: FlankSequences
    #: Plus-strand genomic sequence of the whole locus, keyed by coordinate.
    locus: dict[int, str]


def _build_transcript_sequence(seed: int, protein_length: int) -> str:
    rng = random.Random(seed)
    utr5 = "".join(rng.choice("ACGT") for _ in range(TOY_UTR5))
    codons = ["ATG"]
    # A GCC codon is planted on purpose. GCC(Ala) becomes ACC(Thr) and GTC(Val)
    # under the two single substitutions, but ATC(Ile) under both together --
    # the MNV third-amino-acid case of spec section 11.
    codons.append("GCC")
    while len(codons) < protein_length:
        codons.append(rng.choice(_NON_STOP_CODONS))
    codons.append("TAA")
    utr3 = "".join(rng.choice("ACGT") for _ in range(TOY_UTR3))
    return utr5 + "".join(codons) + utr3


def _exon_intervals(
    strand: str, anchor: int, exon_lengths: tuple[int, ...], intron_lengths: tuple[int, ...]
) -> list[tuple[int, int]]:
    """Genomic intervals in *transcript* order."""
    intervals: list[tuple[int, int]] = []
    cursor = anchor
    for i, length in enumerate(exon_lengths):
        if strand == "+":
            intervals.append((cursor, cursor + length - 1))
            cursor += length
            if i < len(intron_lengths):
                cursor += intron_lengths[i]
        else:
            intervals.append((cursor - length + 1, cursor))
            cursor -= length
            if i < len(intron_lengths):
                cursor -= intron_lengths[i]
    return intervals


def build_synthetic_gene(
    *,
    gene: str = "TOY1",
    transcript_id: str = "NM_999999.1",
    chrom: str = "chr99",
    strand: str = "+",
    seed: int = 20260803,
    protein_length: int = TOY_PROTEIN_LENGTH,
    exon_labels: tuple[str, ...] | None = None,
) -> SyntheticGene:
    """Build a toy gene on either strand.

    ``exon_labels`` defaults to ``1,2,3,5,6`` -- a deliberate discontinuity
    mirroring BRCA1's missing exon 4, so any code tempted to do arithmetic on an
    exon label breaks here rather than in production.
    """
    if strand not in ("+", "-"):
        raise ValueError("strand must be '+' or '-'")
    exon_lengths = TOY_EXON_LENGTHS
    intron_lengths = TOY_INTRON_LENGTHS
    labels = exon_labels or ("1", "2", "3", "5", "6")

    tx_seq = _build_transcript_sequence(seed, protein_length)
    if len(tx_seq) != sum(exon_lengths):
        raise ValueError(f"toy transcript length {len(tx_seq)} != exon total {sum(exon_lengths)}")

    span = sum(exon_lengths) + sum(intron_lengths)
    if strand == "+":
        anchor = 1_000_000
        low, high = anchor - 100, anchor + span + 100
    else:
        anchor = 1_000_000 + span
        low, high = 1_000_000 - 100, anchor + 100

    rng = random.Random(seed + 1)
    locus = {pos: rng.choice("ACGT") for pos in range(low, high + 1)}

    intervals = _exon_intervals(strand, anchor, exon_lengths, intron_lengths)
    exons = tuple(
        Exon(label=label, start=start, end=end)
        for label, (start, end) in zip(labels, intervals, strict=True)
    )

    cursor = 0
    for exon in exons:
        chunk = tx_seq[cursor : cursor + exon.length]
        cursor += exon.length
        genomic_chunk = chunk if strand == "+" else reverse_complement(chunk)
        for offset, base in enumerate(genomic_chunk):
            locus[exon.start + offset] = base

    transcript = Transcript(
        transcript_id=transcript_id,
        gene=gene,
        chrom=chrom,
        strand=strand,
        exons=exons,
        cds_start_tx=TOY_UTR5 + 1,
        cds_end_tx=TOY_UTR5 + (protein_length + 1) * 3,
        sequence=tx_seq,
    )

    # Canonical GT..AG at every junction, written in transcript orientation.
    step = transcript.step
    for intron in transcript.introns:
        for n, base in ((1, "G"), (2, "T"), (intron.length - 1, "A"), (intron.length, "G")):
            genomic = intron.upstream_last_genomic + step * n
            locus[genomic] = base if strand == "+" else complement(base)

    flank_bases = {
        transcript.intron_genomic(intron, n): locus[transcript.intron_genomic(intron, n)]
        for intron in transcript.introns
        for n in range(1, intron.length + 1)
    }
    return SyntheticGene(
        transcript=transcript, flanks=FlankSequences(bases=flank_bases), locus=locus
    )


def synthetic_gene_config(
    synthetic: SyntheticGene,
    *,
    spec: str = "toy_v0.1.0",
    lof_mechanism: str = "established",
) -> GeneConfig:
    """A :class:`GeneConfig` for a synthetic gene, with two functional regions."""
    transcript = synthetic.transcript
    return GeneConfig(
        gene=transcript.gene,
        assembly="TOY",
        spec=spec,
        lof_mechanism=lof_mechanism,
        transcript=TranscriptConfig(
            id=transcript.transcript_id,
            chrom=transcript.chrom,
            strand=transcript.strand,
            cds_start_tx=transcript.cds_start_tx,
            cds_end_tx=transcript.cds_end_tx,
            cds_length=transcript.cds_length,
            protein_length=transcript.protein_length,
            exon_count=len(transcript.exons),
            exon_labels=tuple(e.label for e in transcript.exons),
            exons=transcript.exons,
        ),
        functional_regions=(
            FunctionalRegion(
                name="toy critical domain", start_aa=5, end_aa=12, tags=("critical", "domain")
            ),
            FunctionalRegion(name="toy repeat region", start_aa=25, end_aa=32, tags=("repeat",)),
        ),
    )


def build_demo_runner(
    spec_path: str | Path,
    *,
    strand: str = "-",
    computed_at: datetime | None = None,
    clinvar_snapshot_path: str | Path | None = None,
    clinvar_snapshot_date: date | None = None,
) -> tuple[MapRunner, list[Variant]]:
    """A complete, reference-free run: synthetic gene, toy spec, no snapshots.

    Every file-backed adapter is deliberately absent by default, so the demo
    shows the map in its most informative state -- almost everything
    ``NOT_EVALUABLE`` for want of data, which is precisely the picture the
    ``available_uningested`` report is meant to surface.

    A ClinVar snapshot may be supplied so that running this twice against
    different dates and diffing the results demonstrates the section 4
    observation on a machine with no reference data at all.
    """
    from .adapters.clinvar import ClinVarSnapshotAdapter
    from .enumeration import EnumerationClass, enumerate_all

    synthetic = build_synthetic_gene(strand=strand)
    config = synthetic_gene_config(synthetic)
    spec = load_spec(spec_path)
    extra = []
    if clinvar_snapshot_path is not None:
        extra.append(
            ClinVarSnapshotAdapter.from_path(
                clinvar_snapshot_path, snapshot_date=clinvar_snapshot_date
            )
        )
    runner = MapRunner(
        transcript=synthetic.transcript,
        gene=config,
        spec=spec,
        adapters=default_registry(config, extra=extra),
        computed_at=computed_at or datetime(1970, 1, 1),
        clinvar_snapshot=clinvar_snapshot_date,
    )
    variants = list(
        enumerate_all(
            synthetic.transcript,
            classes=(
                EnumerationClass.CODING_SNV,
                EnumerationClass.INTRONIC_SNV,
                EnumerationClass.FRAMESHIFT_CLASS,
                EnumerationClass.INFRAME_DELETION,
            ),
            flanks=synthetic.flanks,
        )
    )
    return runner, variants


def demo_adapters(synthetic: SyntheticGene, tmp_path: Path) -> list:
    """File-backed adapters over tiny snapshots, for tests that need evidence."""
    frequency = tmp_path / "frequency.tsv"
    frequency.write_text("grch38_pos\tgnomad.faf95_popmax\tabraom.af\n", encoding="utf-8")
    predictor = tmp_path / "predictor.tsv"
    predictor.write_text("grch38_pos\tbayesdel\n", encoding="utf-8")
    splice = tmp_path / "splice.tsv"
    splice.write_text("grch38_pos\tds_max\n", encoding="utf-8")
    functional = tmp_path / "functional.tsv"
    functional.write_text("hgvs_p\tclassification\tscore\tdataset\n", encoding="utf-8")
    return [
        FrequencyAdapter(frequency, "test", default_payload={"gnomad": {"faf95_popmax": 0.0}}),
        PredictorAdapter(predictor, "test"),
        SpliceAdapter(splice, "test"),
        FunctionalAdapter(functional, "test"),
    ]
