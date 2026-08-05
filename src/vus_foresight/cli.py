"""Command line interface.

``vus-foresight --help`` for the list. The commands that matter:

    enumerate   how many variants exist, by class -- the layer-1 invariants
    map         the gap map itself, to Parquet
    report      the aggregations: blocking reasons, available-uningested
    timeline    diff maps computed at different dates, and attribute what moved
    validate    the section 10 protocol against vus-hindsight outcomes
    selftest    the whole pipeline on a synthetic gene, no reference needed
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Optional

import typer

from .adapters import (
    ClinVarSnapshotAdapter,
    FrequencyAdapter,
    FunctionalAdapter,
    PredictorAdapter,
    SpliceAdapter,
)
from .adapters.base import Adapter
from .adapters.clinvar_import import build_snapshot
from .engine.cnv_scoring import CNVScoringConfig, cnv_row
from .engine.pipeline import MapRunner, default_registry
from .engine.spec import load_spec
from .enumeration import (
    EnumerationClass,
    build_frameshift_classes,
    coding_snv_count,
    enumerate_all,
    enumerate_exon_cnvs,
    exon_interval_count,
    mnv_count,
)
from .genome.importer import import_reference
from .genome.reference import (
    ReferenceUnavailable,
    build_transcript,
    load_flanks,
    load_gene_config,
)
from .output import (
    GapMapSource,
    available_uningested_report,
    blocking_summary,
    equivalence_summary,
    read_parquet,
    summarise_counts,
    write_parquet,
)

app = typer.Typer(
    add_completion=False,
    help="Evidence gap map for clinical genomics. No patient data, no LLM, no network.",
)
reference_app = typer.Typer(help="Manage MANE Select reference resources.")
app.add_typer(reference_app, name="reference")
clinvar_app = typer.Typer(help="Build dated ClinVar snapshots for PS1 and PM5.")
app.add_typer(clinvar_app, name="clinvar")
data_app = typer.Typer(help="Inspect the external snapshots. See docs/data-acquisition.md.")
app.add_typer(data_app, name="data")

#: Determinism: the default stamp is a fixed epoch, not the wall clock. A run
#: that wants a real timestamp must say so, because two runs that differ only in
#: their timestamp are not reproducible and the CI check would catch it as noise.
FIXED_EPOCH = datetime(1970, 1, 1, tzinfo=UTC).replace(tzinfo=None)


def _load_reference(gene_path: Path, data_root: Path):
    """Load a gene and its reference, or fail readably.

    A missing reference is an expected state, not a crash: the resources are not
    vendored on purpose (see docs/reference-data.md), so the message has to say
    what to run next rather than print a stack trace at a curator.
    """
    config = load_gene_config(gene_path)
    try:
        transcript = build_transcript(config, data_root=data_root)
        flanks = load_flanks(config, data_root=data_root)
    except ReferenceUnavailable as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from None
    if config.flanks is not None and flanks.is_empty:
        # Degraded, not failed: every intronic position is still enumerated,
        # each marked reference_base_unknown rather than given a guessed allele.
        typer.secho(
            f"  note: {config.flanks.path} is not present; intronic rows will carry "
            "reference_base_unknown",
            fg=typer.colors.YELLOW,
            err=True,
        )
    return config, transcript, flanks


def _load(gene_path: Path, spec_dir: Path, data_root: Path):
    """As above, plus the specification. Only commands that classify need it.

    Enumeration deliberately does not: how many variants exist is a property of
    the transcript, not of anyone's interpretation rules, and coupling the two
    would make counting fail because a specification file was elsewhere.
    """
    config, transcript, flanks = _load_reference(gene_path, data_root)
    spec = load_spec(spec_dir / f"{config.spec}.yaml")
    return config, spec, transcript, flanks


@app.command("enumerate")
def enumerate_command(
    gene: Path = typer.Option(..., help="Path to a gene YAML."),
    data_root: Path = typer.Option(Path("."), help="Root for reference resources."),
    out: Optional[Path] = typer.Option(None, help="Write the enumerated variants as TSV."),
) -> None:
    """Count and optionally dump every possible variant for a gene.

    Needs no specification: how many variants exist is a property of the
    transcript, not of anyone's interpretation rules.
    """
    config, transcript, flanks = _load_reference(gene, data_root)

    counts = {
        "coding_snv": coding_snv_count(transcript),
        "intracodon_mnv": mnv_count(transcript),
        "exon_cnv_intervals_per_direction": exon_interval_count(transcript),
    }
    frameshift = build_frameshift_classes(transcript)
    counts["frameshift_classes_reachable"] = len(frameshift.classes)
    counts["frameshift_codon_positions_unreachable"] = len(frameshift.unreachable)

    typer.echo(f"{config.gene} / {transcript.transcript_id}")
    typer.echo(f"  CDS {transcript.cds_length} nt, {transcript.protein_length} residues")
    for key, value in counts.items():
        typer.echo(f"  {key:44s} {value:>9,}")

    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        written = 0
        with out.open("w", encoding="utf-8") as handle:
            handle.write("variant_id\thgvs_c\thgvs_p\tconsequence\tkind\tgrch38_pos\n")
            for variant in enumerate_all(
                transcript, flanks=flanks, flank_bp=config.intronic_flank_bp
            ):
                handle.write(
                    f"{variant.variant_id}\t{variant.hgvs_c}\t{variant.hgvs_p or ''}\t"
                    f"{variant.consequence.value}\t{variant.kind.value}\t"
                    f"{variant.grch38_pos or ''}\n"
                )
                written += 1
        typer.echo(f"  wrote {written:,} rows to {out}")


@app.command("map")
def map_command(
    gene: Path = typer.Option(..., help="Path to a gene YAML."),
    data_root: Path = typer.Option(Path(".")),
    spec_dir: Path = typer.Option(Path("config/specs")),
    out_dir: Path = typer.Option(Path("out"), help="Directory for the Parquet output."),
    frequency: Optional[Path] = typer.Option(None, help="gnomAD/ABraOM snapshot TSV."),
    predictor: Optional[Path] = typer.Option(None, help="dbNSFP snapshot TSV."),
    splice: Optional[Path] = typer.Option(None, help="SpliceAI snapshot TSV."),
    functional: Optional[Path] = typer.Option(None, help="MAVE/SGE snapshot TSV."),
    clinvar: Optional[Path] = typer.Option(
        None, help="Normalised ClinVar snapshot from 'vus-foresight clinvar build'."
    ),
    clinvar_date: Optional[str] = typer.Option(None, help="ISO date of the ClinVar snapshot."),
    clinvar_min_stars: int = typer.Option(
        1,
        help=(
            "Review-status floor a neighbouring ClinVar record must clear before it "
            "can support PS1 or PM5."
        ),
    ),
    gnomad_version: str = typer.Option("v4", help="Recorded in every row's provenance."),
    computed_at: Optional[str] = typer.Option(
        None, help="ISO timestamp stamped on every row. Defaults to a fixed epoch."
    ),
    tiers: str = typer.Option("1,2", help="Enumeration tiers to include."),
    allow_unverified: bool = typer.Option(
        False, help="Write a map from a specification whose thresholds are uncurated."
    ),
    reuse_by_signature: bool = typer.Option(
        True,
        help=(
            "Evaluate once per distinct evidence profile instead of once per variant. "
            "Turn off to run the reference implementation."
        ),
    ),
) -> None:
    """Compute the gap map for one gene and write it to Parquet."""
    config, spec, transcript, flanks = _load(gene, spec_dir, data_root)

    if not spec.verified and not allow_unverified:
        raise typer.BadParameter(
            f"{spec.spec_version} is marked verified: false -- its numeric thresholds "
            "have not been curated against the published document. Pass "
            "--allow-unverified to produce a demonstration map anyway."
        )

    classes: list[EnumerationClass] = []
    if "1" in tiers:
        classes += [EnumerationClass.CODING_SNV, EnumerationClass.INTRONIC_SNV]
    if "2" in tiers:
        classes += [
            EnumerationClass.INTRACODON_MNV,
            EnumerationClass.FRAMESHIFT_CLASS,
            EnumerationClass.INFRAME_DELETION,
            EnumerationClass.EXON_CNV,
        ]
    if not classes:
        raise typer.BadParameter("no enumeration tiers selected")

    extra: list[Adapter] = []
    if frequency is not None:
        extra.append(
            FrequencyAdapter(
                frequency,
                gnomad_version,
                default_payload={"gnomad": {"observed": False, "faf95_popmax": 0.0}},
            )
        )
    if predictor is not None:
        extra.append(PredictorAdapter(predictor, "dbnsfp-snapshot"))
    if splice is not None:
        extra.append(SpliceAdapter(splice, "spliceai-snapshot"))
    if functional is not None:
        extra.append(
            FunctionalAdapter(
                functional,
                "mave-snapshot",
                assayed_regions=tuple(
                    (int(d["start_aa"]), int(d["end_aa"]), str(d["name"]))
                    for d in config.mave_datasets
                ),
            )
        )
    if clinvar is not None:
        extra.append(
            ClinVarSnapshotAdapter.from_path(
                clinvar,
                snapshot_date=date.fromisoformat(clinvar_date) if clinvar_date else None,
                min_stars=clinvar_min_stars,
            )
        )

    stamp = datetime.fromisoformat(computed_at) if computed_at else FIXED_EPOCH
    runner = MapRunner(
        transcript=transcript,
        gene=config,
        spec=spec,
        adapters=default_registry(config, extra=extra),
        computed_at=stamp,
        clinvar_snapshot=date.fromisoformat(clinvar_date) if clinvar_date else None,
        gnomad_version=gnomad_version if frequency is not None else None,
        reuse_by_signature=reuse_by_signature,
    )

    # Streamed, not materialised. Every row carries its full evaluation trace,
    # so holding a whole gene's worth in memory before writing needs tens of
    # gigabytes -- measured, on real BRCA1, before this was a generator.
    variants = enumerate_all(
        transcript,
        classes=tuple(classes),
        flanks=flanks,
        flank_bp=config.intronic_flank_bp,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    sequence_path = out_dir / f"gene={config.gene}" / "gap_map.parquet"
    written = write_parquet((result.row for result in runner.run(variants)), sequence_path)
    typer.echo(f"wrote {written:,} sequence-level rows to {sequence_path}")
    if reuse_by_signature:
        typer.echo(f"  evaluated {runner.cache.summary()}")

    if EnumerationClass.EXON_CNV in classes:
        cnv_config = CNVScoringConfig()
        cnv_rows = [
            cnv_row(
                variant,
                transcript,
                config,
                cnv_config,
                computed_at=stamp,
                clinvar_snapshot=runner.clinvar_snapshot,
                gnomad_version=runner.gnomad_version,
            )
            for variant in enumerate_exon_cnvs(transcript)
        ]
        cnv_path = out_dir / f"gene={config.gene}" / "gap_map_cnv.parquet"
        write_parquet(cnv_rows, cnv_path)
        typer.echo(
            f"wrote {len(cnv_rows):,} copy-number rows to {cnv_path} "
            f"(scored by {cnv_config.version_string}, a different scale)"
        )


@app.command("report")
def report_command(
    parquet: Path = typer.Argument(..., help="A gap map Parquet file."),
    top: int = typer.Option(20, help="Rows to show per section."),
    out_dir: Optional[Path] = typer.Option(
        None,
        help=(
            "Also write the aggregations as TSV. These are small and "
            "deterministic, so they are the part of a run worth committing; the "
            "Parquet they came from is not."
        ),
    ),
) -> None:
    """Print the aggregations the map exists to support."""
    frame = read_parquet(parquet)
    typer.echo(f"{frame.height:,} rows\n")

    blocking = blocking_summary(frame)
    typer.echo("== blocking reason ==")
    typer.echo(str(blocking))

    typer.echo("\n== available_uningested: public data, not yet loaded ==")
    report = available_uningested_report(frame)
    if report.height == 0:
        typer.echo("(none -- every remaining gap needs new observations)")
    else:
        typer.echo(str(report.head(top)))

    classes = equivalence_summary(frame)
    typer.echo("\n== largest equivalence classes ==")
    typer.echo(str(classes.head(top)))

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        written = {
            "blocking_summary.tsv": blocking,
            "available_uningested.tsv": report,
            "consequence_class.tsv": summarise_counts(
                frame, ["gene", "consequence", "class_current", "blocking_reason"]
            ),
            # Truncated on purpose: a whole gene has upwards of a hundred
            # thousand classes, and the tail of size-one classes carries no
            # information the per-variant map does not already have.
            f"equivalence_top{top}.tsv": classes.head(top),
        }
        for name, table in written.items():
            table.write_csv(out_dir / name, separator="\t")
        typer.echo(f"\nwrote {len(written)} TSV aggregations to {out_dir}")


@reference_app.command("from-mane")
def reference_from_mane(
    gene: Path = typer.Option(..., help="Path to the gene YAML."),
    gtf: Path = typer.Option(..., help="MANE genomic GTF (.gtf or .gtf.gz)."),
    fasta: Path = typer.Option(..., help="MANE RefSeq RNA FASTA (.fna or .fna.gz)."),
    out_dir: Path = typer.Option(Path("data/reference")),
    import_after: bool = typer.Option(
        True, help="Run 'reference import' on the extracted files immediately."
    ),
) -> None:
    """Extract one transcript's FASTA and exon table from a MANE Select release.

    The exon table comes out in **transcript order**, taken from the GTF's
    ``exon_number`` and then checked against the strand. Sorting by coordinate
    instead -- the obvious thing -- reverses a minus-strand gene into a
    transcript that is internally consistent and biologically wrong.

    Exon labels come from the gene config, not the GTF: BRCA1's GTF numbers its
    exons 1..23 while the clinical numbering runs 1,2,3,5,...,24.
    """
    from .genome.mane import extract_transcript

    config = load_gene_config(gene)
    try:
        extract = extract_transcript(config, gtf_path=gtf, fasta_path=fasta)
    except ReferenceUnavailable as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from None

    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_out = out_dir / f"{config.gene}_{extract.transcript_id}.fa"
    exons_out = out_dir / f"{config.gene}_exons.tsv"
    fasta_out.write_text(extract.fasta_text(), encoding="ascii")
    exons_out.write_text(extract.exon_table(), encoding="utf-8")

    typer.echo(
        f"{config.gene}: {len(extract.fasta):,} nt, {len(extract.exons)} exons on "
        f"{extract.chrom}{extract.strand}"
    )
    for warning in extract.warnings:
        typer.secho(f"  note: {warning}", fg=typer.colors.YELLOW)
    typer.echo(f"  wrote {fasta_out}")
    typer.echo(f"  wrote {exons_out}")

    if import_after:
        block = import_reference(
            config, sequence_path=fasta_out, exon_table_path=exons_out, config_path=gene
        )
        typer.echo(
            f"  imported: CDS c.1 at transcript position {block['cds_start_tx']}, "
            f"sha256 {block['sequence_sha256'][:16]}..."
        )


@reference_app.command("import")
def reference_import(
    gene: Path = typer.Option(..., help="Path to the gene YAML to populate."),
    sequence: Path = typer.Option(..., help="Spliced transcript FASTA (single record)."),
    exons: Path = typer.Option(..., help="TSV of label/start/end in transcript order."),
    write: bool = typer.Option(True, help="Write the derived block back into the YAML."),
) -> None:
    """Derive exon coordinates and CDS bounds from MANE Select resources."""
    config = load_gene_config(gene)
    try:
        block = import_reference(
            config,
            sequence_path=sequence,
            exon_table_path=exons,
            config_path=gene if write else None,
        )
    except ReferenceUnavailable as exc:
        typer.secho(f"import failed: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from None
    typer.echo(
        f"{config.gene}: CDS c.1 at transcript position {block['cds_start_tx']}, "
        f"{len(block['exons'])} exons, sha256 {block['sequence_sha256'][:16]}..."
    )
    if write:
        typer.echo(f"wrote the transcript block into {gene}")


@app.command("validate")
def validate_command(
    map_at_t: Path = typer.Option(
        ..., "--map-at-t", help="Gap map computed with the evidence state of date T."
    ),
    map_at_t_plus_n: Optional[Path] = typer.Option(
        None,
        "--map-at-t-plus-n",
        help=(
            "Gap map recomputed at T+n. With it, metric 2 comes from the map's own "
            "time series and no curated outcome table is needed."
        ),
    ),
    clinvar_at_t: Optional[Path] = typer.Option(
        None, help="ClinVar snapshot at T. The ground truth, which no criterion reads."
    ),
    clinvar_at_t_plus_n: Optional[Path] = typer.Option(None, help="ClinVar snapshot at T+n."),
    outcomes: Optional[Path] = typer.Option(
        None,
        help=(
            "Curated outcome table. Takes precedence over the snapshots: somebody "
            "who read the submission records knows what the sources cannot say."
        ),
    ),
    spec_dir: Path = typer.Option(Path("config/specs")),
    spec_name: Optional[str] = typer.Option(
        None, help="Specification id, defaulting to the one recorded in the map."
    ),
    reference_date: Optional[str] = typer.Option(
        None, help="ISO date T, for the temporal calibration metric."
    ),
    min_stars: int = typer.Option(
        1, help="Review-status floor for a ClinVar record to count as ground truth."
    ),
    show: int = typer.Option(10, help="How many example variants to list per section."),
    out_dir: Optional[Path] = typer.Option(
        None, help="Also write the study as TSV: metrics, predictions, misses."
    ),
) -> None:
    """Run the section 10 protocol.

    This is a validation study, not a test. It answers whether the map's claim is
    true: did the variants it called resolvable get resolved, for the reasons it
    predicted, in the direction it pointed, in the order it implied?

    Two modes. With a curated outcome table, it runs exactly as section 10
    describes. With two maps and two ClinVar snapshots it runs with no external
    table at all, because the ground truth is already there:
    ``clinvar.self.classification`` is published by the adapter and read by no
    criterion, so the oracle sits in the same snapshots the engine consumes and
    is firewalled from what is being measured.
    """
    from .adapters.clinvar import ClinVarSnapshot
    from .validation import read_outcomes, run_time_series_study
    from .validation import validate as run_validation

    # Streamed rather than read: a whole gene's rows cost about 9 GB as objects
    # and the time-series mode walks two maps. GapMapSource re-reads from disk
    # on each pass instead.
    rows = GapMapSource(map_at_t)
    first = next(iter(rows), None)
    if first is None:
        typer.secho(f"{map_at_t} has no rows", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2)
    when = date.fromisoformat(reference_date) if reference_date else None
    curated = read_outcomes(outcomes) if outcomes is not None else None

    # Narrowed directly rather than through a boolean: both paths below use
    # these two as non-optional, and an intermediate flag hides that from
    # the reader as much as from the type checker.
    if map_at_t_plus_n is None or clinvar_at_t is None:
        if curated is None:
            raise typer.BadParameter(
                "supply either --outcomes, or --map-at-t-plus-n with --clinvar-at-t "
                "and --clinvar-at-t-plus-n"
            )
        result = run_validation(rows, curated, reference_date=when)
        study = None
    else:
        if clinvar_at_t_plus_n is None:
            raise typer.BadParameter("--clinvar-at-t-plus-n is required in time-series mode")
        spec_id = spec_name or first.spec_version.split("@")[0]
        candidates = sorted(spec_dir.glob(f"{spec_id}*.yaml"))
        if not candidates:
            raise typer.BadParameter(f"no specification matching {spec_id!r} in {spec_dir}")
        study = run_time_series_study(
            rows,
            GapMapSource(map_at_t_plus_n),
            ClinVarSnapshot.read(clinvar_at_t),
            ClinVarSnapshot.read(clinvar_at_t_plus_n),
            load_spec(candidates[0]),
            transcript_id=first.transcript,
            reference_date=when,
            min_stars=min_stars,
            curated_outcomes=curated,
        )
        result = study.result

    typer.echo(study.summary() if study is not None else result.summary())

    if result.misses:
        typer.echo(f"\n{len(result.misses)} resolution(s) the map did not anticipate:")
        for line in result.misses[:show]:
            typer.echo(f"  {line}")
    if result.cause_confusion:
        typer.echo("\npredicted blocking_reason -> evidence that actually arrived:")
        for (predicted, observed), count in sorted(result.cause_confusion.items()):
            marker = "  " if predicted == observed else "! "
            typer.echo(f"  {marker}{predicted:28s} {observed:28s} {count:>6,}")
    if study is not None and study.unmoved:
        typer.echo(
            f"\n{len(study.unmoved)} resolved in ClinVar but not moved here "
            "(evidence never became public data, or a source is missing):"
        )
        for variant_id in study.unmoved[:show]:
            typer.echo(f"  {variant_id}")
    if study is not None and study.ahead_of_clinvar:
        typer.echo(
            f"\n{len(study.ahead_of_clinvar)} moved here and not yet in the archive "
            "-- these are the map's live predictions:"
        )
        for variant_id in study.ahead_of_clinvar[:show]:
            typer.echo(f"  {variant_id}")

    if out_dir is not None:
        _write_study(out_dir, result, study)
        typer.echo(f"\nwrote the study to {out_dir}")


def _write_study(out_dir: Path, result, study) -> None:
    """Persist a validation study as TSV.

    Metrics as name/value pairs rather than one wide row: the set of metrics has
    already grown once, and a long table survives that without every previous
    run's file having a different header from the next.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics: list[tuple[str, str]] = [
        ("considered", str(result.considered)),
        ("resolved", str(result.resolved)),
        ("resolvability_recall", _ratio(result.resolvability_recall)),
        ("cause_accuracy", _ratio(result.cause_accuracy)),
        ("cause_scored", str(result.cause_scored)),
        ("unobserved_cause", str(result.unobserved_cause)),
        ("direction_accuracy", _ratio(result.direction_accuracy)),
        ("direction_scored", str(result.direction_scored)),
        ("direction_unpredicted", str(result.direction_unpredicted)),
        (
            "temporal_correlation",
            "" if result.temporal_correlation is None else f"{result.temporal_correlation:.6f}",
        ),
    ]
    if study is not None:
        metrics += [
            ("anticipated", str(len(study.anticipated))),
            ("saw_evidence_only", str(len(study.saw_evidence_only))),
            ("unmoved", str(len(study.unmoved))),
            ("ahead_of_clinvar", str(len(study.ahead_of_clinvar))),
            ("anticipation_rate", _ratio(study.anticipation_rate)),
            ("evidence_visibility_rate", _ratio(study.evidence_visibility_rate)),
        ]
    _write_tsv(out_dir / "metrics.tsv", ("metric", "value"), metrics)

    _write_tsv(
        out_dir / "cause_confusion.tsv",
        ("predicted_blocking_reason", "observed_evidence", "variants"),
        [
            (predicted, observed, str(count))
            for (predicted, observed), count in sorted(result.cause_confusion.items())
        ],
    )
    _write_tsv(out_dir / "misses.tsv", ("miss",), [(line,) for line in result.misses])
    if study is not None:
        # The whole lists, not a sample: 'ahead_of_clinvar' is the set of live,
        # still-unfalsified predictions, and truncating it would throw away the
        # only output of this command that a future run can be scored against.
        for name, values in (
            ("predictions_ahead_of_clinvar.tsv", study.ahead_of_clinvar),
            ("unmoved.tsv", study.unmoved),
            ("saw_evidence_only.tsv", study.saw_evidence_only),
        ):
            _write_tsv(out_dir / name, ("variant_id",), [(v,) for v in values])


def _ratio(value: float | None) -> str:
    return "" if value is None else f"{value:.6f}"


def _write_tsv(path: Path, header: tuple[str, ...], rows: Sequence[tuple[str, ...]]) -> None:
    lines = ["\t".join(header)]
    lines += ["\t".join(field.replace("\t", " ") for field in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@data_app.command("check")
def data_check(
    gene: Path = typer.Option(..., help="Path to a gene YAML."),
    data_root: Path = typer.Option(Path(".")),
    spec_dir: Path = typer.Option(Path("config/specs")),
    frequency: Optional[Path] = typer.Option(None),
    predictor: Optional[Path] = typer.Option(None),
    splice: Optional[Path] = typer.Option(None),
    functional: Optional[Path] = typer.Option(None),
    clinvar: Optional[Path] = typer.Option(None),
    sample: int = typer.Option(500, help="Variants to sample when measuring coverage."),
) -> None:
    """Report what has been materialised, and how much of the gene it covers.

    Presence is not the useful question. A snapshot built against the wrong
    transcript version is empty, and a region slice with the wrong coordinates
    covers nothing -- both look exactly like "this variant has no data" once
    they reach the engine, which is why this runs before a map rather than
    during one.
    """
    from .adapters import ClinVarSnapshotAdapter
    from .datacheck import DataReport, check_reference, check_sources
    from .enumeration import enumerate_coding_snvs

    config = load_gene_config(gene)
    reference = check_reference(config, data_root=data_root)
    report = DataReport(gene=config.gene, reference=reference)

    if reference.present:
        transcript = build_transcript(config, data_root=data_root)
        variants = list(enumerate_coding_snvs(transcript))
        step = max(1, len(variants) // sample)
        sampled = variants[::step][:sample]
        report.sampled = len(sampled)

        adapters: list[Adapter] = []
        if frequency is not None:
            adapters.append(FrequencyAdapter(frequency, "check"))
        if predictor is not None:
            adapters.append(PredictorAdapter(predictor, "check"))
        if splice is not None:
            adapters.append(SpliceAdapter(splice, "check"))
        if functional is not None:
            adapters.append(FunctionalAdapter(functional, "check"))
        if clinvar is not None:
            adapters.append(ClinVarSnapshotAdapter.from_path(clinvar))
        report.sources = check_sources(adapters, transcript, sampled)

    typer.echo(report.summary())
    if not report.usable:
        raise typer.Exit(code=2)
    if report.empty_but_present:
        raise typer.Exit(code=1)


@clinvar_app.command("build")
def clinvar_build(
    source: Path = typer.Option(..., help="ClinVar variant_summary.txt(.gz) release."),
    out: Path = typer.Option(..., help="Normalised snapshot TSV to write."),
    transcript: Optional[str] = typer.Option(
        None, help="MANE transcript with its version, e.g. NM_007294.4."
    ),
    gene: Optional[Path] = typer.Option(
        None, help="Gene YAML to take the transcript from, instead of --transcript."
    ),
    assembly: str = typer.Option("GRCh38"),
) -> None:
    """Build a dated ClinVar snapshot for one transcript.

    The transcript is matched *with its version*: a classification made against
    NM_007294.3 is not evidence about a coordinate in NM_007294.4 unless
    somebody has checked that the two agree.
    """
    if (transcript is None) == (gene is None):
        raise typer.BadParameter("pass exactly one of --transcript or --gene")
    if gene is not None:
        transcript = load_gene_config(gene).transcript.id
    assert transcript is not None
    stats = build_snapshot(source, out, transcript=transcript, assembly=assembly)
    typer.echo(f"{out}: {stats.summary()}")
    if stats.kept == 0:
        typer.secho(
            "no records kept -- check the transcript version and the assembly",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)


@app.command("timeline")
def timeline_command(
    snapshots: list[str] = typer.Argument(
        ...,
        help=(
            "Two or more 'label=path.parquet' pairs, oldest first, e.g. "
            "2020-01-01=out/T2020/gap_map.parquet"
        ),
    ),
    out: Optional[Path] = typer.Option(None, help="Write the transitions as TSV."),
    show: int = typer.Option(10, help="Example transitions to print per interval."),
) -> None:
    """Diff maps computed at different dates and attribute what moved.

    The finding this exists for: variants that leave VUS because a *neighbour*
    was classified, with no new evidence about themselves. That is what makes
    the semi-intrinsic criteria worth separating from the intrinsic ones.
    """
    from .engine.timeline import compare_series

    parsed: list[tuple[str, Path]] = []
    for entry in snapshots:
        label, _, path = entry.partition("=")
        # The label is checked for looking like a path, not merely for
        # existing. The map writes its output to a Hive-style directory --
        # out/<date>/gene=BRCA1/gap_map.parquet -- so a path passed without a
        # label still contains an '=' and still splits, into a nonsense pair
        # whose only symptom is a FileNotFoundError about half a path.
        if not path or not label or "/" in label or os.sep in label:
            raise typer.BadParameter(
                f"expected 'label=path', got {entry!r}. Note that a map's own "
                "output path contains an '=' (gene=BRCA1), so passing one "
                "without a label prefix reads as a label."
            )
        parsed.append((label, Path(path)))
    if len(parsed) < 2:
        raise typer.BadParameter("a timeline needs at least two snapshots")

    # Streamed. compare_series reads each snapshot twice -- as the later map of
    # one interval and the earlier map of the next -- and GapMapSource makes
    # that a re-read rather than a second copy in memory.
    series = [(label, GapMapSource(path)) for label, path in parsed]
    diffs = compare_series(series)

    lines: list[str] = []
    for diff in diffs:
        typer.echo(diff.summary())
        for transition in diff.neighbour_driven_resolutions[:show]:
            typer.echo(
                f"    {transition.hgvs_c} {transition.hgvs_p or '':16s} "
                f"{transition.class_before.value} -> {transition.class_after.value} "
                f"via {'+'.join(code for code, _ in transition.criteria_gained)}"
            )
        typer.echo("")
        for transition in diff.transitions:
            lines.append(
                "\t".join(
                    [
                        diff.label_before,
                        diff.label_after,
                        transition.variant_id,
                        transition.hgvs_p or "",
                        transition.consequence,
                        transition.class_before.value,
                        transition.class_after.value,
                        str(transition.points_before),
                        str(transition.points_after),
                        transition.blocking_before.value,
                        transition.blocking_after.value,
                        ";".join(f"{c}:{s}" for c, s in transition.criteria_gained),
                        ";".join(f"{c}:{s}" for c, s in transition.criteria_lost),
                        transition.cause.value,
                    ]
                )
            )

    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        header = "\t".join(  # noqa: FLY002 -- a column list diffs one line per column
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
        out.write_text("\n".join([header, *lines]) + "\n", encoding="utf-8")
        typer.echo(f"wrote {len(lines):,} transitions to {out}")


@app.command("selftest")
def selftest(
    spec_dir: Path = typer.Option(Path("config/specs")),
    spec_name: str = typer.Option("toy_v0.1.0"),
    out: Optional[Path] = typer.Option(None, help="Write the synthetic map here."),
    clinvar: Optional[Path] = typer.Option(
        None,
        help=(
            "Optional ClinVar snapshot. Running selftest twice with different "
            "snapshots and diffing the results with 'timeline' demonstrates the "
            "section 4 observation without any reference data."
        ),
    ),
    clinvar_date: Optional[str] = typer.Option(None, help="ISO date of that snapshot."),
) -> None:
    """Run the whole pipeline on a synthetic gene, with no reference data.

    This is the domain-isolation check of spec section 11 in executable form: if
    it works on a gene that does not exist, no BRCA constant is hiding in the
    engine.
    """
    from .testing import build_demo_runner

    snapshot_date = date.fromisoformat(clinvar_date) if clinvar_date else None
    runner, variants = build_demo_runner(
        spec_dir / f"{spec_name}.yaml",
        clinvar_snapshot_path=clinvar,
        clinvar_snapshot_date=snapshot_date,
    )
    rows = [result.row for result in runner.run(variants)]
    typer.echo(f"{runner.gene.gene}: {len(rows):,} rows under {runner.spec.spec_version}")
    typer.echo(f"  {runner.cache.summary()}")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.blocking_reason.value] = counts.get(row.blocking_reason.value, 0) + 1
    for reason, count in sorted(counts.items()):
        typer.echo(f"  {reason:28s} {count:>7,}")
    if out is not None:
        write_parquet(rows, out)
        typer.echo(f"wrote {out}")


if __name__ == "__main__":  # pragma: no cover
    app()
