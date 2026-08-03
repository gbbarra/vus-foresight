"""Command line interface.

``vus-foresight --help`` for the list. The commands that matter:

    enumerate   how many variants exist, by class -- the layer-1 invariants
    map         the gap map itself, to Parquet
    report      the aggregations: blocking reasons, available-uningested
    selftest    the whole pipeline on a synthetic gene, no reference needed
"""

from __future__ import annotations

from datetime import date, datetime, timezone
from pathlib import Path
from typing import Optional

import typer

from .adapters import (
    ClinVarAdapter,
    FrequencyAdapter,
    FunctionalAdapter,
    PredictorAdapter,
    SpliceAdapter,
)
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
    available_uningested_report,
    blocking_summary,
    equivalence_summary,
    read_parquet,
    write_parquet,
)

app = typer.Typer(
    add_completion=False,
    help="Evidence gap map for clinical genomics. No patient data, no LLM, no network.",
)
reference_app = typer.Typer(help="Manage MANE Select reference resources.")
app.add_typer(reference_app, name="reference")

#: Determinism: the default stamp is a fixed epoch, not the wall clock. A run
#: that wants a real timestamp must say so, because two runs that differ only in
#: their timestamp are not reproducible and the CI check would catch it as noise.
FIXED_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc).replace(tzinfo=None)


def _load(gene_path: Path, spec_dir: Path, data_root: Path):
    """Load a gene, its specification and its reference, or fail readably.

    A missing reference is an expected state, not a crash: the resources are not
    vendored on purpose (see docs/reference-data.md), so the message has to say
    what to run next rather than print a stack trace at a curator.
    """
    config = load_gene_config(gene_path)
    spec = load_spec(spec_dir / f"{config.spec}.yaml")
    try:
        transcript = build_transcript(config, data_root=data_root)
        flanks = load_flanks(config, data_root=data_root)
    except ReferenceUnavailable as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from None
    return config, spec, transcript, flanks


@app.command("enumerate")
def enumerate_command(
    gene: Path = typer.Option(..., help="Path to a gene YAML."),
    data_root: Path = typer.Option(Path("data"), help="Root for reference resources."),
    spec_dir: Path = typer.Option(Path("config/specs")),
    out: Optional[Path] = typer.Option(None, help="Write the enumerated variants as TSV."),
) -> None:
    """Count and optionally dump every possible variant for a gene."""
    config, _spec, transcript, flanks = _load(gene, spec_dir, data_root)

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
    data_root: Path = typer.Option(Path("data")),
    spec_dir: Path = typer.Option(Path("config/specs")),
    out_dir: Path = typer.Option(Path("out"), help="Directory for the Parquet output."),
    frequency: Optional[Path] = typer.Option(None, help="gnomAD/ABraOM snapshot TSV."),
    predictor: Optional[Path] = typer.Option(None, help="dbNSFP snapshot TSV."),
    splice: Optional[Path] = typer.Option(None, help="SpliceAI snapshot TSV."),
    functional: Optional[Path] = typer.Option(None, help="MAVE/SGE snapshot TSV."),
    clinvar: Optional[Path] = typer.Option(None, help="Dated ClinVar snapshot TSV."),
    clinvar_date: Optional[str] = typer.Option(None, help="ISO date of the ClinVar snapshot."),
    gnomad_version: str = typer.Option("v4", help="Recorded in every row's provenance."),
    computed_at: Optional[str] = typer.Option(
        None, help="ISO timestamp stamped on every row. Defaults to a fixed epoch."
    ),
    tiers: str = typer.Option("1,2", help="Enumeration tiers to include."),
    allow_unverified: bool = typer.Option(
        False, help="Write a map from a specification whose thresholds are uncurated."
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

    extra = []
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
        extra.append(ClinVarAdapter(clinvar, clinvar_date or "unknown", snapshot_date=clinvar_date))

    stamp = datetime.fromisoformat(computed_at) if computed_at else FIXED_EPOCH
    runner = MapRunner(
        transcript=transcript,
        gene=config,
        spec=spec,
        adapters=default_registry(config, extra=extra),
        computed_at=stamp,
        clinvar_snapshot=date.fromisoformat(clinvar_date) if clinvar_date else None,
        gnomad_version=gnomad_version if frequency is not None else None,
    )

    variants = list(
        enumerate_all(
            transcript,
            classes=tuple(classes),
            flanks=flanks,
            flank_bp=config.intronic_flank_bp,
        )
    )
    rows = [result.row for result in runner.run(variants)]
    out_dir.mkdir(parents=True, exist_ok=True)
    sequence_path = out_dir / f"gene={config.gene}" / "gap_map.parquet"
    write_parquet(rows, sequence_path)
    typer.echo(f"wrote {len(rows):,} sequence-level rows to {sequence_path}")

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
) -> None:
    """Print the aggregations the map exists to support."""
    frame = read_parquet(parquet)
    typer.echo(f"{frame.height:,} rows\n")

    typer.echo("== blocking reason ==")
    typer.echo(str(blocking_summary(frame)))

    typer.echo("\n== available_uningested: public data, not yet loaded ==")
    report = available_uningested_report(frame)
    if report.height == 0:
        typer.echo("(none -- every remaining gap needs new observations)")
    else:
        typer.echo(str(report.head(top)))

    typer.echo("\n== largest equivalence classes ==")
    typer.echo(str(equivalence_summary(frame).head(top)))


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
    parquet: Path = typer.Argument(..., help="A gap map computed with the evidence of date T."),
    outcomes: Path = typer.Argument(..., help="TSV of what ClinVar recorded by T+n."),
    reference_date: Optional[str] = typer.Option(
        None, help="ISO date T, for the temporal calibration metric."
    ),
    show_misses: int = typer.Option(10, help="How many unexplained resolutions to list."),
) -> None:
    """Run the section 10 protocol against vus-hindsight outcomes.

    This is a validation study, not a test. It answers whether the map's claim is
    true: did the variants it called resolvable get resolved, for the reasons it
    predicted, in the direction it pointed, in the order it implied?
    """
    from .validation import read_outcomes, validate as run_validation

    frame = read_parquet(parquet)
    rows = _frame_to_rows(frame)
    result = run_validation(
        rows,
        read_outcomes(outcomes),
        reference_date=date.fromisoformat(reference_date) if reference_date else None,
    )
    typer.echo(result.summary())
    if result.misses:
        typer.echo(f"\n{len(result.misses)} resolution(s) the map did not anticipate:")
        for line in result.misses[:show_misses]:
            typer.echo(f"  {line}")
    if result.cause_confusion:
        typer.echo("\npredicted blocking_reason -> observed evidence type:")
        for (predicted, observed), count in sorted(result.cause_confusion.items()):
            marker = "  " if predicted == observed else "! "
            typer.echo(f"  {marker}{predicted:28s} {observed:28s} {count:>6,}")


def _frame_to_rows(frame):
    """Rehydrate the subset of each row that the validation protocol reads."""
    from .acmg import ACMGClass, BlockingReason, FeasibilityTag
    from .gapmap import EvidenceSet, GapMapRow
    from .variant import Consequence, VariantKind

    rows = []
    for record in frame.iter_rows(named=True):
        sets = tuple(
            EvidenceSet(
                target=ACMGClass(s["target"]),
                requirements=(),
                total_points=s["total_points"],
                feasibility=FeasibilityTag(s["feasibility"]),
                acquisition_cost=s["acquisition_cost"],
            )
            for s in record["minimum_sufficient_sets"]
        )
        rows.append(
            GapMapRow(
                gene=record["gene"],
                transcript=record["transcript"],
                hgvs_c=record["hgvs_c"],
                hgvs_p=record["hgvs_p"],
                grch38_pos=record["grch38_pos"],
                consequence=Consequence(record["consequence"]),
                variant_kind=VariantKind(record["variant_kind"]),
                equivalence_class_id=record["equivalence_class_id"],
                mutational_distance=record["mutational_distance"],
                criteria_applied=(),
                criteria_evaluated_not_applied=(),
                points_current=record["points_current"],
                class_current=ACMGClass(record["class_current"]),
                points_ceiling_intrinsic=record["points_ceiling_intrinsic"],
                class_ceiling_intrinsic=ACMGClass(record["class_ceiling_intrinsic"]),
                gap_to_LP=record["gap_to_LP"],
                gap_to_LB=record["gap_to_LB"],
                minimum_sufficient_sets=sets,
                blocking_reason=BlockingReason(record["blocking_reason"]),
                spec_version=record["spec_version"],
                clinvar_snapshot=record["clinvar_snapshot"],
                gnomad_version=record["gnomad_version"],
                computed_at=record["computed_at"],
            )
        )
    return rows


@app.command("selftest")
def selftest(
    spec_dir: Path = typer.Option(Path("config/specs")),
    spec_name: str = typer.Option("toy_v0.1.0"),
    out: Optional[Path] = typer.Option(None, help="Write the synthetic map here."),
) -> None:
    """Run the whole pipeline on a synthetic gene, with no reference data.

    This is the domain-isolation check of spec section 11 in executable form: if
    it works on a gene that does not exist, no BRCA constant is hiding in the
    engine.
    """
    from .testing import build_demo_runner

    runner, variants = build_demo_runner(spec_dir / f"{spec_name}.yaml")
    rows = [result.row for result in runner.run(variants)]
    typer.echo(f"{runner.gene.gene}: {len(rows):,} rows under {runner.spec.spec_version}")
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
