"""Layer 3: class-level evaluation must be a saving, never a shortcut.

The cache assumes a specification's declared field footprint covers everything
an evaluation reads. If that assumption were ever wrong, two variants that
genuinely differ would silently share a verdict -- the worst possible failure
for this system, because it would be invisible in the output.

So the assumption is not argued, it is tested three ways:

* every row produced with reuse on is identical to the row produced with it off,
  across every enumeration class;
* the paths an evaluation actually reads, recorded by the context's audit mode,
  are contained in the declared footprint;
* perturbing any single footprint path changes the signature.
"""

from __future__ import annotations

from datetime import datetime

import pytest

from vus_foresight.engine.context import EvidenceContext
from vus_foresight.engine.pipeline import MapRunner, default_registry
from vus_foresight.engine.signature import canonical, evidence_signature
from vus_foresight.enumeration import EnumerationClass, enumerate_all
from vus_foresight.variant import Consequence

CLASSES = (
    EnumerationClass.CODING_SNV,
    EnumerationClass.INTRONIC_SNV,
    EnumerationClass.INTRACODON_MNV,
    EnumerationClass.FRAMESHIFT_CLASS,
    EnumerationClass.INFRAME_DELETION,
)


def _runner(gene, config, spec, **kwargs):
    return MapRunner(
        transcript=gene.transcript,
        gene=config,
        spec=spec,
        adapters=default_registry(config),
        computed_at=datetime(1970, 1, 1),
        **kwargs,
    )


def _variants(gene):
    return list(enumerate_all(gene.transcript, classes=CLASSES, flanks=gene.flanks))


def _fully_evidenced_adapters(tmp_path, gene, variants):
    """Snapshots that leave no criterion unevaluable.

    Needed by the "footprint is not padded" test: a path that only a rendered
    evidence template reads is never touched while its criterion sits at
    NOT_EVALUABLE, so a run with no data cannot distinguish a padded footprint
    from a correct one.
    """
    from datetime import date

    from vus_foresight.adapters import (
        ClinVarSnapshotAdapter,
        FrequencyAdapter,
        FunctionalAdapter,
        PredictorAdapter,
        SpliceAdapter,
    )

    coding = [v for v in variants if v.grch38_pos and v.codon_index]
    proteins = sorted({v.hgvs_p for v in variants if v.hgvs_p})

    frequency = tmp_path / "frequency.tsv"
    frequency.write_text(
        "grch38_pos\tgnomad.faf95_popmax\n"
        + "".join(
            f"{v.grch38_pos}\t{0.002 if i % 3 == 0 else 0.0005 if i % 3 == 1 else 0.0}\n"
            for i, v in enumerate(coding)
        ),
        encoding="utf-8",
    )
    predictor = tmp_path / "predictor.tsv"
    predictor.write_text(
        "grch38_pos\tbayesdel\n"
        + "".join(f"{v.grch38_pos}\t{0.62 if i % 2 else 0.05}\n" for i, v in enumerate(coding)),
        encoding="utf-8",
    )
    splice = tmp_path / "splice.tsv"
    splice.write_text(
        "grch38_pos\tds_max\n"
        + "".join(
            f"{v.grch38_pos}\t{0.9 if i % 5 == 0 else 0.01}\n"
            for i, v in enumerate(v for v in variants if v.grch38_pos)
        ),
        encoding="utf-8",
    )
    functional = tmp_path / "functional.tsv"
    functional.write_text(
        "hgvs_p\tclassification\tscore\tdataset\n"
        + "".join(
            f"{p}\t{'abnormal' if i % 2 else 'normal'}\t{-2.1 if i % 2 else 0.1}\tTOY-SGE\n"
            for i, p in enumerate(proteins)
        ),
        encoding="utf-8",
    )
    # Every ClinVar record here is a *neighbour*: a different nucleotide change
    # reaching the same protein change (PS1), or a different protein change at
    # the same codon (PM5). A record matching the variant itself would be
    # excluded by the adapter, so a fixture built from the variants' own
    # descriptions would silently test nothing.
    clinvar = tmp_path / "clinvar.tsv"
    rows = ["hgvs_c\thgvs_p\tcodon\tclassification\tstars\tlast_evaluated"]
    for i, v in enumerate(coding):
        if v.hgvs_p and i % 4 == 0:
            rows.append(f"c.{900000 + i}A>G\t{v.hgvs_p}\t{v.codon_index}\tpathogenic\t3\t")
        if i % 3 == 0:
            rows.append(
                f"c.{800000 + i}A>G\tp.Ala{v.codon_index}Trp\t{v.codon_index}"
                f"\tlikely_pathogenic\t2\t"
            )
    clinvar.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return [
        FrequencyAdapter(frequency, "test", default_payload={"gnomad": {"faf95_popmax": 0.0}}),
        PredictorAdapter(predictor, "test"),
        SpliceAdapter(splice, "test", default_payload={"ds_max": 0.0}),
        FunctionalAdapter(
            functional, "test", assayed_regions=((1, gene.transcript.protein_length, "TOY-SGE"),)
        ),
        ClinVarSnapshotAdapter.from_path(clinvar, snapshot_date=date(2026, 1, 1)),
    ]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_reuse_produces_exactly_the_reference_rows(both_strands, strand, toy_spec):
    """The uncached path is the reference implementation, not a fallback."""
    gene = both_strands[strand]
    config = __import__(
        "vus_foresight.testing", fromlist=["synthetic_gene_config"]
    ).synthetic_gene_config(gene)
    variants = _variants(gene)

    reference = [
        r.row for r in _runner(gene, config, toy_spec, reuse_by_signature=False).run(variants)
    ]
    cached_runner = _runner(gene, config, toy_spec, reuse_by_signature=True)
    cached = [r.row for r in cached_runner.run(variants)]

    assert len(reference) == len(cached) == len(variants)
    assert reference == cached
    assert cached_runner.cache.hits > 0, "the fixture should exercise the cache"


def test_the_declared_footprint_covers_every_path_an_evaluation_reads(
    minus_gene, minus_config, toy_spec
):
    """Proof by audit rather than by reading the call sites.

    A future code path that reaches for an undeclared field fails here, at the
    moment it is added, instead of quietly merging two variants that differ.
    """
    runner = _runner(
        minus_gene, minus_config, toy_spec, audit_context_reads=True, reuse_by_signature=False
    )
    list(runner.run(_variants(minus_gene)))

    declared = set(toy_spec.field_footprint)
    undeclared = runner.audited_paths - declared
    assert not undeclared, (
        f"the evaluation reads {sorted(undeclared)}, which the specification's "
        "field footprint does not declare. Class-level reuse would merge variants "
        "that differ on those paths."
    )
    assert runner.audited_paths, "the audit recorded nothing, so it proves nothing"


def test_the_footprint_is_not_padded_with_paths_nobody_reads(
    tmp_path, minus_gene, minus_config, toy_spec
):
    """The converse guard: a footprint of everything would be sound and useless.

    An entry nobody reads still enters the signature, so it can only split
    classes that should have merged.

    This has to run with every source populated. Paths that only a rendered
    evidence template touches -- ``functional.dataset`` is one -- are never read
    while their criterion sits at NOT_EVALUABLE, so a run with no data cannot
    tell a padded footprint from a correct one. Making every criterion evaluable
    also proves each one is reachable at all.
    """
    from vus_foresight.adapters import AdapterRegistry
    from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter

    variants = _variants(minus_gene)
    registry = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(minus_config)])
    for adapter in _fully_evidenced_adapters(tmp_path, minus_gene, variants):
        registry.add(adapter)

    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=registry,
        computed_at=datetime(1970, 1, 1),
        audit_context_reads=True,
        reuse_by_signature=False,
    )
    rows = [r.row for r in runner.run(variants)]

    unread = set(toy_spec.field_footprint) - runner.audited_paths
    assert not unread, f"declared but never read: {sorted(unread)}"

    # Every criterion with a rule should be reachable under some evidence state.
    applied = {c.code for row in rows for c in row.criteria_applied}
    rule_bearing = {c.code for c in toy_spec.criteria if c.rule is not None}
    assert rule_bearing - applied == set(), (
        f"criteria that never fired even with every source populated: "
        f"{sorted(rule_bearing - applied)}"
    )


def test_reuse_matches_the_reference_with_every_source_populated(
    tmp_path, minus_gene, minus_config, toy_spec
):
    """The equality test again, but where the evidence actually varies per variant.

    With no adapters loaded almost everything shares one profile, so a cache bug
    could hide behind the uniformity. Here the frequency, predictor, splice,
    assay and ClinVar values all differ across variants.
    """
    from vus_foresight.adapters import AdapterRegistry
    from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter

    variants = _variants(minus_gene)

    def build(reuse):
        registry = AdapterRegistry(
            [VariantAdapter(), TranscriptAdapter(), RegionAdapter(minus_config)]
        )
        for adapter in _fully_evidenced_adapters(tmp_path, minus_gene, variants):
            registry.add(adapter)
        return MapRunner(
            transcript=minus_gene.transcript,
            gene=minus_config,
            spec=toy_spec,
            adapters=registry,
            computed_at=datetime(1970, 1, 1),
            reuse_by_signature=reuse,
        )

    reference = [r.row for r in build(False).run(variants)]
    cached_runner = build(True)
    cached = [r.row for r in cached_runner.run(variants)]
    assert reference == cached
    assert cached_runner.cache.size > 1, "the evidence should not be uniform here"


def test_changing_any_footprint_path_changes_the_signature(minus_gene, minus_config, toy_spec):
    runner = _runner(minus_gene, minus_config, toy_spec)
    variant = next(v for v in _variants(minus_gene) if v.consequence is Consequence.MISSENSE)
    context = runner.build_context(variant)
    footprint = toy_spec.field_footprint
    baseline = evidence_signature(variant, context, footprint)

    for path in footprint:
        if path == "variant.consequence_terms":
            continue
        namespace, *rest = path.split(".")
        perturbed = EvidenceContext(
            data={key: dict(value) for key, value in context.data.items()},
            sources=dict(context.sources),
        )
        node = perturbed.data.setdefault(namespace, {})
        for segment in rest[:-1]:
            node = node.setdefault(segment, {})
            if not isinstance(node, dict):  # pragma: no cover - defensive
                break
        node[rest[-1]] = "\x01perturbed"
        assert evidence_signature(variant, perturbed, footprint) != baseline, (
            f"perturbing {path} did not change the signature"
        )


def test_consequence_terms_are_part_of_the_signature(minus_gene, minus_config, toy_spec):
    """The applies_to gate reads them from the variant, not from the context."""
    runner = _runner(minus_gene, minus_config, toy_spec)
    variants = _variants(minus_gene)
    missense = next(v for v in variants if v.consequence is Consequence.MISSENSE)
    synonymous = next(v for v in variants if v.consequence is Consequence.SYNONYMOUS)
    footprint = toy_spec.field_footprint
    assert evidence_signature(
        missense, runner.build_context(missense), footprint
    ) != evidence_signature(synonymous, runner.build_context(synonymous), footprint)


def test_missense_shares_computation_but_never_shares_a_class(minus_gene, minus_config, toy_spec):
    """Section 5 forbids collapsing missense in the report, not in the arithmetic.

    With no predictor or assay loaded, every missense variant has the same
    evidence profile and is evaluated once -- but each still gets its own
    equivalence class, because the moment a score lands they diverge.
    """
    runner = _runner(minus_gene, minus_config, toy_spec)
    rows = [r.row for r in runner.run(_variants(minus_gene))]
    missense = [row for row in rows if row.consequence is Consequence.MISSENSE]
    assert len(missense) > 50
    assert len({row.equivalence_class_id for row in missense}) == len(missense)
    assert runner.cache.size < len(rows) / 10


def test_truncating_variants_are_not_over_merged(minus_gene, minus_config, toy_spec):
    """A real difference must survive the signature.

    PVS1's evidence string names the PTC codon, so two nonsense variants
    terminating at different codons read different values on a footprint path
    and are correctly kept apart.
    """
    runner = _runner(minus_gene, minus_config, toy_spec)
    rows = [r.row for r in runner.run(_variants(minus_gene))]
    nonsense = [row for row in rows if row.consequence is Consequence.NONSENSE]
    assert nonsense
    evidence = {c.evidence for row in nonsense for c in row.criteria_applied if c.code == "PVS1"}
    assert len(evidence) > 1, "distinct PTC positions collapsed into one evidence string"


def test_canonical_is_order_independent_for_mappings_only():
    assert canonical({"b": 1, "a": 2}) == canonical({"a": 2, "b": 1})
    # Lists keep their order: a different order means a different adapter output.
    assert canonical(["a", "b"]) != canonical(["b", "a"])


def test_every_shipped_specification_declares_a_usable_footprint():
    """Structural, so it holds for a specification whose thresholds are uncurated.

    A spec whose footprint came out empty would make every variant share one
    signature -- the failure mode that turns the cache from a saving into a bug.
    """
    from pathlib import Path

    from vus_foresight.engine.spec import load_spec

    specs = sorted((Path(__file__).resolve().parents[2] / "config" / "specs").glob("*.yaml"))
    assert specs
    for path in specs:
        spec = load_spec(path)
        footprint = spec.field_footprint
        assert footprint, f"{path.name} declares no readable fields"
        assert len(set(footprint)) == len(footprint)
        assert footprint == tuple(sorted(footprint))
        assert "pvs1.strength" in footprint
        # Every path a criterion can read must be reachable from the footprint.
        for criterion in spec.criteria:
            assert set(criterion.iter_fields()) <= set(footprint), criterion.code


def test_cache_counters_describe_the_run(minus_gene, minus_config, toy_spec):
    runner = _runner(minus_gene, minus_config, toy_spec)
    variants = _variants(minus_gene)
    list(runner.run(variants))
    assert runner.cache.lookups == len(variants)
    assert runner.cache.hits + runner.cache.misses == len(variants)
    assert runner.cache.size == runner.cache.misses
    assert 0.0 < runner.cache.hit_rate < 1.0
    assert "distinct evidence profiles" in runner.cache.summary()
