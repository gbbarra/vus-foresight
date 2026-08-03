"""Layer 3: no domain constant may leak into the engine.

Spec section 9 and principle 5: adding ATM or PALB2 must be adding a YAML and a
transcript, not writing code. The suite already runs on a synthetic gene under a
toy specification, which is most of the argument. This file makes the remainder
mechanical: the *executable* code of the engine, enumeration, annotation and
genome packages is parsed, its docstrings removed, and the result checked for
any gene, transcript, expert-panel or dataset name.

Docstrings are stripped rather than ignored because the point is that BRCA may
be *discussed* anywhere and *depended on* nowhere.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[2] / "src" / "vus_foresight"

#: Packages that must be free of domain constants. ``testing.py`` and ``cli.py``
#: are excluded: the first builds fixtures, the second names config paths.
ENGINE_PACKAGES = ("engine", "enumeration", "annotate", "genome", "adapters")

FORBIDDEN = (
    "BRCA1",
    "BRCA2",
    "ENIGMA",
    "NM_007294",
    "NM_000059",
    "ATM",
    "PALB2",
    "MaveDB",
    "gnomAD",
    "ABraOM",
)


def _strip_docstrings(tree: ast.AST) -> ast.AST:
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        first = body[0]
        if (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
        ):
            node.body = body[1:] or [ast.Pass()]
    return tree


def _source_files() -> list[Path]:
    files: list[Path] = []
    for package in ENGINE_PACKAGES:
        files.extend(sorted((PACKAGE / package).rglob("*.py")))
    files.append(PACKAGE / "acmg.py")
    files.append(PACKAGE / "variant.py")
    files.append(PACKAGE / "gapmap.py")
    return files


@pytest.mark.parametrize("path", _source_files(), ids=lambda p: p.name)
def test_engine_code_contains_no_domain_constant(path: Path):
    tree = _strip_docstrings(ast.parse(path.read_text(encoding="utf-8")))
    code = ast.unparse(tree)
    offenders = [token for token in FORBIDDEN if token in code]
    assert not offenders, (
        f"{path.relative_to(PACKAGE)} depends on the domain constant(s) {offenders}. "
        "Move it to config/genes or config/specs."
    )


def test_the_whole_pipeline_runs_on_a_gene_that_does_not_exist(
    minus_gene, minus_config, toy_spec
):
    """The executable form of the isolation guarantee."""
    from datetime import datetime

    from vus_foresight.engine.pipeline import MapRunner, default_registry
    from vus_foresight.enumeration import EnumerationClass, enumerate_all

    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=default_registry(minus_config),
        computed_at=datetime(1970, 1, 1),
    )
    variants = list(
        enumerate_all(
            minus_gene.transcript,
            classes=(
                EnumerationClass.CODING_SNV,
                EnumerationClass.INTRONIC_SNV,
                EnumerationClass.FRAMESHIFT_CLASS,
                EnumerationClass.INFRAME_DELETION,
            ),
            flanks=minus_gene.flanks,
        )
    )
    rows = [result.row for result in runner.run(variants)]
    assert rows
    assert {row.gene for row in rows} == {minus_config.gene}
    assert {row.spec_version for row in rows} == {toy_spec.spec_version}


def test_a_second_synthetic_gene_needs_no_code_change(toy_spec):
    """Adding a gene is adding a config. Demonstrated, not asserted in prose."""
    from datetime import datetime

    from vus_foresight.engine.pipeline import MapRunner, default_registry
    from vus_foresight.enumeration import enumerate_coding_snvs
    from vus_foresight.testing import build_synthetic_gene, synthetic_gene_config

    other = build_synthetic_gene(
        gene="TOY2",
        transcript_id="NM_999003.1",
        chrom="chr98",
        strand="+",
        seed=7,
        exon_labels=("a", "b", "c", "d", "e"),
    )
    config = synthetic_gene_config(other, lof_mechanism="not_established")
    runner = MapRunner(
        transcript=other.transcript,
        gene=config,
        spec=toy_spec,
        adapters=default_registry(config),
        computed_at=datetime(1970, 1, 1),
    )
    rows = [r.row for r in runner.run(list(enumerate_coding_snvs(other.transcript))[:50])]
    assert rows
    assert {row.gene for row in rows} == {"TOY2"}
