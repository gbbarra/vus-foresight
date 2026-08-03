"""Shared fixtures.

The suite runs entirely on synthetic genes. Tests that need a real MANE Select
transcript are marked ``requires_reference`` and skip cleanly when the resources
have not been imported -- see ``docs/reference-data.md``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vus_foresight.engine.spec import load_spec
from vus_foresight.genome.reference import ReferenceUnavailable, build_transcript, load_gene_config
from vus_foresight.testing import build_synthetic_gene, synthetic_gene_config

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = REPO_ROOT / "config"
DATA_ROOT = REPO_ROOT / "data"


@pytest.fixture(scope="session")
def toy_spec():
    return load_spec(CONFIG_DIR / "specs" / "toy_v0.1.0.yaml")


@pytest.fixture(scope="session")
def plus_gene():
    return build_synthetic_gene(strand="+", gene="TOYP", transcript_id="NM_999001.1")


@pytest.fixture(scope="session")
def minus_gene():
    return build_synthetic_gene(strand="-", gene="TOYM", transcript_id="NM_999002.1")


@pytest.fixture(scope="session")
def both_strands(plus_gene, minus_gene):
    return {"+": plus_gene, "-": minus_gene}


@pytest.fixture(scope="session")
def plus_config(plus_gene):
    return synthetic_gene_config(plus_gene)


@pytest.fixture(scope="session")
def minus_config(minus_gene):
    return synthetic_gene_config(minus_gene)


def _real_transcript(name: str):
    config_path = CONFIG_DIR / "genes" / f"{name}.yaml"
    if not config_path.exists():
        pytest.skip(f"{name}.yaml is not present")
    config = load_gene_config(config_path)
    try:
        return config, build_transcript(config, data_root=DATA_ROOT)
    except ReferenceUnavailable as exc:
        pytest.skip(f"{name} reference not imported: {exc}")


@pytest.fixture(scope="session")
def brca1():
    return _real_transcript("BRCA1")


@pytest.fixture(scope="session")
def brca2():
    return _real_transcript("BRCA2")
