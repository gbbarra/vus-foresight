"""Re-export the packaged synthetic-gene fixtures.

The builders live in :mod:`vus_foresight.testing` because the CLI's ``selftest``
command uses them too, and because the domain-isolation guarantee they
demonstrate is a property of the library rather than of the test suite.
"""

from vus_foresight.testing import (
    TOY_EXON_LENGTHS,
    TOY_INTRON_LENGTHS,
    TOY_PROTEIN_LENGTH,
    TOY_UTR3,
    TOY_UTR5,
    SyntheticGene,
    build_demo_runner,
    build_synthetic_gene,
    demo_adapters,
    synthetic_gene_config,
)

__all__ = [
    "TOY_EXON_LENGTHS",
    "TOY_INTRON_LENGTHS",
    "TOY_PROTEIN_LENGTH",
    "TOY_UTR3",
    "TOY_UTR5",
    "SyntheticGene",
    "build_demo_runner",
    "build_synthetic_gene",
    "demo_adapters",
    "synthetic_gene_config",
]
