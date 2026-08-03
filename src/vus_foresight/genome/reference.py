"""Gene configuration and reference resource loading (spec section 9).

The engine receives ``(variant, gene_config, spec)`` and contains no domain
constant of its own. This module is where the domain constants live, and they
live in YAML.

Reference *sequence* is deliberately not vendored. Exon coordinates and the
transcript sequence must come from MANE Select itself, so a gene YAML declares
where to find them and what they must hash to; :func:`build_transcript` refuses
to proceed on a mismatch. Making up coordinates would produce HGVS that looks
correct and is wrong -- exactly the failure mode spec section 11 calls the
riskiest in the system.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .transcript import Exon, FlankSequences, Transcript

__all__ = [
    "ReferenceUnavailable",
    "SequenceResource",
    "TranscriptConfig",
    "FunctionalRegion",
    "GeneConfig",
    "load_gene_config",
    "build_transcript",
    "load_flanks",
    "read_fasta",
]


class ReferenceUnavailable(RuntimeError):
    """Raised when a declared reference resource is absent or fails its digest.

    Callers are expected to let this propagate. Silently substituting a
    placeholder sequence would turn a missing-data error into a wrong-answer
    error, which is the one trade this project never makes.
    """


class SequenceResource(BaseModel):
    """A sequence file plus the digest it must have."""

    model_config = ConfigDict(frozen=True)

    path: str
    sha256: str | None = None
    expected_length: int | None = None
    #: Where a curator can obtain the file. Never fetched at engine runtime.
    provenance: str | None = None

    def resolve(self, root: Path) -> Path:
        candidate = Path(self.path)
        return candidate if candidate.is_absolute() else root / candidate

    def read(self, root: Path) -> str:
        target = self.resolve(root)
        if not target.exists():
            raise ReferenceUnavailable(
                f"reference sequence {target} is missing. "
                + (f"Obtain it from: {self.provenance}. " if self.provenance else "")
                + "See docs/reference-data.md."
            )
        seq = read_fasta(target)
        if self.expected_length is not None and len(seq) != self.expected_length:
            raise ReferenceUnavailable(
                f"{target}: expected {self.expected_length} nt, found {len(seq)}"
            )
        if self.sha256 is not None:
            digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
            if digest != self.sha256:
                raise ReferenceUnavailable(
                    f"{target}: sha256 mismatch. Declared {self.sha256}, found {digest}. "
                    "The reference changed under a pinned config; refusing to run."
                )
        return seq


class TranscriptConfig(BaseModel):
    """MANE Select transcript structure, as declared in the gene YAML."""

    model_config = ConfigDict(frozen=True)

    id: str
    mane_select: bool = True
    chrom: str
    strand: Literal["+", "-"]
    #: 1-based transcript position of the initiation codon's first base, and of
    #: the terminator's last. Both are ``None`` in a config whose reference has
    #: not been imported yet: the UTR lengths of a RefSeq transcript are a
    #: property of the record, not something to be assumed.
    cds_start_tx: int | None = None
    cds_end_tx: int | None = None
    #: Declared, and cross-checked against the built transcript. Present so a
    #: config error is caught by arithmetic rather than by a reviewer.
    cds_length: int
    protein_length: int
    exon_count: int
    #: Exon labels in transcript order. BRCA1 legitimately skips "4".
    exon_labels: tuple[str, ...] = ()
    #: Inline genomic exon coordinates, written by ``vus-foresight reference import``.
    exons: tuple[Exon, ...] = ()

    @model_validator(mode="after")
    def _consistent(self) -> "TranscriptConfig":
        if self.cds_start_tx is not None and self.cds_end_tx is not None:
            span = self.cds_end_tx - self.cds_start_tx + 1
            if self.cds_length != span:
                raise ValueError(
                    f"{self.id}: declared cds_length {self.cds_length} disagrees with "
                    f"cds_start_tx/cds_end_tx span {span}"
                )
        elif (self.cds_start_tx is None) != (self.cds_end_tx is None):
            raise ValueError(f"{self.id}: cds_start_tx and cds_end_tx must be set together")
        if self.cds_length % 3 != 0:
            raise ValueError(f"{self.id}: cds_length {self.cds_length} is not a multiple of 3")
        if self.cds_length // 3 - 1 != self.protein_length:
            raise ValueError(
                f"{self.id}: cds_length {self.cds_length} implies "
                f"{self.cds_length // 3 - 1} residues, not {self.protein_length}"
            )
        if self.exon_labels and len(self.exon_labels) != self.exon_count:
            raise ValueError(
                f"{self.id}: {len(self.exon_labels)} exon labels for "
                f"exon_count {self.exon_count}"
            )
        if self.exons and len(self.exons) != self.exon_count:
            raise ValueError(
                f"{self.id}: {len(self.exons)} exon records for exon_count {self.exon_count}"
            )
        return self

    @property
    def has_coordinates(self) -> bool:
        return bool(self.exons) and self.cds_start_tx is not None


class FunctionalRegion(BaseModel):
    """A protein region the specification can refer to by name.

    Coordinates are in residues, 1-based inclusive. The engine never looks at
    ``name``; the spec YAML does, via the ``region_in`` predicate.
    """

    model_config = ConfigDict(frozen=True)

    name: str
    start_aa: int
    end_aa: int
    #: Free-form tags the spec can test, e.g. ``critical``, ``hotspot``,
    #: ``repeat`` -- kept open so a new VCEP does not require a code change.
    tags: tuple[str, ...] = ()
    source: str | None = None

    def contains(self, residue: int) -> bool:
        return self.start_aa <= residue <= self.end_aa


class GeneConfig(BaseModel):
    """Everything the engine needs to know about one gene."""

    model_config = ConfigDict(frozen=True)

    gene: str
    assembly: str = "GRCh38"
    transcript: TranscriptConfig
    #: Which VCEP specification applies. Resolved against ``config/specs``.
    spec: str
    #: Whether loss of function is an established disease mechanism for this
    #: gene. This is the gene-level fact PVS1 turns on, and it belongs to the
    #: gene, not to the specification.
    lof_mechanism: Literal["established", "not_established", "unknown"] = "unknown"
    sequence: SequenceResource | None = None
    flanks: SequenceResource | None = None
    functional_regions: tuple[FunctionalRegion, ...] = ()
    #: Named MAVE/SGE datasets and the residue ranges they cover, so the engine
    #: can tell "assayed and benign" from "never assayed".
    mave_datasets: tuple[dict[str, Any], ...] = ()
    #: Width of the intronic window enumerated on each side of every junction.
    intronic_flank_bp: int = 50
    notes: str | None = None

    @property
    def spec_path_hint(self) -> str:
        return f"{self.spec}.yaml"


def read_fasta(path: Path) -> str:
    """Read a single-record FASTA and return the sequence, uppercased.

    Multi-record files are rejected: a transcript resource with two records is
    ambiguous, and guessing which one is meant is how the wrong isoform ends up
    in production.
    """
    records: list[list[str]] = []
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append([])
                continue
            if not records:
                raise ReferenceUnavailable(f"{path}: sequence data before any header")
            records[-1].append(line)
    if len(records) != 1:
        raise ReferenceUnavailable(
            f"{path}: expected exactly one FASTA record, found {len(records)}"
        )
    return "".join(records[0]).upper()


def load_gene_config(path: str | Path) -> GeneConfig:
    """Load and validate a gene YAML."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: gene config must be a mapping")
    return GeneConfig.model_validate(raw)


def build_transcript(config: GeneConfig, *, data_root: str | Path = ".") -> Transcript:
    """Materialise a :class:`Transcript` from a gene config plus its reference.

    Raises :class:`ReferenceUnavailable` when coordinates or sequence are not
    present. That is a deliberate hard failure -- see the module docstring.
    """
    tc = config.transcript
    if not tc.has_coordinates:
        raise ReferenceUnavailable(
            f"{config.gene}: gene config declares no exon coordinates. Run "
            f"'vus-foresight reference import' to populate them from MANE Select; "
            f"see docs/reference-data.md."
        )
    if config.sequence is None:
        raise ReferenceUnavailable(
            f"{config.gene}: gene config declares no transcript sequence resource."
        )
    sequence = config.sequence.read(Path(data_root))
    transcript = Transcript(
        transcript_id=tc.id,
        gene=config.gene,
        chrom=tc.chrom,
        strand=tc.strand,
        exons=tc.exons,
        cds_start_tx=tc.cds_start_tx,
        cds_end_tx=tc.cds_end_tx,
        sequence=sequence,
    )
    if transcript.cds_length != tc.cds_length:
        raise ReferenceUnavailable(
            f"{config.gene}: built CDS length {transcript.cds_length} disagrees with "
            f"the declared {tc.cds_length}"
        )
    return transcript


def load_flanks(
    config: GeneConfig, *, data_root: str | Path = ".", required: bool = False
) -> FlankSequences:
    """Load intronic flanking bases, if the gene config declares them.

    The file is a two-column TSV of ``genomic_position<TAB>plus_strand_base``,
    or a JSON object of the same mapping.

    **Absent flanks are not an error.** The intronic enumeration still yields
    every position -- its cardinality does not depend on the reference base --
    it simply cannot name a reference allele, and marks those rows
    ``reference_base_unknown`` instead of guessing one. A gene config declaring
    a flank resource that has not been materialised is therefore a degraded run,
    not a failed one, and refusing to enumerate at all would contradict the
    design the marker exists to serve.

    Pass ``required=True`` where the caller genuinely cannot proceed without
    them.
    """
    if config.flanks is None:
        return FlankSequences()
    target = config.flanks.resolve(Path(data_root))
    if not target.exists():
        if required:
            raise ReferenceUnavailable(f"declared flank resource {target} is missing")
        return FlankSequences()
    if target.suffix == ".json":
        with target.open("r", encoding="ascii") as handle:
            raw = json.load(handle)
        return FlankSequences(bases={int(k): v.upper() for k, v in raw.items()})
    bases: dict[int, str] = {}
    with target.open("r", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pos, base = line.split("\t")[:2]
            bases[int(pos)] = base.upper()
    return FlankSequences(bases=bases)
