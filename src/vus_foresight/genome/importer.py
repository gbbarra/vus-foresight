"""Populate a gene config's coordinates from MANE Select resources.

Run once per gene, by a curator, with the reference files in hand. The result is
written back into the gene YAML as an explicit, reviewable block -- coordinates
belong in version control where a change to them shows up in a diff, not in a
cache that silently refreshes.

The CDS boundaries are *derived and checked*, not asked for: given the declared
CDS length, there is normally exactly one open reading frame in the transcript
that starts with ATG, ends with a terminator, and contains no internal stop. If
there is not exactly one, the import fails rather than picking.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import yaml

from .reference import GeneConfig, ReferenceUnavailable, read_fasta
from .sequence import translate_codon
from .transcript import Exon

__all__ = ["ImportedReference", "read_exon_table", "locate_cds", "import_reference"]


class ImportedReference(dict):
    """The block written back into the gene YAML."""


def read_exon_table(path: str | Path) -> tuple[Exon, ...]:
    """Read a TSV of ``label<TAB>start<TAB>end`` in **transcript** order.

    Transcript order, not genomic order: for a minus-strand gene the coordinates
    descend, and silently sorting them would produce a transcript that is
    internally consistent and biologically wrong.
    """
    path = Path(path)
    if not path.exists():
        raise ReferenceUnavailable(f"exon table {path} is missing")
    exons: list[Exon] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0].lower() in {"label", "exon"}:
                continue
            label, start, end = row[0], int(row[1]), int(row[2])
            exons.append(Exon(label=label, start=start, end=end))
    if not exons:
        raise ReferenceUnavailable(f"exon table {path} contains no exons")
    return tuple(exons)


def locate_cds(sequence: str, cds_length: int) -> tuple[int, int]:
    """Find the unique ORF of exactly ``cds_length`` bases. 1-based, inclusive.

    Raises when zero or more than one candidate exists -- both mean the declared
    length disagrees with the sequence, and guessing would be exactly the kind of
    plausible-looking error this project is built to avoid.
    """
    candidates: list[int] = []
    for start in range(1, len(sequence) - cds_length + 2):
        if sequence[start - 1 : start + 2] != "ATG":
            continue
        window = sequence[start - 1 : start - 1 + cds_length]
        if translate_codon(window[-3:]) != "*":
            continue
        internal = any(
            translate_codon(window[i : i + 3]) == "*" for i in range(0, cds_length - 3, 3)
        )
        if internal:
            continue
        candidates.append(start)
    if len(candidates) != 1:
        raise ReferenceUnavailable(
            f"expected exactly one open reading frame of {cds_length} nt, found "
            f"{len(candidates)} at {candidates[:5]}. Check the declared cds_length "
            "and that the FASTA is the spliced transcript, not the genomic locus."
        )
    start = candidates[0]
    return start, start + cds_length - 1


def import_reference(
    config: GeneConfig,
    *,
    sequence_path: str | Path,
    exon_table_path: str | Path,
    config_path: str | Path | None = None,
) -> ImportedReference:
    """Derive the transcript block and, optionally, write it back to the YAML."""
    sequence = read_fasta(Path(sequence_path))
    exons = read_exon_table(exon_table_path)

    tc = config.transcript
    if len(exons) != tc.exon_count:
        raise ReferenceUnavailable(
            f"{config.gene}: exon table has {len(exons)} exons, config declares "
            f"{tc.exon_count}"
        )
    if tc.exon_labels and tuple(e.label for e in exons) != tc.exon_labels:
        raise ReferenceUnavailable(
            f"{config.gene}: exon labels in the table do not match the config. "
            f"Table: {[e.label for e in exons]}"
        )
    spliced_length = sum(e.length for e in exons)
    if spliced_length != len(sequence):
        raise ReferenceUnavailable(
            f"{config.gene}: exons sum to {spliced_length} nt but the FASTA holds "
            f"{len(sequence)} nt"
        )

    cds_start, cds_end = locate_cds(sequence, tc.cds_length)

    block = ImportedReference(
        {
            "cds_start_tx": cds_start,
            "cds_end_tx": cds_end,
            "exons": [
                {"label": e.label, "start": e.start, "end": e.end} for e in exons
            ],
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "sequence_length": len(sequence),
        }
    )

    if config_path is not None:
        path = Path(config_path)
        with path.open("r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle)
        raw["transcript"]["cds_start_tx"] = block["cds_start_tx"]
        raw["transcript"]["cds_end_tx"] = block["cds_end_tx"]
        raw["transcript"]["exons"] = block["exons"]
        raw.setdefault("sequence", {})["sha256"] = block["sequence_sha256"]
        raw["sequence"]["expected_length"] = block["sequence_length"]
        with path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(raw, handle, sort_keys=False, allow_unicode=True)
    return block
