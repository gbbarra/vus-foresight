"""Extract one transcript's reference from a MANE Select release.

The MANE release ships a genome-wide RNA FASTA and a genome-wide GTF. This turns
the two into the pair ``reference import`` wants: a single-record FASTA and an
exon table in **transcript order**.

The transcript-order requirement is the whole reason this is code rather than an
``awk`` one-liner in a README. For a minus-strand gene the genomic coordinates
descend, and the obvious ``sort -k2,2n`` produces a transcript that is internally
consistent and biologically wrong -- the failure mode spec section 11 names as
the riskiest in the system. GTF ``exon_number`` is already in transcript order,
so the ordering is taken from there and then *checked* against the strand.

Exon labels come from the gene config, not from the GTF. BRCA1's GTF numbers its
exons 1..23 while the clinical numbering runs 1,2,3,5,...,24, and the difference
is exactly the kind of thing that silently mislabels every downstream report.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import IO

from .reference import GeneConfig, ReferenceUnavailable
from .transcript import Exon

__all__ = ["ManeExtract", "extract_transcript", "normalise_chrom", "open_maybe_gzip"]

_ATTRIBUTE = re.compile(r'(\w+)\s+"([^"]*)"')

#: RefSeq chromosome accessions for the two assemblies we care about, so a GTF
#: keyed by ``NC_000017.11`` still matches a config that says ``chr17``.
_ACCESSION_TO_CHROM = {f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)}
_ACCESSION_TO_CHROM.update({"NC_000023": "chrX", "NC_000024": "chrY", "NC_012920": "chrM"})


def open_maybe_gzip(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def normalise_chrom(seqname: str) -> str:
    """Reduce the several spellings of a chromosome to ``chrN``.

    MANE ships GTFs keyed by RefSeq accession, by bare number and by ``chrN``
    depending on which file you took. Comparing them literally would reject a
    perfectly good reference.
    """
    name = seqname.strip()
    accession = name.split(".")[0]
    if accession in _ACCESSION_TO_CHROM:
        return _ACCESSION_TO_CHROM[accession]
    if name.startswith("chr"):
        return name
    return f"chr{name}"


def _attributes(field: str) -> dict[str, str]:
    return dict(_ATTRIBUTE.findall(field))


@dataclass(frozen=True, slots=True)
class ManeExtract:
    """What a curator needs on disk, plus what was checked to get it."""

    transcript_id: str
    fasta: str
    exons: tuple[Exon, ...]
    chrom: str
    strand: str
    #: Notes worth printing: things that were tolerated rather than rejected.
    warnings: tuple[str, ...] = ()

    def fasta_text(self, *, line_width: int = 60) -> str:
        body = "\n".join(
            self.fasta[i : i + line_width] for i in range(0, len(self.fasta), line_width)
        )
        return f">{self.transcript_id}\n{body}\n"

    def exon_table(self) -> str:
        rows = "".join(f"{e.label}\t{e.start}\t{e.end}\n" for e in self.exons)
        return f"label\tstart\tend\n{rows}"


def _read_fasta_record(path: Path, transcript_id: str) -> str:
    """Pull one record out of a multi-record FASTA by accession."""
    lines: list[str] = []
    keeping = False
    found = False
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if keeping:
                    break
                keeping = line[1:].split()[0] == transcript_id
                found = found or keeping
                continue
            if keeping:
                lines.append(line.strip())
    if not found:
        raise ReferenceUnavailable(
            f"{path}: no FASTA record for {transcript_id!r}. Check that the release "
            "carries this transcript version -- the version is part of the identity."
        )
    return "".join(lines).upper()


def _iter_exon_features(path: Path, transcript_id: str) -> Iterator[tuple[str, str, int, int, int]]:
    """``(seqname, strand, start, end, exon_number)`` for one transcript."""
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attributes = _attributes(fields[8])
            if transcript_id not in attributes.values():
                continue
            number = attributes.get("exon_number")
            if number is None:
                raise ReferenceUnavailable(
                    f"{path}: an exon of {transcript_id} carries no exon_number "
                    "attribute, so transcript order cannot be recovered"
                )
            yield fields[0], fields[6], int(fields[3]), int(fields[4]), int(number)


def extract_transcript(
    config: GeneConfig,
    *,
    gtf_path: str | Path,
    fasta_path: str | Path,
) -> ManeExtract:
    """Extract the transcript named by ``config`` from a MANE release.

    Everything the gene config already declares -- strand, chromosome, exon
    count, exon labels, CDS length -- is used as a cross-check rather than
    trusted from the release. A disagreement is an error, because the two
    disagreeing is precisely the situation where guessing produces coordinates
    that look right.
    """
    transcript_id = config.transcript.id
    gtf_path, fasta_path = Path(gtf_path), Path(fasta_path)

    features = list(_iter_exon_features(gtf_path, transcript_id))
    if not features:
        raise ReferenceUnavailable(
            f"{gtf_path}: no exon features for {transcript_id!r}. If the release is "
            "keyed by Ensembl identifiers, use the RefSeq genomic GTF instead."
        )

    seqnames = {normalise_chrom(f[0]) for f in features}
    strands = {f[1] for f in features}
    if len(seqnames) != 1 or len(strands) != 1:
        raise ReferenceUnavailable(
            f"{transcript_id}: exons span {sorted(seqnames)} on strands "
            f"{sorted(strands)}; a transcript must sit on one strand of one sequence"
        )
    chrom, strand = seqnames.pop(), strands.pop()

    warnings: list[str] = []
    if strand != config.transcript.strand:
        raise ReferenceUnavailable(
            f"{transcript_id}: release says strand {strand!r}, config says "
            f"{config.transcript.strand!r}"
        )
    if chrom != normalise_chrom(config.transcript.chrom):
        raise ReferenceUnavailable(
            f"{transcript_id}: release says {chrom}, config says {config.transcript.chrom}"
        )

    features.sort(key=lambda f: f[4])
    numbers = [f[4] for f in features]
    if numbers != list(range(1, len(numbers) + 1)):
        raise ReferenceUnavailable(
            f"{transcript_id}: exon_number values are {numbers}, not a contiguous "
            "1..n run; transcript order cannot be trusted"
        )

    if len(features) != config.transcript.exon_count:
        raise ReferenceUnavailable(
            f"{transcript_id}: release has {len(features)} exons, config declares "
            f"{config.transcript.exon_count}"
        )

    # exon_number is defined in transcript order, so on the minus strand the
    # genomic coordinates must descend. Checking it here is what stops a
    # differently-conventioned release from silently reversing the gene.
    genomic_starts = [f[2] for f in features]
    descending = all(a > b for a, b in pairwise(genomic_starts))
    ascending = all(a < b for a, b in pairwise(genomic_starts))
    if strand == "+" and not ascending:
        raise ReferenceUnavailable(
            f"{transcript_id}: plus strand but exon_number order is not ascending "
            "in genomic coordinates"
        )
    if strand == "-" and not descending:
        raise ReferenceUnavailable(
            f"{transcript_id}: minus strand but exon_number order is not descending "
            "in genomic coordinates"
        )

    labels = config.transcript.exon_labels or tuple(str(n) for n in numbers)
    if len(labels) != len(features):
        raise ReferenceUnavailable(
            f"{transcript_id}: {len(labels)} exon labels for {len(features)} exons"
        )
    if config.transcript.exon_labels and list(labels) != [str(n) for n in numbers]:
        warnings.append(
            f"clinical exon labels differ from GTF exon_number "
            f"({labels[0]}..{labels[-1]} vs 1..{numbers[-1]}); using the clinical "
            "labels from the gene config, which is the intended behaviour"
        )

    exons = tuple(
        Exon(label=label, start=start, end=end)
        # strict: the count is validated a few lines above, so this is a
        # redundant assertion -- and a free one against a future edit that
        # drops that check. Silently truncating would build a transcript
        # missing exons, which is the failure this module exists to prevent.
        for label, (_seq, _strand, start, end, _n) in zip(labels, features, strict=True)
    )

    sequence = _read_fasta_record(fasta_path, transcript_id)
    spliced = sum(e.length for e in exons)
    if spliced != len(sequence):
        raise ReferenceUnavailable(
            f"{transcript_id}: exons sum to {spliced} nt but the FASTA record holds "
            f"{len(sequence)} nt. The GTF and the FASTA are from different releases."
        )
    if len(sequence) < config.transcript.cds_length:
        raise ReferenceUnavailable(
            f"{transcript_id}: transcript is {len(sequence)} nt, shorter than the "
            f"declared CDS length {config.transcript.cds_length}"
        )

    return ManeExtract(
        transcript_id=transcript_id,
        fasta=sequence,
        exons=exons,
        chrom=chrom,
        strand=strand,
        warnings=tuple(warnings),
    )
