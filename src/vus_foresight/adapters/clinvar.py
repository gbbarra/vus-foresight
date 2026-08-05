"""Dated ClinVar snapshots for the semi-intrinsic criteria (spec section 4).

PS1 and PM5 are the only criteria whose answer depends on the state of a public
database on a date rather than on the variant itself. They are the direct hook
into ``vus-hindsight``: a variant can leave VUS without anything new being
learned *about it*, purely because a different variant in the same codon got
classified.

Two exclusions make the difference between a criterion and a tautology, and
neither is expressible in a flat keyed table:

* **PS1 must exclude the variant's own record.** The criterion is "the same
  amino acid change, reached by a *different* nucleotide change, is already
  established pathogenic". A variant that is itself pathogenic in ClinVar
  satisfying PS1 from its own record is circular, and it would silently
  manufacture four points for every already-classified variant in the gene.
* **PM5 must exclude the same protein change**, which is PS1's business, and
  the neighbouring record has to be a missense -- "a different *missense*
  change at this residue". A nonsense variant in the same codon is a PVS1
  observation, not a PM5 one.

So the snapshot is indexed three ways and resolved per variant, rather than
looked up by a single key.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

from ..genome.transcript import Transcript
from ..variant import Variant
from .base import Adapter

__all__ = [
    "CLINVAR_REVIEW_STARS",
    "ClinVarRecord",
    "ClinVarSnapshot",
    "ClinVarSnapshotAdapter",
    "normalise_classification",
    "protein_change_kind",
    "review_stars",
]

#: ClinVar's documented review-status strings, mapped to the star rating the
#: archive displays. Anything unrecognised scores zero rather than raising: a
#: new status string must not be able to promote a record.
CLINVAR_REVIEW_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "no interpretation for the single variant": 0,
}

#: Normalised significance values, strongest assertion first.
_CLASSIFICATION_ALIASES: dict[str, str] = {
    "pathogenic": "pathogenic",
    "pathogenic/likely pathogenic": "pathogenic",
    "likely pathogenic": "likely_pathogenic",
    "likely pathogenic, low penetrance": "likely_pathogenic",
    "uncertain significance": "uncertain",
    "uncertain risk allele": "uncertain",
    "conflicting classifications of pathogenicity": "conflicting",
    "conflicting interpretations of pathogenicity": "conflicting",
    "likely benign": "likely_benign",
    "benign/likely benign": "benign",
    "benign": "benign",
}

#: Precedence when several records back the same protein change or codon. The
#: strongest *pathogenic* assertion wins, because both PS1 and PM5 ask whether a
#: pathogenic precedent exists, not what the consensus is.
_CLASSIFICATION_RANK: dict[str, int] = {
    "pathogenic": 0,
    "likely_pathogenic": 1,
    "benign": 2,
    "likely_benign": 3,
    "conflicting": 4,
    "uncertain": 5,
}


def normalise_classification(value: str) -> str | None:
    """Map a ClinVar significance string onto the vocabulary a spec can test.

    Unrecognised values return ``None`` and the record is dropped: a rule must
    never fire on a string nobody has looked at.
    """
    return _CLASSIFICATION_ALIASES.get(value.strip().lower())


def review_stars(value: str) -> int:
    return CLINVAR_REVIEW_STARS.get(value.strip().lower(), 0)


def protein_change_kind(hgvs_p: str | None) -> str:
    """Classify a ``p.`` description without re-deriving it from the sequence.

    Only the distinctions PS1 and PM5 need: a terminator, a synonymous change,
    an unpredictable start loss, or a missense.
    """
    if not hgvs_p:
        return "unknown"
    text = hgvs_p.removeprefix("p.").removeprefix("(").removesuffix(")")
    if text.endswith("="):
        return "synonymous"
    if text.endswith("?"):
        return "unpredictable"
    if "fs" in text or "ext" in text:
        return "truncating"
    if text.endswith("Ter") or text.endswith("*"):
        return "nonsense"
    if "del" in text or "ins" in text or "dup" in text:
        return "indel"
    return "missense"


@dataclass(frozen=True, slots=True)
class ClinVarRecord:
    """One classified variant in a snapshot, on the transcript of interest."""

    hgvs_c: str
    hgvs_p: str | None
    classification: str
    stars: int
    codon: int | None = None
    last_evaluated: date | None = None

    @property
    def kind(self) -> str:
        return protein_change_kind(self.hgvs_p)


def _best(records: Iterable[ClinVarRecord]) -> ClinVarRecord | None:
    ordered = sorted(
        records, key=lambda r: (_CLASSIFICATION_RANK.get(r.classification, 9), -r.stars, r.hgvs_c)
    )
    return ordered[0] if ordered else None


@dataclass(slots=True)
class ClinVarSnapshot:
    """A dated set of records, indexed for the two questions PS1 and PM5 ask."""

    snapshot_date: date | None
    by_nucleotide: dict[str, ClinVarRecord]
    by_protein: dict[str, list[ClinVarRecord]]
    by_codon: dict[int, list[ClinVarRecord]]

    @classmethod
    def from_records(
        cls, records: Iterable[ClinVarRecord], snapshot_date: date | None = None
    ) -> ClinVarSnapshot:
        by_nucleotide: dict[str, ClinVarRecord] = {}
        by_protein: dict[str, list[ClinVarRecord]] = {}
        by_codon: dict[int, list[ClinVarRecord]] = {}
        for record in records:
            by_nucleotide[record.hgvs_c] = record
            if record.hgvs_p:
                by_protein.setdefault(record.hgvs_p, []).append(record)
            if record.codon is not None:
                by_codon.setdefault(record.codon, []).append(record)
        return cls(
            snapshot_date=snapshot_date,
            by_nucleotide=by_nucleotide,
            by_protein=by_protein,
            by_codon=by_codon,
        )

    @classmethod
    def read(cls, path: str | Path, snapshot_date: date | None = None) -> ClinVarSnapshot:
        """Read the normalised snapshot TSV written by ``clinvar build``.

        Columns: ``hgvs_c``, ``hgvs_p``, ``codon``, ``classification``,
        ``stars``, ``last_evaluated``.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(
                f"ClinVar snapshot {path} is missing. Adapters never fetch at "
                "runtime; build the snapshot first with 'vus-foresight clinvar build'."
            )
        records: list[ClinVarRecord] = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                codon = row.get("codon") or ""
                evaluated = row.get("last_evaluated") or ""
                records.append(
                    ClinVarRecord(
                        hgvs_c=row["hgvs_c"],
                        hgvs_p=(row.get("hgvs_p") or "") or None,
                        classification=row["classification"],
                        stars=int(row.get("stars") or 0),
                        codon=int(codon) if codon else None,
                        last_evaluated=date.fromisoformat(evaluated) if evaluated else None,
                    )
                )
        return cls.from_records(records, snapshot_date)

    def __len__(self) -> int:
        return len(self.by_nucleotide)

    def resolve(self, variant: Variant, *, min_stars: int) -> dict[str, Any]:
        """Everything the semi-intrinsic criteria may read about one variant."""
        payload: dict[str, Any] = {}
        if self.snapshot_date is not None:
            payload["snapshot_date"] = self.snapshot_date.isoformat()

        own = self.by_nucleotide.get(variant.hgvs_c)
        if own is not None:
            # Reported for provenance and for the hindsight join. No criterion
            # reads it: classifying a variant because ClinVar already did would
            # make the map a mirror rather than a measurement.
            payload["self"] = {
                "classification": own.classification,
                "stars": own.stars,
            }

        if variant.hgvs_p:
            peers = [
                record
                for record in self.by_protein.get(variant.hgvs_p, ())
                if record.hgvs_c != variant.hgvs_c and record.stars >= min_stars
            ]
            best = _best(peers)
            if best is not None:
                payload["same_protein_change"] = {
                    "classification": best.classification,
                    "stars": best.stars,
                    "count": len(peers),
                    "example": best.hgvs_c,
                }

        if variant.codon_index is not None:
            neighbours = [
                record
                for record in self.by_codon.get(variant.codon_index, ())
                if record.hgvs_p != variant.hgvs_p
                and record.stars >= min_stars
                and record.kind == "missense"
            ]
            best = _best(neighbours)
            if best is not None:
                payload["codon"] = {
                    "other_change_classification": best.classification,
                    "stars": best.stars,
                    "count": len(neighbours),
                    "example": best.hgvs_p,
                }
        return payload


class ClinVarSnapshotAdapter(Adapter):
    """Serves the ``clinvar`` namespace from a dated snapshot.

    ``min_stars`` is the review-status floor a neighbouring record must clear
    before it can support PS1 or PM5. It defaults to 1 -- assertion criteria
    provided -- because a submission with no stated criteria is not a precedent
    anyone would cite.
    """

    namespace = "clinvar"
    name = "clinvar"

    def __init__(
        self,
        snapshot: ClinVarSnapshot,
        version: str,
        *,
        min_stars: int = 1,
    ) -> None:
        super().__init__(version)
        self.snapshot = snapshot
        self.min_stars = min_stars

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        snapshot_date: date | None = None,
        min_stars: int = 1,
    ) -> ClinVarSnapshotAdapter:
        version = snapshot_date.isoformat() if snapshot_date else Path(path).stem
        return cls(ClinVarSnapshot.read(path, snapshot_date), version, min_stars=min_stars)

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        return self.snapshot.resolve(variant, min_stars=self.min_stars)

    def describe(self) -> dict[str, str]:
        return {
            **super().describe(),
            "records": str(len(self.snapshot)),
            "min_stars": str(self.min_stars),
        }
