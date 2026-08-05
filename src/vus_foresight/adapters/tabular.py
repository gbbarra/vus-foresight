"""File-backed adapters for the external sources of spec section 8.

All of them share one shape: a local, versioned, tab-separated snapshot keyed by
something the enumerator already knows about a variant. The key differs per
source and that is the only thing the subclasses change.

Column values are parsed once at load time into the nested structure a rule
addresses -- ``frequency.gnomad.popmax_af``, ``functional.classification`` --
so a specification never sees a column name.

Frequency deserves a note. **ABraOM is reported in its own sub-namespace and
never substituted for gnomAD.** A Brazilian allele frequency answers a different
question from a global one; silently merging them would let a variant common in
Brazil and absent from gnomAD be scored as if the world had looked and found
nothing.
"""

from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path
from typing import Any, ClassVar

from ..genome.transcript import Transcript
from ..variant import Variant
from .base import Adapter

__all__ = [
    "FrequencyAdapter",
    "FunctionalAdapter",
    "PredictorAdapter",
    "SpliceAdapter",
    "TableAdapter",
]


def _number(value: str) -> Any:
    """Parse a cell into the narrowest sensible type; empty means absent."""
    text = value.strip()
    if text == "" or text.upper() in {"NA", "NULL", "."}:
        return None
    lowered = text.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


class TableAdapter(Adapter):
    """A TSV snapshot keyed by one variant attribute.

    The file's first column is the key; every other column becomes a field.
    A column named ``a.b`` nests as ``{"a": {"b": ...}}``, which is how a single
    flat file populates ``frequency.gnomad.popmax_af`` and
    ``frequency.abraom.af`` without either shadowing the other.
    """

    #: Attribute of :class:`~vus_foresight.variant.Variant` used as the key.
    key_attribute: ClassVar[str] = "grch38_pos"

    def __init__(
        self,
        path: str | Path | None,
        version: str,
        *,
        key_fn: Callable[[Variant, Transcript], str | None] | None = None,
        default_payload: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(version)
        self._key_fn = key_fn
        self._rows: dict[str, dict[str, Any]] = {}
        #: What this source asserts about a variant it does *not* list.
        #:
        #: The distinction matters more than it looks. gnomAD not listing a
        #: variant is a positive statement -- it was looked for and not seen,
        #: which is exactly the evidence PM2 rests on. A functional assay not
        #: listing a variant is the opposite: nothing was measured. Sources that
        #: mean the first pass a default payload; sources that mean the second
        #: pass nothing and their criteria come out NOT_EVALUABLE.
        self.default_payload = default_payload or {}
        self.path = Path(path) if path is not None else None
        if self.path is not None:
            self._load(self.path)

    def _load(self, path: Path) -> None:
        if not path.exists():
            raise FileNotFoundError(
                f"{type(self).__name__}: snapshot {path} is missing. Adapters never "
                "fetch at runtime; materialise the snapshot first."
            )
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                raise ValueError(f"{path}: snapshot has no header row")
            key_column, *value_columns = header
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                payload: dict[str, Any] = {}
                for column, cell in zip(value_columns, row[1:]):
                    value = _number(cell)
                    if value is None:
                        continue
                    node = payload
                    *parents, leaf = column.split(".")
                    for parent in parents:
                        node = node.setdefault(parent, {})
                    node[leaf] = value
                self._rows[row[0]] = payload

    def key_for(self, variant: Variant, transcript: Transcript) -> str | None:
        if self._key_fn is not None:
            return self._key_fn(variant, transcript)
        return getattr(variant, self.key_attribute, None)

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        key = self.key_for(variant, transcript)
        if key is None:
            return {}
        row = self._rows.get(key)
        if row is None:
            return _deep_copy(self.default_payload)
        merged = _deep_copy(self.default_payload)
        _deep_merge(merged, row)
        return merged

    def __len__(self) -> int:
        return len(self._rows)


def _deep_copy(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        key: _deep_copy(value) if isinstance(value, dict) else value
        for key, value in payload.items()
    }


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> None:
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = value


class FrequencyAdapter(TableAdapter):
    """Population frequency. Columns like ``gnomad.popmax_af``, ``abraom.af``."""

    namespace = "frequency"
    name = "gnomad"
    key_attribute = "grch38_pos"


class PredictorAdapter(TableAdapter):
    """Computational predictors, keyed genomically.

    The specification names which predictor to use and at which threshold; this
    adapter has no opinion. Loading REVEL and then scoring with a BayesDel
    threshold is a specification error, and it is caught by the spec declaring
    the field it requires.
    """

    namespace = "predictor"
    name = "dbnsfp"
    key_attribute = "grch38_pos"


class SpliceAdapter(TableAdapter):
    """Pre-computed splicing predictions, keyed genomically.

    Kept separate from :class:`PredictorAdapter` because it covers every
    possible SNV rather than only missense, and because switching it off is the
    single most informative counterfactual the map supports.
    """

    namespace = "splice"
    name = "spliceai"
    key_attribute = "grch38_pos"


class FunctionalAdapter(TableAdapter):
    """MAVE/SGE results, keyed by protein change.

    Provenance is per variant: BRCA1 and BRCA2 have separate, independent
    datasets, and a row must record which one produced its score. The adapter
    also publishes ``region_assayed`` so the gap engine can tell "assayed, no
    effect" from "nobody has ever assayed this residue" -- the difference
    between ``available_uningested`` and ``assay_feasible``.
    """

    namespace = "functional"
    name = "mavedb"
    key_attribute = "hgvs_p"

    def __init__(
        self,
        path: str | Path | None,
        version: str,
        *,
        assayed_regions: tuple[tuple[int, int, str], ...] = (),
        key_fn: Callable[[Variant, Transcript], str | None] | None = None,
        default_payload: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(path, version, key_fn=key_fn, default_payload=default_payload)
        #: ``(first_residue, last_residue, dataset)`` triples covered by an assay.
        self.assayed_regions = assayed_regions

    def lookup(self, variant: Variant, transcript: Transcript) -> dict[str, Any]:
        payload = super().lookup(variant, transcript)
        residue = variant.codon_index
        if residue is None:
            payload.setdefault("region_assayed", False)
            return payload
        covering = [
            dataset for first, last, dataset in self.assayed_regions if first <= residue <= last
        ]
        payload["region_assayed"] = bool(covering)
        if covering:
            payload.setdefault("assay_datasets", sorted(covering))
        return payload


# ClinVar deliberately has no TableAdapter subclass. PS1 and PM5 cannot be
# answered by a single-key lookup: both require excluding the variant's own
# record, and PM5 additionally requires the neighbouring record to be a
# missense. A flat table keyed by protein change would let a variant satisfy
# PS1 from itself. See :mod:`vus_foresight.adapters.clinvar`.
