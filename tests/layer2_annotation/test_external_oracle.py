"""Layer 2: the external oracles of spec section 11.

Two fixtures, both optional, both skipped when absent so that CI stays green on
a clone with no curated data:

``tests/fixtures/clinvar_hgvs_oracle.tsv``
    A few thousand real BRCA1/2 variants with curated ``c.`` and ``p.``
    descriptions on the MANE transcript. The system must reproduce them
    character for character. Any divergence is a bug until proven otherwise.

``tests/fixtures/reference_validator_sample.tsv``
    A smaller sample cross-checked against VariantValidator or Mutalyzer.
    *Systematic* divergence indicates a convention error -- parentheses around a
    predicted effect, ``Ter`` versus ``*``, frameshift numbering. *Sparse*
    divergence indicates a coordinate error. The test reports which, because the
    two have completely different fixes.

Format for both: a header row, then ``gene<TAB>hgvs_c<TAB>hgvs_p``.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vus_foresight.annotate.annotator import CodingEdit, annotate_coding_edits

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
CLINVAR_ORACLE = FIXTURES / "clinvar_hgvs_oracle.tsv"
VALIDATOR_SAMPLE = FIXTURES / "reference_validator_sample.tsv"


def _read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")]


def _annotate_substitution(transcript, hgvs_c: str):
    """Re-derive an annotation from a ``c.NNNr>a`` description."""
    body = hgvs_c.removeprefix("c.")
    if ">" not in body:
        return None
    left, alt = body.split(">")
    position_text, ref = left[:-1], left[-1]
    if not position_text.isdigit():
        return None  # intronic or UTR; protein oracle does not apply
    position = int(position_text)
    if not 1 <= position <= transcript.cds_length:
        return None
    if transcript.cds[position - 1] != ref:
        pytest.fail(
            f"{hgvs_c}: oracle says the reference base is {ref}, the transcript has "
            f"{transcript.cds[position - 1]}. This is a coordinate or reference error."
        )
    return annotate_coding_edits(
        transcript, (CodingEdit(cds_position=position, ref=ref, alt=alt),)
    )


@pytest.mark.requires_reference
def test_reproduces_curated_clinvar_hgvs_character_for_character(brca1, brca2):
    if not CLINVAR_ORACLE.exists():
        pytest.skip(f"{CLINVAR_ORACLE.name} not present; see the module docstring")
    transcripts = {"BRCA1": brca1[1], "BRCA2": brca2[1]}

    mismatches: list[tuple[str, str, str]] = []
    compared = 0
    for row in _read(CLINVAR_ORACLE):
        transcript = transcripts.get(row["gene"])
        if transcript is None:
            continue
        variant = _annotate_substitution(transcript, row["hgvs_c"])
        if variant is None:
            continue
        compared += 1
        if variant.hgvs_c != row["hgvs_c"]:
            mismatches.append((row["hgvs_c"], "hgvs_c", variant.hgvs_c))
        expected_p = row.get("hgvs_p", "").strip()
        if expected_p and variant.hgvs_p != expected_p:
            mismatches.append((row["hgvs_c"], expected_p, variant.hgvs_p or ""))

    assert compared > 0, "the oracle file contains no usable coding substitutions"
    assert not mismatches, (
        f"{len(mismatches)} of {compared} curated descriptions were not reproduced. "
        f"First ten: {mismatches[:10]}"
    )


@pytest.mark.requires_reference
def test_divergence_against_a_reference_validator_is_diagnosed(brca1, brca2):
    if not VALIDATOR_SAMPLE.exists():
        pytest.skip(f"{VALIDATOR_SAMPLE.name} not present; see the module docstring")
    transcripts = {"BRCA1": brca1[1], "BRCA2": brca2[1]}

    compared = 0
    protein_mismatches: list[tuple[str, str, str]] = []
    for row in _read(VALIDATOR_SAMPLE):
        transcript = transcripts.get(row["gene"])
        if transcript is None:
            continue
        variant = _annotate_substitution(transcript, row["hgvs_c"])
        if variant is None or not row.get("hgvs_p"):
            continue
        compared += 1
        if variant.hgvs_p != row["hgvs_p"]:
            protein_mismatches.append((row["hgvs_c"], row["hgvs_p"], variant.hgvs_p or ""))

    assert compared > 0
    rate = len(protein_mismatches) / compared
    diagnosis = (
        "SYSTEMATIC -- suspect an HGVS convention (parentheses, Ter vs *, "
        "frameshift numbering)"
        if rate > 0.2
        else "SPARSE -- suspect individual coordinate errors"
    )
    assert not protein_mismatches, (
        f"{len(protein_mismatches)}/{compared} ({rate:.1%}) protein descriptions "
        f"diverge. Diagnosis: {diagnosis}. First ten: {protein_mismatches[:10]}"
    )
