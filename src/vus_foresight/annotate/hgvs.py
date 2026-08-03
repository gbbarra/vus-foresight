"""HGVS string construction.

Conventions fixed here once, because a systematic convention error (``Ter`` vs
``*``, parentheses around a predicted effect, frameshift numbering) shows up as
a *systematic* divergence against an external reference and is therefore
diagnosable -- see the secondary-oracle test in spec section 11.

Choices made, and why:

* Three-letter amino acids (``p.Leu100Pro``), matching ClinVar's display form
  for MANE transcripts.
* No parentheses around predicted protein effects. Every protein consequence
  this system emits is predicted, so parenthesising all of them would carry no
  information while guaranteeing a character-level mismatch against ClinVar.
* Synonymous rendered ``p.Leu100=``.
* Frameshift rendered ``p.Arg100SerfsTer12``, where 12 counts the new codon
  through to and including the terminator.
"""

from __future__ import annotations

from ..genome.sequence import three_letter

__all__ = [
    "format_substitution",
    "format_delins",
    "format_deletion",
    "format_protein_substitution",
    "format_protein_synonymous",
    "format_protein_start_lost",
    "format_protein_stop_lost",
    "format_protein_deletion",
    "format_protein_delins",
    "format_protein_frameshift",
    "format_exon_cnv",
]


def format_substitution(c_pos: str, ref: str, alt: str) -> str:
    """``c.100A>G``"""
    return f"c.{c_pos}{ref}>{alt}"


def format_delins(c_start: str, c_end: str, alt: str) -> str:
    """``c.100_102delinsGCT`` (or ``c.100delinsGC`` for a single-base span)."""
    span = c_start if c_start == c_end else f"{c_start}_{c_end}"
    return f"c.{span}delins{alt}"


def format_deletion(c_start: str, c_end: str) -> str:
    """``c.100_102del``"""
    span = c_start if c_start == c_end else f"{c_start}_{c_end}"
    return f"c.{span}del"


def format_protein_substitution(ref_aa: str, residue: int, alt_aa: str) -> str:
    """``p.Leu100Pro`` / ``p.Arg1699Ter``"""
    return f"p.{three_letter(ref_aa)}{residue}{three_letter(alt_aa)}"


def format_protein_synonymous(ref_aa: str, residue: int) -> str:
    """``p.Leu100=``"""
    return f"p.{three_letter(ref_aa)}{residue}="


def format_protein_start_lost() -> str:
    """``p.Met1?`` -- the effect of losing the initiation codon is not predictable."""
    return "p.Met1?"


def format_protein_stop_lost(residue: int, alt_aa: str) -> str:
    """``p.Ter1864LeuextTer?``

    The extension length is unknown without the 3' UTR reading frame, so it is
    reported as ``?`` rather than guessed.
    """
    return f"p.Ter{residue}{three_letter(alt_aa)}extTer?"


def format_protein_deletion(
    ref_aa_first: str, residue_first: int, ref_aa_last: str, residue_last: int
) -> str:
    """``p.Leu100del`` / ``p.Leu100_Pro101del``"""
    if residue_first == residue_last:
        return f"p.{three_letter(ref_aa_first)}{residue_first}del"
    return (
        f"p.{three_letter(ref_aa_first)}{residue_first}_"
        f"{three_letter(ref_aa_last)}{residue_last}del"
    )


def format_protein_delins(
    ref_aa_first: str,
    residue_first: int,
    ref_aa_last: str,
    residue_last: int,
    inserted: str,
) -> str:
    """``p.Leu100_Pro101delinsGln``"""
    ins = "".join(three_letter(aa) for aa in inserted)
    if residue_first == residue_last:
        return f"p.{three_letter(ref_aa_first)}{residue_first}delins{ins}"
    return (
        f"p.{three_letter(ref_aa_first)}{residue_first}_"
        f"{three_letter(ref_aa_last)}{residue_last}delins{ins}"
    )


def format_protein_frameshift(
    ref_aa: str, residue: int, alt_aa: str, new_codons_to_stop: int | None
) -> str:
    """``p.Arg100SerfsTer12``, or ``p.Arg100fs`` when the new residue is unknown.

    ``new_codons_to_stop`` counts the first changed codon through the
    terminator inclusive, per HGVS. ``None`` renders ``fsTer?``.
    """
    if alt_aa is None:
        return f"p.{three_letter(ref_aa)}{residue}fs"
    tail = "?" if new_codons_to_stop is None else str(new_codons_to_stop)
    return f"p.{three_letter(ref_aa)}{residue}{three_letter(alt_aa)}fsTer{tail}"


def format_exon_cnv(gene: str, first_label: str, last_label: str, kind: str) -> str:
    """Exon-level CNV pseudo-HGVS, e.g. ``BRCA1 exon 3-5 deletion``.

    Exon-level copy number changes are scored by a different framework
    (ClinGen CNV, Riggs et al.) and have no single ``c.`` description without
    breakpoint resolution, which is precisely the information a curated exon
    call does not carry. Naming them explicitly is more honest than inventing
    coordinates.
    """
    span = first_label if first_label == last_label else f"{first_label}-{last_label}"
    return f"{gene} exon {span} {kind}"
