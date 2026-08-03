"""Nucleotide and amino-acid primitives.

Everything here is pure sequence arithmetic with no notion of a gene. The
standard genetic code is hard-coded because it is not domain configuration --
but the *selenocysteine* and alternative-code cases are deliberately absent: no
BRCA-family gene uses one, and silently supporting them would hide the fact that
a future gene needs an explicit decision.
"""

from __future__ import annotations

__all__ = [
    "BASES",
    "CODON_TABLE",
    "AA_THREE_LETTER",
    "complement",
    "reverse_complement",
    "translate_codon",
    "translate",
    "alternatives",
    "three_letter",
]

#: Canonical base ordering. Every enumeration that iterates alternative alleles
#: uses this order, which is one of the reasons the output is byte-reproducible.
BASES: tuple[str, ...] = ("A", "C", "G", "T")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_THREE_LETTER: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}


def complement(seq: str) -> str:
    """Base-wise complement, preserving case and passing ``N`` through."""
    return seq.translate(_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    """Reverse complement. The only place strand is ever flipped."""
    return seq.translate(_COMPLEMENT)[::-1]


def translate_codon(codon: str) -> str:
    """Translate exactly three bases.

    Returns ``"X"`` for any codon containing a character outside ``ACGT``.
    Raises on a wrong-length input, because a two-base "codon" is always a
    coordinate bug upstream and must not be silently tolerated.
    """
    if len(codon) != 3:
        raise ValueError(f"codon must be 3 bases, got {len(codon)}: {codon!r}")
    return CODON_TABLE.get(codon.upper(), "X")


def translate(seq: str, *, stop_at_terminator: bool = False) -> str:
    """Translate a coding sequence in frame 0.

    A trailing partial codon is ignored. This is the *independent* translation
    path used by the layer-2 tests: it never consults the annotator.
    """
    peptide: list[str] = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = translate_codon(seq[i : i + 3])
        peptide.append(aa)
        if stop_at_terminator and aa == "*":
            break
    return "".join(peptide)


def alternatives(ref: str) -> tuple[str, ...]:
    """The three alternative bases for a reference base, in canonical order."""
    ref = ref.upper()
    if ref not in BASES:
        raise ValueError(f"reference base must be one of {BASES}, got {ref!r}")
    return tuple(b for b in BASES if b != ref)


def three_letter(aa: str) -> str:
    """One-letter amino acid to HGVS three-letter form."""
    try:
        return AA_THREE_LETTER[aa]
    except KeyError as exc:  # pragma: no cover - defensive
        raise ValueError(f"unknown amino acid code {aa!r}") from exc
