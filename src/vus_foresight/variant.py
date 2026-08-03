"""Variant identity and molecular consequence.

A :class:`Variant` here is a *possible* variant, not an observed one. Nothing in
this package ever sees patient data (spec section 13).
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, field_validator

__all__ = [
    "Consequence",
    "VariantKind",
    "Variant",
    "PROTEIN_TRUNCATING",
    "SPLICE_CONSEQUENCES",
]


class Consequence(str, Enum):
    """Molecular consequence, one primary term per variant.

    Tier 1 of spec section 3 needs the first seven. The remainder are required
    by Tier 2 (MNV, frameshift equivalence classes, in-frame deletions, CNV).
    """

    MISSENSE = "missense"
    NONSENSE = "nonsense"
    SYNONYMOUS = "synonymous"
    SPLICE_DONOR = "splice_donor"
    SPLICE_ACCEPTOR = "splice_acceptor"
    SPLICE_REGION = "splice_region"
    INTRONIC = "intronic"
    START_LOST = "start_lost"
    STOP_LOST = "stop_lost"
    STOP_RETAINED = "stop_retained"
    FRAMESHIFT = "frameshift"
    INFRAME_DELETION = "inframe_deletion"
    INFRAME_INSERTION = "inframe_insertion"
    EXON_DELETION = "exon_deletion"
    EXON_DUPLICATION = "exon_duplication"
    UTR5 = "5_prime_utr"
    UTR3 = "3_prime_utr"


#: Consequences that produce a premature termination or loss of the start codon,
#: i.e. the entry gate for PVS1. Membership is a property of the vocabulary, not
#: of any gene.
PROTEIN_TRUNCATING: frozenset[Consequence] = frozenset(
    {
        Consequence.NONSENSE,
        Consequence.FRAMESHIFT,
        Consequence.SPLICE_DONOR,
        Consequence.SPLICE_ACCEPTOR,
        Consequence.START_LOST,
        Consequence.EXON_DELETION,
    }
)

SPLICE_CONSEQUENCES: frozenset[Consequence] = frozenset(
    {
        Consequence.SPLICE_DONOR,
        Consequence.SPLICE_ACCEPTOR,
        Consequence.SPLICE_REGION,
    }
)


class VariantKind(str, Enum):
    """Which enumeration tier produced this row.

    ``FRAMESHIFT_CLASS`` and the CNV kinds are *classes*, not single variants:
    they stand for a set of concrete alleles that share an evidence profile.
    """

    SNV = "snv"
    MNV = "mnv"
    FRAMESHIFT_CLASS = "frameshift_class"
    INFRAME_DELETION = "inframe_deletion"
    CNV = "cnv"


class Variant(BaseModel):
    """One enumerated variant (or variant class) with its annotation.

    ``grch38_pos`` is a VCF-style ``chrom-pos-ref-alt`` string on the *genomic
    plus strand*, regardless of the transcript's strand. For a minus-strand gene
    such as BRCA1, ``ref``/``alt`` are therefore the reverse complement of the
    bases named in the ``c.`` description -- which is exactly the asymmetry the
    strand test of spec section 11 exploits.

    Two documented departures from strict VCF:

    * deletions use a literal ``del`` as the alt field
      (``chr17-43094500-CTG-del``) because the anchor base a VCF deletion needs
      may be intronic and therefore absent from a transcript sequence;
    * the field is ``None`` when the change spans an exon junction, since such a
      change is one variant in ``c.`` space and two records in genomic space.
      ``attributes["genomic_discontiguous"]`` marks those rows.
    """

    model_config = ConfigDict(frozen=True)

    gene: str
    transcript: str
    kind: VariantKind
    hgvs_c: str
    hgvs_p: str | None = None
    grch38_pos: str | None = None
    consequence: Consequence
    #: Every applicable SO-style term, primary term first. ``consequence`` is
    #: ``consequence_terms[0]``; the rest carry secondary facts such as a
    #: missense that also sits in the splice region.
    consequence_terms: tuple[Consequence, ...] = ()
    #: 1 for an SNV, 2-3 for an intra-codon MNV, ``0`` for class-level rows.
    mutational_distance: int = 1
    #: 1-based codon index within the CDS, when the variant touches coding
    #: sequence. ``None`` for purely intronic variants.
    codon_index: int | None = None
    #: 1-based ``c.`` position of the (first) affected nucleotide, coding only.
    cds_position: int | None = None
    ref_aa: str | None = None
    alt_aa: str | None = None
    #: Position of the resulting premature termination codon, for truncating
    #: variants and frameshift classes. 1-based codon index.
    ptc_codon: int | None = None
    #: Free-form, deterministic extras used by class-level rows (exon spans for
    #: CNVs, representative indels for frameshift classes). Sorted keys only.
    attributes: dict[str, str] = Field(default_factory=dict)

    @field_validator("consequence_terms", mode="after")
    @classmethod
    def _primary_term_first(
        cls, value: tuple[Consequence, ...], info
    ) -> tuple[Consequence, ...]:
        primary = info.data.get("consequence")
        if value and primary is not None and value[0] != primary:
            raise ValueError("consequence_terms[0] must equal consequence")
        return value

    @property
    def variant_id(self) -> str:
        """Stable identifier: transcript plus HGVS ``c.`` description."""
        return f"{self.transcript}:{self.hgvs_c}"
