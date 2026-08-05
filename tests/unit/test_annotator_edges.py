"""Edge cases in the coding annotator: its input guards and its severity ranking.

Two oracles, because the two halves of this file fail in two different ways.

* The guards at the top of :func:`annotate_coding_edits` are a *contract*, and
  the oracle is that contract as the function's own docstring states it: the
  edits must be non-empty, must not repeat a CDS position, must all lie in one
  codon, and each ``alt`` must differ from its ``ref``. Every one of those
  describes a variant that cannot exist, so an edit reaching the annotator in
  that shape means the enumerator upstream is broken. The expectations below
  therefore pin the *message* and not merely the exception type: a guard that
  fires quoting the wrong coordinate sends whoever reads the traceback to the
  wrong place in the enumerator, which is worse than no coordinate at all. The
  way these tests can be wrong is that the contract was transcribed wrongly
  here -- something a second reader checks against the docstring in a minute.

* :func:`most_severe` and ``CONSEQUENCE_RANK`` are ranked to agree with Ensembl
  VEP, so that "primary consequence" means the same thing here as in the
  external annotation this project is compared against. The oracle is that
  published ordering, and each expectation below is read off it (nonsense
  outranks missense, a canonical splice term outranks a nonsense one, a splice
  region term outranks a synonymous one) rather than off the dict under test --
  otherwise the test would only prove the dict equals itself. One pair does not
  agree with VEP, and that disagreement is recorded rather than smoothed over:
  see ``test_an_intronic_term_currently_outranks_a_utr_term``.

Not repeated here, because ``tests/layer2_annotation`` already pins it: the
whole-codon translation rule, the minimal ``delins`` span, and the strand
symmetry of the genomic coordinate. The guards are pure CDS-offset arithmetic
and carry no strand reasoning at all, so they are exercised on one strand only.
"""

from __future__ import annotations

import pytest

from vus_foresight.annotate.annotator import (
    CONSEQUENCE_RANK,
    CodingEdit,
    annotate_coding_edits,
    most_severe,
)
from vus_foresight.variant import Consequence

#: CDS positions 4, 5, 6 of the synthetic gene are the ``GCC`` codon that
#: ``build_synthetic_gene`` plants on purpose (codon index 2). Hard-coding it
#: keeps the expected guard messages literal instead of reconstructing them
#: from the same expression the code under test uses.
GCC_CODON_FIRST_CDS = 4


def _edit_matching_reference(transcript, cds_position: int) -> CodingEdit:
    """A well-formed edit at ``cds_position``: correct ref, some different alt.

    Used where the point of the test is a *different* guard, so this edit must
    be innocent of everything except the thing under test.
    """
    ref = transcript.cds[cds_position - 1]
    alt = next(base for base in "ACGT" if base != ref)
    return CodingEdit(cds_position=cds_position, ref=ref, alt=alt)


# ----------------------------------------------------------------------
# severity ranking
# ----------------------------------------------------------------------


def test_every_consequence_term_has_exactly_one_severity_rank():
    """A missing or shared rank makes the primary term non-deterministic.

    ``most_severe`` is ``min`` over the ranks, and ``min`` breaks a tie by
    input order -- so two terms sharing a rank would make the emitted primary
    consequence depend on the order a caller happened to list them in, and this
    project guarantees byte-identical output for identical input. A term
    missing from the dict is worse still: a ``KeyError`` at annotation time.
    """
    ranks = sorted(CONSEQUENCE_RANK.values())

    assert set(CONSEQUENCE_RANK) == set(Consequence)
    assert ranks == list(range(len(Consequence)))


@pytest.mark.parametrize(
    ("terms", "expected"),
    [
        ([Consequence.MISSENSE, Consequence.NONSENSE], Consequence.NONSENSE),
        ([Consequence.SYNONYMOUS, Consequence.SPLICE_REGION], Consequence.SPLICE_REGION),
        ([Consequence.MISSENSE, Consequence.SPLICE_REGION], Consequence.MISSENSE),
        ([Consequence.NONSENSE, Consequence.SPLICE_DONOR], Consequence.SPLICE_DONOR),
        ([Consequence.FRAMESHIFT, Consequence.STOP_LOST], Consequence.FRAMESHIFT),
        ([Consequence.EXON_DELETION, Consequence.SPLICE_ACCEPTOR], Consequence.EXON_DELETION),
        ([Consequence.UTR3, Consequence.UTR5], Consequence.UTR5),
        ([Consequence.STOP_RETAINED], Consequence.STOP_RETAINED),
    ],
    ids=[
        "nonsense-beats-missense",
        "splice-region-beats-synonymous",
        "missense-beats-splice-region",
        "donor-beats-nonsense",
        "frameshift-beats-stop-lost",
        "exon-deletion-beats-acceptor",
        "utr5-beats-utr3",
        "single-term",
    ],
)
def test_the_primary_term_of_a_set_is_the_one_vep_ranks_most_severe(terms, expected):
    """Downstream, a criterion fires on the primary term alone.

    Picking the weaker of two terms would silently drop the evidence gate the
    stronger one opens -- e.g. a donor variant annotated as a nonsense one is
    scored by the wrong branch of the spec.
    """
    result = most_severe(terms)

    assert result == expected


def test_an_intronic_term_currently_outranks_a_utr_term():
    """Characterisation of a deviation from the oracle, not an endorsement of it.

    My expectation here was ``UTR5``, read off Ensembl VEP's published table,
    where ``5_prime_UTR_variant`` and ``3_prime_UTR_variant`` both outrank
    ``intron_variant``. ``CONSEQUENCE_RANK`` inverts that pair, so the
    expectation is what changed -- and the reason it changed is that the pair is
    unreachable today, not that VEP was read wrongly: no code path in the
    package ever puts ``UTR5`` or ``UTR3`` on a variant, so no emitted map can
    contain a term set where this ordering decides anything. It is pinned here
    so that whoever adds UTR enumeration sees the discrepancy as a failing test
    instead of shipping a primary consequence that disagrees with the external
    annotation this project validates itself against.
    """
    result = most_severe([Consequence.INTRONIC, Consequence.UTR5, Consequence.UTR3])

    assert result == Consequence.INTRONIC


def test_the_primary_term_does_not_depend_on_the_order_the_terms_are_listed_in():
    """Order independence is what makes the emitted map reproducible.

    The caller assembles terms in a ``set``, whose iteration order is not part
    of the contract; if the answer moved with it, two runs on identical input
    could disagree and the byte-comparison test would fail intermittently.
    """
    terms = [Consequence.SPLICE_REGION, Consequence.NONSENSE, Consequence.SYNONYMOUS]

    forward = most_severe(terms)
    backward = most_severe(list(reversed(terms)))

    assert (forward, backward) == (Consequence.NONSENSE, Consequence.NONSENSE)


def test_the_exported_ranking_agrees_with_the_term_order_the_annotator_emits(plus_gene):
    """``most_severe`` is exported but unused inside the package, so it can drift.

    An outside consumer calling it on ``consequence_terms`` must land on the
    same primary term the annotator itself chose; if the two orderings ever
    diverged, the published map and its own helper would disagree. c.28 is the
    first base of exon 2, so its codon carries a splice-region term alongside
    the coding one and there is genuinely something to rank.
    """
    transcript = plus_gene.transcript

    variant = annotate_coding_edits(transcript, (CodingEdit(cds_position=28, ref="C", alt="T"),))

    assert variant.hgvs_c == "c.28C>T"
    assert variant.hgvs_p == "p.Arg10Cys"
    assert variant.consequence_terms == (Consequence.MISSENSE, Consequence.SPLICE_REGION)
    assert most_severe(list(variant.consequence_terms)) == variant.consequence


# ----------------------------------------------------------------------
# input guards
# ----------------------------------------------------------------------


def test_an_empty_edit_tuple_is_refused_instead_of_annotated(plus_gene):
    """With no edits there is no codon to translate and no span to report.

    Annotating it anyway would emit a row describing the reference allele as if
    it were a variant, which is a fabricated entry in the map.
    """
    transcript = plus_gene.transcript

    with pytest.raises(ValueError, match="at least one edit is required"):
        annotate_coding_edits(transcript, ())


def test_repeating_a_cds_position_is_refused_rather_than_silently_resolved(plus_gene):
    """Two alts at one base is not a variant, it is two variants.

    Both edits below are individually well-formed -- correct ref, real alt --
    so the repetition is the only thing wrong. Without the guard the second
    edit would just overwrite the first in the codon buffer and the annotator
    would emit one of the two, chosen by sort order, with no trace of the other.
    """
    transcript = plus_gene.transcript
    edits = (
        CodingEdit(cds_position=GCC_CODON_FIRST_CDS, ref="G", alt="A"),
        CodingEdit(cds_position=GCC_CODON_FIRST_CDS, ref="G", alt="T"),
    )

    with pytest.raises(ValueError, match=r"duplicate CDS positions in edits: \[4, 4\]"):
        annotate_coding_edits(transcript, edits)


@pytest.mark.parametrize(
    ("positions", "expected"),
    [
        ((3, 4), r"edits span 2 codons \(\[1, 2\]\)"),
        ((1, 4, 7), r"edits span 3 codons \(\[1, 2, 3\]\)"),
    ],
    ids=["straddles-one-codon-boundary", "three-codons"],
)
def test_edits_that_leave_a_single_codon_are_refused(plus_gene, positions, expected):
    """This function's whole contract is "one codon, translated as a unit".

    Given positions in two codons it would pick one codon index out of the set
    and annotate that codon only, producing a ``c.`` description of part of the
    change and a ``p.`` description of a protein that no allele encodes. The
    message names every codon involved so the caller can see which of its
    positions disagreed.
    """
    transcript = plus_gene.transcript
    edits = tuple(_edit_matching_reference(transcript, p) for p in positions)

    with pytest.raises(ValueError, match=expected):
        annotate_coding_edits(transcript, edits)


@pytest.mark.parametrize(
    ("declared_ref", "expected"),
    [
        ("A", r"reference mismatch at CDS 4: config says 'G', edit says 'A'"),
        ("C", r"reference mismatch at CDS 4: config says 'G', edit says 'C'"),
        ("T", r"reference mismatch at CDS 4: config says 'G', edit says 'T'"),
    ],
    ids=["ref-A", "ref-C", "ref-T"],
)
def test_an_edit_whose_reference_disagrees_with_the_transcript_is_refused(
    plus_gene, declared_ref, expected
):
    """A wrong ref means the caller and the transcript disagree about the locus.

    The annotator builds the alt codon from the *transcript* sequence, so it
    would happily produce a well-formed ``c.`` and ``p.`` for a codon the caller
    never meant -- the exact "looks right, is wrong" failure this layer exists
    to catch. The message quotes both bases so the disagreement is visible
    without opening the reference.
    """
    transcript = plus_gene.transcript
    edit = CodingEdit(cds_position=GCC_CODON_FIRST_CDS, ref=declared_ref, alt="A")

    with pytest.raises(ValueError, match=expected):
        annotate_coding_edits(transcript, (edit,))


@pytest.mark.parametrize(
    ("edits", "expected"),
    [
        (
            (CodingEdit(cds_position=4, ref="G", alt="G"),),
            r"no-op edit at CDS 4: G>G",
        ),
        (
            (
                CodingEdit(cds_position=4, ref="G", alt="A"),
                CodingEdit(cds_position=5, ref="C", alt="C"),
            ),
            r"no-op edit at CDS 5: C>C",
        ),
    ],
    ids=["only-edit-is-a-no-op", "no-op-riding-along-with-a-real-edit"],
)
def test_an_edit_that_changes_no_base_is_refused(plus_gene, edits, expected):
    """A no-op edit is a bookkeeping error upstream, not a variant.

    The lone case would leave the alt codon identical to the reference, so the
    minimal-differing-span computation would index an empty list and die with
    an opaque ``IndexError`` far from the cause. The second case is the one the
    guard has to be per-edit to catch: the codon *does* change overall, so a
    check phrased on the assembled codon would let it through and the emitted
    ``mutational_distance`` would count a base that never moved.
    """
    transcript = plus_gene.transcript

    with pytest.raises(ValueError, match=expected):
        annotate_coding_edits(transcript, edits)
