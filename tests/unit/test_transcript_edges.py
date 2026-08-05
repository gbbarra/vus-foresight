"""Transcript geometry in the places that only exist between the exons.

``test_transcript_contract.py`` next door asks what the model refuses to be
*built* from. This file asks what a well-formed transcript answers once it
exists, in the three regions where an answer can be internally consistent and
still wrong: the introns, the exonic bases that abut them, and the NMD
boundary derived from the penultimate exon.

How these tests could be wrong -- the oracle, stated so it can be attacked:

* **Independent arithmetic.** Expected coordinates are derived from the
  fixture's *declared* geometry (``TOY_EXON_LENGTHS``, ``TOY_INTRON_LENGTHS``,
  ``TOY_UTR5``) or hard-coded, never recomputed with the same expression the
  module under test uses. A test that re-derives ``penultimate_end - 50 + 1``
  from ``transcript._exon_tx_bounds`` agrees with the code by construction and
  would survive the formula being wrong; the numbers here would not.
* **Strand symmetry.** The plus and minus fixtures are the same gene reflected.
  Every transcript-space answer must be *identical* between them and every
  genomic answer must *mirror*. A reverse-complement bug breaks the minus case
  only, which is the one control this project cannot afford to lose.

Transcripts built locally by ``_transcript`` carry a filler sequence: every
behaviour exercised here is coordinate arithmetic, and none of it reads a base.
"""

from __future__ import annotations

import pytest

from vus_foresight.genome.transcript import (
    CPosition,
    Exon,
    Transcript,
    format_c_position,
)
from vus_foresight.testing import TOY_EXON_LENGTHS, TOY_UTR5

STRANDS = ("+", "-")


def _transcript(
    *,
    strand: str,
    exon_lengths: tuple[int, ...],
    intron_lengths: tuple[int, ...],
    cds_start_tx: int,
    cds_end_tx: int,
) -> Transcript:
    """A minimal transcript with the requested exon/intron geometry.

    Laid out from a fixed anchor and walked in *transcript* direction, so the
    minus-strand version genuinely descends rather than being a plus-strand
    layout with a flag changed.
    """
    total = sum(exon_lengths)
    span = total + sum(intron_lengths)
    cursor = 1_000_000 if strand == "+" else 1_000_000 + span
    intervals: list[tuple[int, int]] = []
    for i, length in enumerate(exon_lengths):
        gap = intron_lengths[i] if i < len(intron_lengths) else 0
        if strand == "+":
            intervals.append((cursor, cursor + length - 1))
            cursor += length + gap
        else:
            intervals.append((cursor - length + 1, cursor))
            cursor -= length + gap
    return Transcript(
        transcript_id="NM_000001.1",
        gene="EDGE",
        chrom="chr1",
        strand=strand,
        exons=tuple(Exon(label=str(i + 1), start=s, end=e) for i, (s, e) in enumerate(intervals)),
        cds_start_tx=cds_start_tx,
        cds_end_tx=cds_end_tx,
        sequence="A" * total,
    )


# --------------------------------------------------------------------------
# What a c. position says about itself.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("position", "expected"),
    [
        (CPosition(base=100), (False, True)),
        (CPosition(base=100, offset=3), (True, False)),
        (CPosition(base=101, offset=-2), (True, False)),
        (CPosition(base=-20), (False, False)),
        (CPosition(base=-20, offset=1), (True, False)),
        (CPosition(base=12, utr3=True), (False, False)),
        (CPosition(base=12, offset=-1, utr3=True), (True, False)),
    ],
    ids=["cds", "donor-side", "acceptor-side", "utr5", "utr5-intronic", "utr3", "utr3-intronic"],
)
def test_a_c_position_reports_intronic_and_coding_independently(position, expected):
    """These two flags gate which criteria a variant is even eligible for.

    The pair matters more than either flag: ``c.-20+1`` is intronic and not
    coding, and code that inferred one from the other -- treating any non-coding
    position as intronic, or any offset-free one as coding -- would misroute the
    5' UTR and the 3' UTR into the coding branch of the spec.
    """
    result = (position.is_intronic, position.is_coding)

    assert result == expected


# --------------------------------------------------------------------------
# Asking a question the spliced transcript cannot answer.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("offset", [-1, 1], ids=["before-the-start", "past-the-end"])
def test_asking_which_exon_holds_a_position_off_the_transcript_is_refused(minus_gene, offset):
    """Unlike ``genomic_at``, this entry point has no up-front bounds check.

    It falls through the exon loop instead, so the refusal has to come from the
    loop exhausting. If it ever returned a clamped index instead, a position off
    the end would be attributed to the last exon and reported under that exon's
    clinical label.
    """
    transcript = minus_gene.transcript
    tx_pos = 0 if offset < 0 else transcript.length + 1

    with pytest.raises(ValueError, match=rf"transcript position {tx_pos} is out of range"):
        transcript.exon_index_at_tx(tx_pos)


def test_the_first_and_last_transcript_positions_name_the_first_and_last_exons(minus_gene):
    """The legal side of the boundary above, and the index is transcript-order.

    On the minus strand the first exon in transcript order is the genomically
    *last* one, so an implementation that indexed by genomic position would
    return 4 here instead of 0.
    """
    transcript = minus_gene.transcript

    first = transcript.exon_index_at_tx(1)
    last = transcript.exon_index_at_tx(transcript.length)

    assert (first, last) == (0, len(transcript.exons) - 1)


@pytest.mark.parametrize(
    ("position", "rendered"),
    [
        (CPosition(base=100, offset=3), r"c\.100\+3"),
        (CPosition(base=101, offset=-2), r"c\.101-2"),
        (CPosition(base=12, offset=-1, utr3=True), r"c\.\*12-1"),
    ],
    ids=["donor-side", "acceptor-side", "utr3-acceptor-side"],
)
def test_an_intronic_c_position_has_no_transcript_position(minus_gene, position, rendered):
    """The transcript sequence is spliced mRNA: intronic bases are not in it.

    Silently dropping the offset would map ``c.100+3`` onto the exonic base
    ``c.100``, three bases away from where the variant actually is and inside a
    different consequence class.
    """
    with pytest.raises(ValueError, match=rf"{rendered} is intronic, not exonic"):
        minus_gene.transcript.tx_at_c(position)


# --------------------------------------------------------------------------
# CDS offset straight to the genome.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("strand", STRANDS)
def test_the_first_cds_base_sits_one_utr_length_into_the_first_exon(both_strands, strand):
    """The shortcut from CDS space to the genome, checked against the 5' UTR.

    ``TOY_UTR5`` is the fixture's *declared* untranslated length, so the
    expectation here is independent of ``cds_start_tx`` and of the module's own
    exon arithmetic. The initiation codon must then read forward on the plus
    strand and backward on the minus one -- the single place where strand
    becomes arithmetic, seen from the CDS side.
    """
    transcript = both_strands[strand].transcript
    first_exon = transcript.exons[0]
    expected_a_of_atg = first_exon.start + TOY_UTR5 if strand == "+" else first_exon.end - TOY_UTR5

    codon_one = [transcript.genomic_at_cds(p) for p in (1, 2, 3)]

    step = 1 if strand == "+" else -1
    assert codon_one == [expected_a_of_atg + step * i for i in range(3)]


# --------------------------------------------------------------------------
# Introns.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("strand", STRANDS)
def test_two_touching_exons_cannot_describe_the_intron_between_them(strand):
    """Adjacency passes validation but leaves nothing to number.

    The contract test asserts that abutting exons are legal -- overlap is the
    thing refused -- so the impossible case survives construction and only
    surfaces here, on first access to ``introns``. A zero-length intron would
    otherwise make ``c.100+1`` and ``c.101-1`` name the same nonexistent base.
    """
    transcript = _transcript(
        strand=strand,
        exon_lengths=(60, 40),
        intron_lengths=(0,),
        cds_start_tx=1,
        cds_end_tx=99,
    )

    with pytest.raises(ValueError, match=r"exons 1/2 are contiguous"):
        _ = transcript.introns


@pytest.mark.parametrize("strand", STRANDS)
@pytest.mark.parametrize(
    "method", ["intron_c_position", "intron_genomic"], ids=["c-position", "genomic"]
)
@pytest.mark.parametrize("n", [0, "past-the-end"], ids=["zero", "past-the-end"])
def test_an_intronic_offset_outside_the_intron_is_refused(both_strands, strand, method, n):
    """Intronic offsets are 1-based from the donor; there is no ``+0``.

    Both entry points guard independently, and both must: ``intron_genomic``
    would happily return ``upstream_last_genomic`` for ``n=0`` -- an exonic base
    reported as intronic -- and one past the length would return the first base
    of the next exon.
    """
    transcript = both_strands[strand].transcript
    intron = transcript.introns[0]
    offset = 0 if n == 0 else intron.length + 1

    with pytest.raises(ValueError, match=rf"intronic offset {offset} outside 1\.\.{intron.length}"):
        getattr(transcript, method)(intron, offset)


@pytest.mark.parametrize("strand", STRANDS)
def test_an_intron_is_numbered_from_the_c_anchors_of_the_exons_it_flanks(both_strands, strand):
    """The anchor, not just the offset, has to be right.

    The existing midpoint test checks only ``.offset``, so a wrong anchor would
    render ``c.999+1`` for the first intronic base and still pass it. Here the
    whole HGVS string is pinned: intron 1 follows CDS base 30 (the first exon is
    60 nt and 30 of them are 5' UTR) and intron 2 follows CDS base 70, and both
    numbering halves must name the exon on their own side.
    """
    transcript = both_strands[strand].transcript
    first, second = transcript.introns[0], transcript.introns[1]

    rendered = {
        "first-donor": format_c_position(transcript.intron_c_position(first, 1)),
        "first-last-of-5prime-half": format_c_position(transcript.intron_c_position(first, 60)),
        "first-first-of-3prime-half": format_c_position(transcript.intron_c_position(first, 61)),
        "first-acceptor": format_c_position(transcript.intron_c_position(first, first.length)),
        "odd-central": format_c_position(transcript.intron_c_position(second, 16)),
        "odd-just-past-central": format_c_position(transcript.intron_c_position(second, 17)),
    }

    assert rendered == {
        "first-donor": "30+1",
        "first-last-of-5prime-half": "30+60",
        "first-first-of-3prime-half": "31-60",
        "first-acceptor": "31-1",
        "odd-central": "70+16",
        "odd-just-past-central": "71-15",
    }


@pytest.mark.parametrize("strand", STRANDS)
def test_the_nth_intronic_base_walks_the_genome_in_transcript_direction(both_strands, strand):
    """Intronic offsets count from the donor, which is strand-relative.

    On the minus strand the first intronic base is genomically *below* the last
    exonic one. Counting up on both strands would place every deep-intronic and
    every splice-site variant on the wrong side of the junction -- and the c.
    string would still look perfectly well formed.
    """
    transcript = both_strands[strand].transcript
    intron = transcript.introns[1]
    step = 1 if strand == "+" else -1

    walked = [transcript.intron_genomic(intron, n) for n in (1, 2, intron.length)]

    assert walked == [
        intron.upstream_last_genomic + step,
        intron.upstream_last_genomic + 2 * step,
        intron.downstream_first_genomic - step,
    ]


# --------------------------------------------------------------------------
# Where NMD stops.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("strand", STRANDS)
def test_the_toy_nmd_boundary_is_the_84th_cds_base(both_strands, strand):
    """The 50-nt rule, counted from the fixture's declared exon lengths.

    Exons of 60/40/30/33 nt put the end of the penultimate exon at transcript
    position 163; the last 50 nt of it therefore start at 114, and with a 30 nt
    5' UTR that is CDS base 84. Codon 28 spans CDS 82-84, so it *begins* before
    the boundary and is degraded; codon 29 is the first that escapes. This is
    the number that decides whether PVS1 applies at full strength, and it is
    hard-coded here on purpose: rederiving it from ``_exon_tx_bounds`` with the
    module's own expression would agree with a wrong formula.
    """
    transcript = both_strands[strand].transcript

    boundary = transcript.nmd_escape_cds_start
    degraded, escaping = transcript.ptc_escapes_nmd(28), transcript.ptc_escapes_nmd(29)

    assert (boundary, degraded, escaping) == (84, False, True)


def test_the_toy_boundary_matches_the_declared_exon_and_utr_lengths():
    """Guards the constants the test above hard-codes.

    If the fixture's geometry is ever retuned, this fails and says so, instead
    of leaving 84 silently stale and the NMD test passing for the wrong reason.
    """
    penultimate_end_tx = sum(TOY_EXON_LENGTHS[:-1])

    boundary_cds = penultimate_end_tx - 50 + 1 - TOY_UTR5

    assert (penultimate_end_tx, boundary_cds) == (163, 84)


@pytest.mark.parametrize("strand", STRANDS)
def test_a_single_exon_transcript_degrades_no_premature_terminator(strand):
    """No exon junction means no exon junction complex, so nothing is degraded.

    NMD is EJC-dependent: a lone exon has no downstream junction to mark, and
    every PTC in it escapes. Returning the usual penultimate-exon boundary here
    would index ``_exon_tx_bounds[-2]``, wrap around to the only exon, and
    declare the first two thirds of the gene NMD-competent -- upgrading PVS1 on
    a transcript where the whole mechanism is absent.
    """
    transcript = _transcript(
        strand=strand,
        exon_lengths=(120,),
        intron_lengths=(),
        cds_start_tx=1,
        cds_end_tx=120,
    )

    escapes = [transcript.ptc_escapes_nmd(c) for c in (1, 2, transcript.n_codons)]

    assert (transcript.nmd_escape_cds_start, escapes) == (1, [True, True, True])


@pytest.mark.parametrize("strand", STRANDS)
def test_a_penultimate_exon_ending_within_50_nt_of_the_start_clamps_the_boundary(strand):
    """The 50-nt window can start before the CDS does, and must not go negative.

    With a 10 nt first exon the window opens at transcript position -39. Left
    unclamped, ``nmd_escape_cds_start`` would be negative, every codon would
    compare greater than it, and the answer would come out right by accident --
    until some caller used the value as a CDS offset. Clamping to 1 says the
    same thing in coordinates that exist.
    """
    transcript = _transcript(
        strand=strand,
        exon_lengths=(10, 30),
        intron_lengths=(80,),
        cds_start_tx=1,
        cds_end_tx=39,
    )

    escapes_everywhere = [transcript.ptc_escapes_nmd(c) for c in (1, transcript.n_codons)]

    assert (transcript.nmd_escape_cds_start, escapes_everywhere) == (1, [True, True])


# --------------------------------------------------------------------------
# Exonic bases adjacent to a junction.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("strand", STRANDS)
def test_a_two_base_exon_contributes_only_its_own_two_bases_to_the_splice_region(strand):
    """An exon shorter than the 3 nt window must not spill into its neighbours.

    The window is clamped to the exon on both sides. Without the clamp a 2 nt
    exon's acceptor window would reach one base into the following exon and its
    donor window one base back into the preceding one, quietly widening the
    splice region across a junction the reader cannot see. The union is asserted
    exactly, not by membership, so a phantom position anywhere fails.
    """
    transcript = _transcript(
        strand=strand,
        exon_lengths=(60, 2, 40),
        intron_lengths=(100, 100),
        cds_start_tx=1,
        cds_end_tx=102,
    )

    region = transcript.exonic_splice_region_tx

    # 58-60 donor of exon 1, 61-62 the whole of exon 2, 63-65 acceptor of exon 3.
    assert region == frozenset({58, 59, 60, 61, 62, 63, 64, 65})
