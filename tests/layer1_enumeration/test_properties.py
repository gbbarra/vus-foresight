"""Layer 1, property-based.

Hypothesis is used where the invariant holds for *every* position rather than
for a chosen few, which is precisely where hand-written cases give false
confidence.
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from vus_foresight.annotate.annotator import CodingEdit, annotate_coding_edits
from vus_foresight.genome.sequence import alternatives
from vus_foresight.testing import build_synthetic_gene

GENE = build_synthetic_gene(strand="-")
TRANSCRIPT = GENE.transcript


@settings(max_examples=200, deadline=None)
@given(st.integers(min_value=1, max_value=TRANSCRIPT.cds_length))
def test_alternatives_at_any_cds_position_are_the_other_three_bases(cds_position):
    ref = TRANSCRIPT.cds[cds_position - 1]
    produced = {
        annotate_coding_edits(
            TRANSCRIPT, (CodingEdit(cds_position=cds_position, ref=ref, alt=alt),)
        ).hgvs_c.split(">")[1]
        for alt in alternatives(ref)
    }
    assert produced == set("ACGT") - {ref}


@settings(max_examples=200, deadline=None)
@given(st.integers(min_value=1, max_value=TRANSCRIPT.cds_length))
def test_every_cds_position_round_trips_through_the_transcript(cds_position):
    tx = TRANSCRIPT.cds_position_to_tx(cds_position)
    genomic = TRANSCRIPT.genomic_at(tx)
    assert TRANSCRIPT.tx_at_genomic(genomic) == tx
    assert TRANSCRIPT.c_at_tx(tx).base == cds_position


@settings(max_examples=100, deadline=None)
@given(st.integers(min_value=1, max_value=TRANSCRIPT.n_codons))
def test_codon_bounds_are_contiguous_and_in_range(codon_index):
    first, middle, last = TRANSCRIPT.codon_bounds_cds(codon_index)
    assert (first, middle, last) == (first, first + 1, first + 2)
    assert 1 <= first and last <= TRANSCRIPT.cds_length
    assert TRANSCRIPT.codon_sequence(codon_index) == TRANSCRIPT.cds[first - 1 : last]
