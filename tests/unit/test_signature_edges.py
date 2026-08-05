"""The edges of the evidence signature: what the digest can and cannot tell apart.

Oracle for this file: hand-built contexts whose expected canonical form is
written out in full, plus digest comparisons between two states that a reader
can see differ. Nothing here goes through the pipeline, so nothing here depends
on a shipped threshold or on the toy gene.

How these tests could nevertheless be wrong. The signature is only meaningful as
a *negative* claim -- "these two variants cannot evaluate differently" -- and a
digest comparison can only ever demonstrate the positive direction: that two
inputs the author believed different do in fact hash apart. A test written
around a pair the author never thought to try proves nothing about the pair that
actually collides. So the pairs below are chosen adversarially rather than
representatively: the sentinel against every falsy value that could stand in for
it, a bool against the int it compares equal to, a secondary consequence term
against its absence, and the sentinel against the ordinary list data that
happens to canonicalise to the same tuple. Two of those pairs do collide, and
they are recorded as characterisations, not as approvals.

The equality direction -- that two variants sharing a signature really do share
an evaluation -- is not testable from here at all. It is established by the
audit-mode tests in ``tests/layer3_classification/test_signature_reuse.py``,
which compare the cached run against the uncached one row for row.
"""

from __future__ import annotations

from datetime import date

import pytest

from vus_foresight.engine.context import MISSING, EvidenceContext
from vus_foresight.engine.signature import SignatureCache, canonical, evidence_signature
from vus_foresight.variant import Consequence, Variant, VariantKind

# --------------------------------------------------------------------------
# Builders. Deliberately thin: they set defaults, never behaviour.
# --------------------------------------------------------------------------


def _variant(*terms: Consequence) -> Variant:
    """A variant carrying nothing but the consequence terms under test."""
    chosen = terms or (Consequence.MISSENSE,)
    return Variant(
        gene="TOY",
        transcript="NM_000000.1",
        kind=VariantKind.SNV,
        hgvs_c="c.1A>G",
        consequence=chosen[0],
        consequence_terms=chosen,
    )


def _context(**namespaces: object) -> EvidenceContext:
    return EvidenceContext(data=dict(namespaces))


class _Stable:
    """A value outside the canonical vocabulary whose repr is deterministic."""

    def __repr__(self) -> str:
        return "Stable(7)"


class _AddressRepr:
    """Equal by value, but with no ``__repr__`` of its own."""

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _AddressRepr)

    def __hash__(self) -> int:
        return 0


# --------------------------------------------------------------------------
# canonical: the absent sentinel
# --------------------------------------------------------------------------


def test_the_absent_sentinel_canonicalises_to_a_control_character_token():
    """Absence has to survive canonicalisation as *something*, and that something
    must not be a value an adapter would plausibly emit."""
    reduced = canonical(MISSING)

    assert reduced == ("\x00missing",)


@pytest.mark.parametrize(
    ("falsy", "expected"),
    [
        (None, None),
        (False, False),
        (0, 0),
        (0.0, 0.0),
        ("", ""),
        ([], ()),
        ({}, ()),
    ],
)
def test_no_falsy_value_canonicalises_to_the_absent_sentinel(falsy, expected):
    """ "No frequency data" and "frequency is zero" are different evidence states.

    Collapsing them is how a system silently invents evidence, so every value
    that is falsy in Python must stay distinguishable from absence.
    """
    reduced = canonical(falsy)

    assert reduced == expected
    assert reduced != canonical(MISSING)


# --------------------------------------------------------------------------
# canonical: structure
# --------------------------------------------------------------------------


def test_mappings_are_ordered_by_key_at_every_depth():
    """Order independence has to be recursive, or a nested adapter payload built
    in a different insertion order would split one equivalence class in two."""
    reduced = canonical({"b": {"y": 1, "x": 2}, "a": 1})

    assert reduced == (("a", 1), ("b", (("x", 2), ("y", 1))))


def test_sequences_keep_their_order_inside_a_mapping():
    """A list is adapter output, not a set: ``region.tags`` arrives sorted, so a
    different order means a different adapter and is worth separating."""
    reduced = canonical({"tags": ["b", "a"]})

    assert reduced == (("tags", ("b", "a")),)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (date(2026, 1, 1), "datetime.date(2026, 1, 1)"),
        (frozenset({"tag"}), "frozenset({'tag'})"),
        (_Stable(), "Stable(7)"),
    ],
)
def test_a_value_outside_the_scalar_vocabulary_becomes_its_repr(value, expected):
    """The fallback is what keeps an unforeseen payload type from raising in the
    middle of a run; it is reached by anything that is not dict, sequence,
    scalar or None."""
    reduced = canonical(value)

    assert reduced == expected


def test_the_repr_fallback_lets_object_identity_reach_the_signature():
    """Characterisation of a latent determinism hazard, reported not fixed.

    An object without its own ``__repr__`` canonicalises to a string containing
    its memory address, so two values that compare equal reduce differently and
    the same input can digest differently between runs. No adapter emits such a
    value today -- payloads are dicts, lists and scalars parsed from TSV -- but
    the byte-for-byte guarantee rests on that remaining true.
    """
    first, second = _AddressRepr(), _AddressRepr()

    reduced_first, reduced_second = canonical(first), canonical(second)

    assert first == second
    assert reduced_first == f"<{__name__}._AddressRepr object at 0x{id(first):x}>"
    assert reduced_first != reduced_second


# --------------------------------------------------------------------------
# evidence_signature: absence versus presence
# --------------------------------------------------------------------------


@pytest.mark.parametrize("present", [0, 0.0, False, "", [], {}])
def test_an_absent_footprint_path_never_digests_like_a_present_falsy_value(present):
    """The distinction canonical() draws has to reach the digest itself; a
    variant with no frequency data must not share a verdict with one whose
    frequency is measured at zero."""
    footprint = ("frequency.gnomad.faf95_popmax",)
    absent = _context(other={"unused": 1})
    populated = _context(frequency={"gnomad": {"faf95_popmax": present}})

    absent_digest = evidence_signature(_variant(), absent, footprint)
    populated_digest = evidence_signature(_variant(), populated, footprint)

    assert absent_digest == "bdf5013faca99d507577a8b6bf7b6fc1"
    assert populated_digest != absent_digest


def test_an_explicit_null_is_indistinguishable_from_an_absent_path():
    """Characterisation. The expectation moved, not the code.

    An adapter that writes ``None`` is saying something different from an
    adapter that writes nothing, so the digest was expected to separate them.
    It does not, because ``EvidenceContext.get`` deliberately folds ``None``
    into MISSING before the signature ever sees it -- absent and null are one
    evidence state by decision upstream, and the digest inherits that decision
    rather than contradicting it.
    """
    footprint = ("frequency.gnomad.faf95_popmax",)
    null = _context(frequency={"gnomad": {"faf95_popmax": None}})
    absent = _context(frequency={"gnomad": {}})

    assert evidence_signature(_variant(), null, footprint) == evidence_signature(
        _variant(), absent, footprint
    )


def test_ordinary_list_data_can_forge_the_absent_sentinel():
    """Characterisation of a collision, reported not fixed.

    ``canonical`` reduces the sentinel to the one-element tuple
    ``("\\x00missing",)`` and reduces the list ``["\\x00missing"]`` to exactly
    the same tuple, so a footprint path holding that list digests as though the
    path were absent. This is the failure this module exists to prevent -- two
    different evidence states sharing one signature -- reachable only by an
    adapter emitting a NUL-prefixed string, which none does.
    """
    footprint = ("region.tags",)
    forged = _context(region={"tags": ["\x00missing"]})
    absent = _context(region={})

    assert evidence_signature(_variant(), forged, footprint) == evidence_signature(
        _variant(), absent, footprint
    )


# --------------------------------------------------------------------------
# evidence_signature: what enters the digest
# --------------------------------------------------------------------------


def test_a_path_outside_the_footprint_cannot_change_the_signature():
    """The bound the cache depends on, in the direction the audit cannot check.

    ``test_signature_reuse`` proves every path read is declared; this proves the
    converse mechanism -- an undeclared path is genuinely invisible to the
    digest, which is why an undeclared read would be unsound rather than merely
    untidy.
    """
    footprint = ("frequency.gnomad.faf95_popmax",)
    lean = _context(frequency={"gnomad": {"faf95_popmax": 0.001}})
    fat = _context(
        frequency={"gnomad": {"faf95_popmax": 0.001, "ac": 12}},
        predictor={"bayesdel": 0.9},
    )

    assert evidence_signature(_variant(), lean, footprint) == evidence_signature(
        _variant(), fat, footprint
    )


def test_the_same_value_reached_by_a_different_path_is_a_different_signature():
    """The path name is digested alongside its value, so a payload landing under
    the wrong key is visible rather than silently equivalent."""
    footprint = ("predictor.bayesdel", "splice.ds_max")
    under_first = _context(predictor={"bayesdel": 0.62}, splice={})
    under_second = _context(predictor={}, splice={"ds_max": 0.62})

    assert evidence_signature(_variant(), under_first, footprint) != evidence_signature(
        _variant(), under_second, footprint
    )


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (True, 1),
        (False, 0),
        (1, 1.0),
    ],
)
def test_values_that_compare_equal_in_python_still_digest_apart(left, right):
    """``canonical`` returns scalars untouched, so ``True`` and ``1`` reduce to
    values Python calls equal. The digest is taken over ``repr``, not over
    equality, which is what keeps a boolean flag from merging with a count."""
    footprint = ("predictor.flag",)

    left_digest = evidence_signature(_variant(), _context(predictor={"flag": left}), footprint)
    right_digest = evidence_signature(_variant(), _context(predictor={"flag": right}), footprint)

    assert left_digest != right_digest


def test_a_secondary_consequence_term_separates_two_variants_of_the_same_class():
    """A missense that also sits in the splice region carries evidence a plain
    missense does not, and the ``applies_to`` gate reads the whole term list --
    so merging them would apply a splice criterion to a variant it never saw."""
    footprint = ("frequency.gnomad.faf95_popmax",)
    context = _context(frequency={"gnomad": {"faf95_popmax": 0.0}})
    plain = _variant(Consequence.MISSENSE)
    also_splice = _variant(Consequence.MISSENSE, Consequence.SPLICE_REGION)

    assert evidence_signature(plain, context, footprint) != evidence_signature(
        also_splice, context, footprint
    )


def test_an_empty_footprint_reduces_every_variant_to_its_consequence_terms():
    """The degenerate case the shipped-spec test guards against, shown directly.

    With nothing declared, two variants differing only in evidence collapse into
    one signature -- which is why an empty ``field_footprint`` turns the cache
    from a saving into a silent merge.
    """
    poor = _context(frequency={"gnomad": {"faf95_popmax": 0.0}})
    rich = _context(frequency={"gnomad": {"faf95_popmax": 0.05}})

    assert evidence_signature(_variant(), poor, ()) == evidence_signature(_variant(), rich, ())
    assert evidence_signature(_variant(Consequence.NONSENSE), poor, ()) != evidence_signature(
        _variant(Consequence.MISSENSE), poor, ()
    )


def test_the_digest_is_a_stable_thirty_two_character_hex_string():
    """A pinned digest, because the recipe is part of the byte-for-byte promise.

    Changing what goes into the signature regroups every equivalence class in
    every shipped map. That is a decision, not an implementation detail, so it
    should have to be made deliberately here first.
    """
    footprint = ("frequency.gnomad.faf95_popmax",)
    context = _context(frequency={"gnomad": {"faf95_popmax": 0.001}})

    digest = evidence_signature(_variant(), context, footprint)

    assert digest == "6749bc52eaf863c09fcccc94e53b64c2"
    assert len(digest) == 32
    assert digest == evidence_signature(_variant(), context, footprint)


# --------------------------------------------------------------------------
# SignatureCache
# --------------------------------------------------------------------------


def test_a_fresh_cache_reports_a_zero_hit_rate_rather_than_dividing_by_zero():
    """The summary line is printed for every run, including one that enumerates
    nothing, so the empty case has to produce a number and not an exception."""
    cache = SignatureCache()

    assert (cache.size, cache.hits, cache.misses, cache.lookups) == (0, 0, 0, 0)
    assert cache.hit_rate == 0.0
    assert cache.summary() == "0 distinct evidence profiles across 0 variants (0.0% reused)"


def test_a_miss_and_a_hit_are_counted_separately():
    cache = SignatureCache()
    cache.put("aaaa", {"points": 4})

    absent = cache.get("bbbb")
    present = cache.get("aaaa")

    assert absent is None
    assert present == {"points": 4}
    assert (cache.hits, cache.misses, cache.lookups) == (1, 1, 2)
    assert cache.hit_rate == 0.5


def test_clearing_forgets_the_entries_and_the_counters_together():
    """A counter that survived the entries would describe a cache that no longer
    exists, and the hit rate is the run's only warning that the footprint has
    grown a per-variant field."""
    cache = SignatureCache()
    cache.put("aaaa", {"points": 4})
    cache.get("aaaa")
    cache.get("bbbb")

    cache.clear()

    assert (cache.size, cache.hits, cache.misses, cache.lookups) == (0, 0, 0, 0)
    assert cache.hit_rate == 0.0
    assert cache.get("aaaa") is None
    assert (cache.hits, cache.misses) == (0, 1)


def test_reputting_a_signature_replaces_the_entry_without_growing_the_cache():
    """``size`` is reported as the number of distinct evidence profiles, so it
    has to count signatures rather than writes."""
    cache = SignatureCache()

    cache.put("aaaa", {"points": 4})
    cache.put("aaaa", {"points": 6})

    assert cache.size == 1
    assert cache.get("aaaa") == {"points": 6}


def test_a_cached_none_is_recounted_as_a_miss_on_every_lookup():
    """Characterisation of a conflation, reported not fixed.

    ``get`` decides hit from miss by testing the stored value against ``None``,
    so a caller that memoises ``None`` gets no hit recorded and no value reused.
    The pipeline never stores ``None``, but the counters are the module's stated
    early warning, and here they would under-report reuse while ``size`` and
    ``misses`` drift apart -- an invariant the layer-3 suite asserts.
    """
    cache = SignatureCache()
    cache.put("aaaa", None)

    cache.get("aaaa")
    cache.get("aaaa")

    assert (cache.hits, cache.misses, cache.size) == (0, 2, 1)


def test_the_summary_groups_large_counts_for_a_human_reader():
    """This line is the operator's whole view of the run; an ungrouped six-digit
    count is what makes a collapsed hit rate easy to skim past."""
    cache = SignatureCache()
    for index in range(1000):
        cache.put(f"sig{index:04d}", index)
    for index in range(1000):
        cache.get(f"sig{index:04d}")

    summary = cache.summary()

    assert summary == "1,000 distinct evidence profiles across 1,000 variants (100.0% reused)"


def test_a_mistyped_counter_cannot_be_assigned_on_the_cache():
    """``__slots__`` is what stops ``cache.hit = 0`` from creating a second,
    silently ignored counter beside ``cache.hits``."""
    cache = SignatureCache()

    with pytest.raises(AttributeError, match=r"'SignatureCache' object has no attribute 'hit'"):
        cache.hit = 1  # type: ignore[attr-defined]
