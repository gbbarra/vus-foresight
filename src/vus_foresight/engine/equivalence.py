"""Equivalence classes (spec section 5).

Do not evaluate forty-seven thousand variants independently when their intrinsic
evidence profile is identical. Every nonsense variant in the same PVS1 zone
carries the same evidence; so does every frameshift whose PTC lands in that
zone; so does every synonymous change with no predicted splice effect.

The class key is derived from the evaluation trace itself rather than from a
hand-written table of "things that are probably the same". Two variants are in
the same class exactly when the engine tested the same criteria and got the same
answers -- which is the definition, not an approximation of it.

Missense never collapses. Its evidence is per-variant (PP3/BP4, PS3/BS3), so
two missense variants that look identical today diverge the moment a predictor
score or an assay result lands. The list of never-collapsing consequences is
specification data, not a constant here.
"""

from __future__ import annotations

from hashlib import blake2b

from ..acmg import CriterionOutcome
from ..variant import Variant
from .evaluator import Evaluation
from .spec import VCEPSpec

__all__ = ["equivalence_class_id", "equivalence_key"]


def equivalence_key(variant: Variant, evaluation: Evaluation, spec: VCEPSpec) -> str:
    """Human-readable class key. The id is a digest of this string."""
    terms = set(variant.consequence_terms) or {variant.consequence}
    never_collapse = terms.intersection(spec.equivalence.never_collapse)
    identity = variant.variant_id if never_collapse else ""

    applied = ";".join(f"{c.code}:{c.strength.value}" for c in evaluation.applied)
    unevaluable = ";".join(
        s.code for s in evaluation.skipped if s.outcome is CriterionOutcome.NOT_EVALUABLE
    )
    return "|".join(
        [
            variant.gene,
            variant.transcript,
            spec.spec_version,
            variant.consequence.value,
            identity,
            applied,
            unevaluable,
        ]
    )


def equivalence_class_id(variant: Variant, evaluation: Evaluation, spec: VCEPSpec) -> str:
    """Stable, short identifier for the class.

    A digest rather than a counter: it must be reproducible across runs, across
    machines and across partial runs of a single gene, which a counter is not.
    """
    digest = blake2b(
        equivalence_key(variant, evaluation, spec).encode("utf-8"), digest_size=8
    ).hexdigest()
    return f"ec_{digest}"
