"""ACMG/AMP vocabulary and the Tavtigian point system (ClinGen SVI).

This module is deliberately free of any gene- or disease-specific knowledge. The
numeric values of the point system and of the classification thresholds are
*defaults* only: a VCEP specification may override every one of them (see
:class:`PointSystem`, loaded from YAML by :mod:`vus_foresight.engine.spec`).

Design note (spec section 2): combining criteria is addition and a gap is
subtraction. Everything downstream -- ceilings, gaps, minimum sufficient
evidence sets -- depends on that being literally true, so the point system is
the single place where strength is turned into a number.
"""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

__all__ = [
    "Direction",
    "Strength",
    "ACMGClass",
    "EvidenceClass",
    "FeasibilityTag",
    "BlockingReason",
    "CriterionOutcome",
    "SkipReason",
    "PointSystem",
    "DEFAULT_POINT_SYSTEM",
]


class Direction(str, Enum):
    """Which way a criterion pushes."""

    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class Strength(str, Enum):
    """Evidence strength labels used by ACMG/AMP.

    ``STAND_ALONE`` is benign-only (BA1); ``VERY_STRONG`` is pathogenic-only in
    the 2015 vocabulary, but the point system treats a benign very-strong as the
    same magnitude as stand-alone, so both are representable.
    """

    SUPPORTING = "supporting"
    MODERATE = "moderate"
    STRONG = "strong"
    VERY_STRONG = "very_strong"
    STAND_ALONE = "stand_alone"


#: Canonical ordering, weakest first. Used for "max plausible strength" logic.
STRENGTH_ORDER: tuple[Strength, ...] = (
    Strength.SUPPORTING,
    Strength.MODERATE,
    Strength.STRONG,
    Strength.VERY_STRONG,
    Strength.STAND_ALONE,
)


class ACMGClass(str, Enum):
    """The five-tier terminology."""

    PATHOGENIC = "P"
    LIKELY_PATHOGENIC = "LP"
    UNCERTAIN = "VUS"
    LIKELY_BENIGN = "LB"
    BENIGN = "B"


class EvidenceClass(str, Enum):
    """Spec section 4 -- the backbone distinction of the whole project.

    ``INTRINSIC``      computable with no patient anywhere.
    ``SEMI_INTRINSIC`` depends on the state of a public database on a date.
    ``EXTRINSIC``      requires observation in a patient, family or cohort.
    """

    INTRINSIC = "intrinsic"
    SEMI_INTRINSIC = "semi_intrinsic"
    EXTRINSIC = "extrinsic"


class FeasibilityTag(str, Enum):
    """How hard it would be to actually acquire a missing piece of evidence.

    ``AVAILABLE_UNINGESTED`` is the headline finding of the whole system: the
    data exists publicly and simply has not been wired in yet. That is
    engineering work, not bench work.
    """

    AVAILABLE_UNINGESTED = "available_uningested"
    ASSAY_FEASIBLE = "assay_feasible"
    NEEDS_FAMILIES = "needs_families"
    NEEDS_CASES = "needs_cases"
    INTRACTABLE = "intractable"


#: Acquisition cost ordering, cheapest first. Drives the ranking of the minimum
#: sufficient sets that the report surfaces.
FEASIBILITY_COST: dict[FeasibilityTag, int] = {
    FeasibilityTag.AVAILABLE_UNINGESTED: 0,
    FeasibilityTag.ASSAY_FEASIBLE: 1,
    FeasibilityTag.NEEDS_CASES: 2,
    FeasibilityTag.NEEDS_FAMILIES: 3,
    FeasibilityTag.INTRACTABLE: 99,
}


class BlockingReason(str, Enum):
    """Why a variant sits at VUS -- categorised so the map is ``GROUP BY``-able.

    The question "how many BRCA2 VUS are blocked for want of functional data?"
    has to be a single aggregation, which is the entire reason this is an enum
    and not free text.
    """

    MISSING_FUNCTIONAL = "MISSING_FUNCTIONAL"
    MISSING_SEGREGATION = "MISSING_SEGREGATION"
    MISSING_CASE_CONTROL = "MISSING_CASE_CONTROL"
    CONFLICTING_INTRINSIC = "CONFLICTING_INTRINSIC"
    PREDICTOR_INDETERMINATE = "PREDICTOR_INDETERMINATE"
    FREQUENCY_UNINFORMATIVE = "FREQUENCY_UNINFORMATIVE"
    RESOLVED_NOT_BLOCKED = "RESOLVED_NOT_BLOCKED"


class CriterionOutcome(str, Enum):
    """The result of testing one criterion against one variant.

    Every criterion in the specification gets exactly one of these for every
    variant. The complete trace *is* the map (spec section 1); the verdict is a
    projection of it.
    """

    APPLIED = "applied"
    NOT_MET = "not_met"
    NOT_EVALUABLE = "not_evaluable"
    NOT_APPLICABLE = "not_applicable"


class SkipReason(str, Enum):
    """Why a criterion was not applied. Machine-readable, never prose."""

    RULE_NOT_MET = "rule_not_met"
    DATA_MISSING = "data_missing"
    EXTRINSIC_EVIDENCE_REQUIRED = "extrinsic_evidence_required"
    CONSEQUENCE_OUT_OF_SCOPE = "consequence_out_of_scope"
    SUPERSEDED_BY_MUTEX = "superseded_by_mutex"
    EXCLUDED_BY_SPEC = "excluded_by_spec"


class PointSystem(BaseModel):
    """Tavtigian points and the classification thresholds, as data.

    Defaults reproduce the table in spec section 2::

        supporting +/-1   moderate +2   strong +/-4   very strong +/-8
        P >= 10 | LP 6..9 | VUS 0..5 | LB -1..-6 | B <= -7
    """

    model_config = ConfigDict(frozen=True)

    pathogenic_points: dict[Strength, int] = Field(
        default_factory=lambda: {
            Strength.SUPPORTING: 1,
            Strength.MODERATE: 2,
            Strength.STRONG: 4,
            Strength.VERY_STRONG: 8,
        }
    )
    benign_points: dict[Strength, int] = Field(
        default_factory=lambda: {
            Strength.SUPPORTING: -1,
            Strength.MODERATE: -2,
            Strength.STRONG: -4,
            Strength.VERY_STRONG: -8,
            Strength.STAND_ALONE: -8,
        }
    )
    threshold_pathogenic: int = 10
    threshold_likely_pathogenic: int = 6
    threshold_likely_benign: int = -1
    threshold_benign: int = -7

    def points_for(self, direction: Direction, strength: Strength) -> int:
        """Signed points contributed by one criterion at one strength."""
        table = self.pathogenic_points if direction is Direction.PATHOGENIC else self.benign_points
        try:
            return table[strength]
        except KeyError as exc:  # pragma: no cover - guarded by spec validation
            raise ValueError(
                f"strength {strength.value!r} is not defined for direction "
                f"{direction.value!r} in this point system"
            ) from exc

    def classify(self, points: int) -> ACMGClass:
        """Map a signed point total onto the five-tier terminology."""
        if points >= self.threshold_pathogenic:
            return ACMGClass.PATHOGENIC
        if points >= self.threshold_likely_pathogenic:
            return ACMGClass.LIKELY_PATHOGENIC
        if points <= self.threshold_benign:
            return ACMGClass.BENIGN
        if points <= self.threshold_likely_benign:
            return ACMGClass.LIKELY_BENIGN
        return ACMGClass.UNCERTAIN

    def gap_to_likely_pathogenic(self, points: int) -> int | None:
        """Points still needed to reach LP, or ``None`` if already there."""
        if points >= self.threshold_likely_pathogenic:
            return None
        return self.threshold_likely_pathogenic - points

    def gap_to_likely_benign(self, points: int) -> int | None:
        """Points still needed (as a positive magnitude) to reach LB."""
        if points <= self.threshold_likely_benign:
            return None
        return points - self.threshold_likely_benign

    def max_strength_points(self, direction: Direction) -> int:
        """Largest magnitude a single criterion can contribute in a direction."""
        table = self.pathogenic_points if direction is Direction.PATHOGENIC else self.benign_points
        return max(abs(v) for v in table.values())


DEFAULT_POINT_SYSTEM = PointSystem()
