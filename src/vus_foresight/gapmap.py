"""The output schema (spec section 6): one row per variant, per gene.

The product is not the classification. It is ``gap_to_LP`` / ``gap_to_LB``,
``minimum_sufficient_sets`` and ``blocking_reason`` -- the answer to "what would
it cost to stop being uncertain about this".
"""

from __future__ import annotations

from datetime import date, datetime

from pydantic import BaseModel, ConfigDict, Field

from .acmg import (
    ACMGClass,
    BlockingReason,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    SkipReason,
    Strength,
)
from .variant import Consequence, VariantKind

__all__ = [
    "GAP_MAP_COLUMNS",
    "AppliedCriterion",
    "EvidenceRequirement",
    "EvidenceSet",
    "GapMapRow",
    "SkippedCriterion",
]


class AppliedCriterion(BaseModel):
    """A criterion that fired, with everything needed to audit why."""

    model_config = ConfigDict(frozen=True)

    code: str
    direction: Direction
    strength: Strength
    points: int
    evidence_class: EvidenceClass
    #: Human-readable summary of the observation that triggered the rule, built
    #: from the spec's own template -- never free text written by the engine.
    evidence: str
    #: Which adapter supplied the data, with its declared version.
    source: str


class SkippedCriterion(BaseModel):
    """A criterion that was tested and did not fire.

    Recording these is the whole point: the trace, not the verdict, is the map.
    """

    model_config = ConfigDict(frozen=True)

    code: str
    direction: Direction
    evidence_class: EvidenceClass
    outcome: CriterionOutcome
    reason: SkipReason
    #: For ``DATA_MISSING``, the dotted context paths that were absent. This is
    #: what makes the ``available_uningested`` report computable.
    missing_fields: tuple[str, ...] = ()
    detail: str | None = None


class EvidenceRequirement(BaseModel):
    """One criterion that would have to be established, and at what strength."""

    model_config = ConfigDict(frozen=True)

    code: str
    direction: Direction
    strength: Strength
    points: int
    evidence_class: EvidenceClass
    feasibility: FeasibilityTag
    #: Number of informative meioses / phenotyped carriers, when the spec states
    #: one for the strength being claimed. ``None`` when not applicable.
    required_observations: int | None = None
    description: str | None = None


class EvidenceSet(BaseModel):
    """A minimal set of evidence sufficient to close the gap (spec section 7)."""

    model_config = ConfigDict(frozen=True)

    target: ACMGClass
    requirements: tuple[EvidenceRequirement, ...]
    total_points: int
    #: Worst (most expensive) feasibility across the requirements -- a set is
    #: only as obtainable as its hardest member.
    feasibility: FeasibilityTag
    #: Sort key: lower is cheaper. Deterministic function of the requirements.
    acquisition_cost: int

    @property
    def codes(self) -> tuple[str, ...]:
        return tuple(r.code for r in self.requirements)


class GapMapRow(BaseModel):
    """One row of the evidence gap map."""

    model_config = ConfigDict(frozen=True)

    # -- Identity ---------------------------------------------------------
    gene: str
    transcript: str
    hgvs_c: str
    hgvs_p: str | None
    grch38_pos: str | None
    consequence: Consequence
    variant_kind: VariantKind
    equivalence_class_id: str
    mutational_distance: int

    # -- Current state ----------------------------------------------------
    criteria_applied: tuple[AppliedCriterion, ...]
    criteria_evaluated_not_applied: tuple[SkippedCriterion, ...]
    points_current: int
    class_current: ACMGClass

    # -- Ceiling ----------------------------------------------------------
    points_ceiling_intrinsic: int
    class_ceiling_intrinsic: ACMGClass

    # -- The gap: the product --------------------------------------------
    gap_to_LP: int | None
    gap_to_LB: int | None
    minimum_sufficient_sets: tuple[EvidenceSet, ...]
    blocking_reason: BlockingReason

    # -- Provenance -------------------------------------------------------
    spec_version: str
    clinvar_snapshot: date | None
    gnomad_version: str | None
    source_versions: dict[str, str] = Field(default_factory=dict)
    computed_at: datetime

    @property
    def variant_id(self) -> str:
        return f"{self.transcript}:{self.hgvs_c}"


#: Explicit, ordered column list for the Parquet writer.
#:
#: Fixing this here rather than deriving it from ``dict`` iteration is part of
#: the determinism guarantee of spec section 11: the schema of the output must
#: not depend on model field ordering surviving a refactor unnoticed.
GAP_MAP_COLUMNS: tuple[str, ...] = (
    "gene",
    "transcript",
    "hgvs_c",
    "hgvs_p",
    "grch38_pos",
    "consequence",
    "variant_kind",
    "equivalence_class_id",
    "mutational_distance",
    "criteria_applied",
    "criteria_evaluated_not_applied",
    "points_current",
    "class_current",
    "points_ceiling_intrinsic",
    "class_ceiling_intrinsic",
    "gap_to_LP",
    "gap_to_LB",
    "minimum_sufficient_sets",
    "blocking_reason",
    "spec_version",
    "clinvar_snapshot",
    "gnomad_version",
    "source_versions",
    "computed_at",
)
