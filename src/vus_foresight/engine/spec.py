"""VCEP specifications as versioned data (spec section 9, principle 2).

The engine interprets; it does not know BRCA. Everything a specification asserts
-- which criteria exist, what they compare against, how strength is modulated,
what a missing piece of evidence would cost to obtain -- lives in YAML under
``config/specs`` and is addressed by ``spec_version`` in every output row.

A change to a threshold is therefore a diff in a versioned file, which is what
makes the golden-snapshot regression test of spec section 11 meaningful.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from functools import cached_property
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..acmg import (
    STRENGTH_ORDER,
    BlockingReason,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    PointSystem,
    Strength,
)
from ..variant import Consequence
from .predicates import RuleNode

#: ``{dotted.path}`` tokens inside an evidence template. Shared with the
#: evaluator so the footprint and the renderer can never disagree about what
#: counts as a token.
TEMPLATE_TOKEN = re.compile(r"\{([A-Za-z0-9_.]+)\}")

__all__ = [
    "TEMPLATE_TOKEN",
    "BlockingConfig",
    "CriterionSpec",
    "EquivalenceConfig",
    "GapConfig",
    "PVS1Config",
    "StrengthRung",
    "VCEPSpec",
    "load_spec",
]


class StrengthRung(BaseModel):
    """One rung of a strength ladder: a condition and the strength it grants."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    strength: Strength
    when: RuleNode


class FeasibilityRung(BaseModel):
    """A condition under which acquiring this evidence is easier than the default.

    The case this exists for: PS3 in a region an SGE screen already covers is
    ``available_uningested`` -- somebody's supplementary table has the answer --
    while PS3 in a region nobody has assayed is ``assay_feasible``. Collapsing
    the two would bury the single most actionable finding the system produces.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    feasibility: FeasibilityTag
    when: RuleNode


class CriterionSpec(BaseModel):
    """One ACMG criterion, exactly as the VCEP defines it."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    code: str
    direction: Direction
    evidence_class: EvidenceClass
    #: Strength granted when the rule fires and no ladder rung matches.
    strength: Strength = Strength.SUPPORTING
    #: Evaluated in order; the first matching rung wins.
    strength_ladder: tuple[StrengthRung, ...] = ()
    #: Highest strength this criterion could plausibly reach if the missing
    #: evidence were obtained. Drives the minimum-sufficient-set search.
    #: Defaults to ``strength``, or the top of the ladder when one is declared.
    max_plausible_strength: Strength | None = None
    #: Criteria in the same group are mutually exclusive; only the strongest
    #: applied member survives. BA1/BS1/PM2 share one, PP3/BP4 another.
    mutex_group: str | None = None
    #: Consequence gate. Empty means "any consequence".
    applies_to: tuple[Consequence, ...] = ()
    #: Context paths that must be present for the criterion to be evaluable at
    #: all. Absent paths yield NOT_EVALUABLE, never NOT_MET.
    requires: tuple[str, ...] = ()
    rule: RuleNode | None = None
    feasibility: FeasibilityTag | None = None
    #: Evaluated in order; the first matching rung overrides ``feasibility``.
    feasibility_ladder: tuple[FeasibilityRung, ...] = ()
    #: Observations needed per strength, e.g. informative meioses for PP1.
    required_observations: dict[Strength, int] = Field(default_factory=dict)
    #: ``str.format``-style template over dotted context paths.
    evidence_template: str | None = None
    #: Which blocking reason this criterion's absence implies.
    blocks_as: BlockingReason | None = None
    description: str | None = None
    notes: str | None = None

    @model_validator(mode="after")
    def _coherent(self) -> CriterionSpec:
        if self.evidence_class is EvidenceClass.EXTRINSIC and self.rule is not None:
            raise ValueError(
                f"{self.code}: extrinsic criteria must not carry a rule -- this system "
                "never sees the observation that would satisfy one"
            )
        if self.evidence_class is not EvidenceClass.EXTRINSIC and self.rule is None:
            raise ValueError(f"{self.code}: intrinsic and semi-intrinsic criteria need a rule")
        if self.evidence_class is EvidenceClass.EXTRINSIC and self.feasibility is None:
            raise ValueError(f"{self.code}: extrinsic criteria must declare a feasibility tag")
        for rung in self.strength_ladder:
            if self.direction is Direction.PATHOGENIC and rung.strength is Strength.STAND_ALONE:
                raise ValueError(f"{self.code}: stand-alone is a benign-only strength")
        return self

    @property
    def ceiling_strength(self) -> Strength:
        """Strongest strength this criterion could ever contribute."""
        if self.max_plausible_strength is not None:
            return self.max_plausible_strength
        candidates = [self.strength, *(rung.strength for rung in self.strength_ladder)]
        return max(candidates, key=STRENGTH_ORDER.index)

    def observations_for(self, strength: Strength) -> int | None:
        return self.required_observations.get(strength)

    def iter_fields(self) -> Iterator[str]:
        """Every context path this criterion can read, from any of its parts.

        Prerequisites, the rule itself, both ladders, and the evidence template
        -- the template counts because it is rendered into the output, so two
        variants that differ only in a templated value are not interchangeable.
        """
        yield from self.requires
        if self.rule is not None:
            yield from self.rule.iter_fields()
        for rung in self.strength_ladder:
            yield from rung.when.iter_fields()
        for rung in self.feasibility_ladder:
            yield from rung.when.iter_fields()
        if self.evidence_template:
            yield from TEMPLATE_TOKEN.findall(self.evidence_template)


class PVS1Config(BaseModel):
    """Parameters of the PVS1 decision tree (Abou Tayoun et al., autoPVS1).

    The tree itself is generic -- NMD zone, fraction of protein removed, overlap
    with a region the *gene config* marks critical. Only these dials are
    specification-dependent, so only these dials are here.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    enabled: bool = True
    #: PTC in an NMD-competent position, in a gene with established LoF mechanism.
    nmd_triggered_strength: Strength = Strength.VERY_STRONG
    #: Same, where the gene's LoF mechanism is not established.
    nmd_triggered_unestablished_strength: Strength = Strength.STRONG
    #: PTC escapes NMD but the truncated region carries a critical tag.
    nmd_escape_critical_strength: Strength = Strength.STRONG
    #: PTC escapes NMD, no critical region, but a large fraction is removed.
    nmd_escape_large_strength: Strength = Strength.STRONG
    #: PTC escapes NMD and little is removed.
    nmd_escape_small_strength: Strength = Strength.MODERATE
    #: Threshold separating "large" from "small", as a fraction of the protein.
    min_fraction_removed: float = 0.10
    #: Region tags that count as critical when checking the truncated segment.
    critical_region_tags: tuple[str, ...] = ("critical",)
    start_lost_strength: Strength = Strength.MODERATE
    #: Whether a canonical splice-site variant needs a splicing prediction
    #: before PVS1 can be assigned. Leaving this true is what makes an
    #: un-ingested SpliceAI table show up as an actionable gap rather than as a
    #: silent assumption.
    splice_requires_prediction: bool = True
    #: Context path holding that prediction, and the threshold above which the
    #: site is taken to be disrupted.
    splice_prediction_field: str = "splice.ds_max"
    splice_prediction_threshold: float = 0.5
    #: Strength when a canonical site is predicted disrupted and the resulting
    #: transcript is presumed NMD-competent.
    splice_disrupted_strength: Strength = Strength.VERY_STRONG


class BlockingConfig(BaseModel):
    """How ``blocking_reason`` is decided, as data rather than as code."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    conflict_reason: BlockingReason = BlockingReason.CONFLICTING_INTRINSIC
    #: Minimum points each direction must contribute before opposing criteria
    #: count as a *conflict*.
    #:
    #: Without this, PM2_Supporting (+1) alongside BP7 (-1) would be reported as
    #: a conflict, and since that pair fires for essentially every rare
    #: synonymous variant it would swamp the aggregation the enum exists for. A
    #: single supporting criterion on each side is routine curation, not a
    #: blocked variant; moderate strength on both sides is a real disagreement.
    conflict_min_points_each: int = 2
    #: First match wins.
    priority: tuple[BlockingReason, ...] = (
        BlockingReason.MISSING_FUNCTIONAL,
        BlockingReason.MISSING_SEGREGATION,
        BlockingReason.MISSING_CASE_CONTROL,
        BlockingReason.PREDICTOR_INDETERMINATE,
        BlockingReason.FREQUENCY_UNINFORMATIVE,
    )
    default: BlockingReason = BlockingReason.PREDICTOR_INDETERMINATE


class EquivalenceConfig(BaseModel):
    """Which consequences may never be collapsed into a class (spec section 5)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    #: Missense is here because PP3/BP4 and PS3/BS3 are per-variant, so two
    #: missense variants sharing an evidence *profile* today can diverge the
    #: moment a predictor score or an assay result lands.
    never_collapse: tuple[Consequence, ...] = (Consequence.MISSENSE,)


class GapConfig(BaseModel):
    """Bounds on the minimum-sufficient-set search (spec section 7)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    #: Largest number of criteria in a proposed set. The search is exhaustive
    #: with pruning over a few dozen criteria, so this is a report-legibility
    #: bound, not a performance one.
    max_cardinality: int = 4
    #: How many sets to keep per target class, cheapest first.
    max_sets_per_target: int = 5
    #: Whether sets containing an intractable requirement are reported at all.
    include_intractable: bool = False


class VCEPSpec(BaseModel):
    """A complete, versioned specification."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    spec_id: str
    name: str
    version: str
    #: Informational only; the engine never gates on it.
    genes: tuple[str, ...] = ()
    source: str | None = None
    #: Whether every numeric threshold in this file has been checked against the
    #: published VCEP document. A specification is structurally valid long before
    #: it is *curated*, and a map computed from unverified thresholds is a
    #: demonstration, not a result. The CLI refuses to write an unverified map
    #: without an explicit override.
    verified: bool = False
    point_system: PointSystem = Field(default_factory=PointSystem)
    pvs1: PVS1Config = Field(default_factory=PVS1Config)
    blocking: BlockingConfig = Field(default_factory=BlockingConfig)
    equivalence: EquivalenceConfig = Field(default_factory=EquivalenceConfig)
    gap: GapConfig = Field(default_factory=GapConfig)
    criteria: tuple[CriterionSpec, ...]

    @model_validator(mode="after")
    def _unique_codes(self) -> VCEPSpec:
        codes = [c.code for c in self.criteria]
        duplicates = sorted({c for c in codes if codes.count(c) > 1})
        if duplicates:
            raise ValueError(f"{self.spec_id}: duplicate criterion codes {duplicates}")
        if not self.criteria:
            raise ValueError(f"{self.spec_id}: a specification with no criteria is useless")
        return self

    @property
    def spec_version(self) -> str:
        return f"{self.spec_id}@{self.version}"

    def by_code(self, code: str) -> CriterionSpec:
        for criterion in self.criteria:
            if criterion.code == code:
                return criterion
        raise KeyError(f"{self.spec_id} has no criterion {code!r}")

    def points_for(self, criterion: CriterionSpec, strength: Strength) -> int:
        return self.point_system.points_for(criterion.direction, strength)

    @cached_property
    def field_footprint(self) -> tuple[str, ...]:
        """Every context path this specification can read, sorted.

        This is the contract behind class-level evaluation: two variants whose
        contexts agree on all of these paths *cannot* evaluate differently,
        because there is nothing else for the engine to look at.

        Getting it wrong would silently merge variants that should differ, so it
        is not left to careful reading of the call sites --
        ``EvidenceContext.audit`` records every path an evaluation actually
        touches, and a test asserts the recorded set is contained here.
        """
        paths: set[str] = set()
        for criterion in self.criteria:
            paths.update(criterion.iter_fields())
        # Read by the PVS1 tree while the context is being assembled, and by the
        # evaluator when it explains why PVS1 could not be determined.
        paths.add(self.pvs1.splice_prediction_field)
        paths.add("pvs1.undetermined_reason")
        return tuple(sorted(paths))


def load_spec(path: str | Path) -> VCEPSpec:
    """Load and validate a specification YAML."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        raw: Any = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: specification must be a mapping")
    return VCEPSpec.model_validate(raw)
