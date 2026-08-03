"""vus-foresight -- an evidence gap map for clinical genomics.

For every *possible* variant in a gene, compute the ACMG classification
reachable using only patient-independent evidence, and what specific evidence
would be needed to resolve the remainder.

The product is not a classification. It is a map of where the uncertainty lives
and what it would cost to remove it.
"""

from .acmg import (
    ACMGClass,
    BlockingReason,
    CriterionOutcome,
    Direction,
    EvidenceClass,
    FeasibilityTag,
    PointSystem,
    SkipReason,
    Strength,
)
from .gapmap import (
    AppliedCriterion,
    EvidenceRequirement,
    EvidenceSet,
    GapMapRow,
    SkippedCriterion,
)
from .variant import Consequence, Variant, VariantKind

__version__ = "0.1.0"

__all__ = [
    "ACMGClass",
    "AppliedCriterion",
    "BlockingReason",
    "Consequence",
    "CriterionOutcome",
    "Direction",
    "EvidenceClass",
    "EvidenceRequirement",
    "EvidenceSet",
    "FeasibilityTag",
    "GapMapRow",
    "PointSystem",
    "SkipReason",
    "SkippedCriterion",
    "Strength",
    "Variant",
    "VariantKind",
    "__version__",
]
