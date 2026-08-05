"""The deterministic classification and gap engine.

Gene-agnostic by construction: it receives ``(variant, gene_config, spec)`` and
holds no domain constant of its own. The test that keeps it that way lives in
``tests/layer3_classification/test_domain_isolation.py`` and runs the whole
suite against a synthetic gene and a toy specification.
"""

from .cnv_scoring import CNVScore, CNVScoringConfig, cnv_row, score_cnv
from .context import MISSING, EvidenceContext, LookupLog
from .equivalence import equivalence_class_id, equivalence_key
from .evaluator import Evaluation, evaluate_variant, render_evidence
from .gap import Candidate, GapAnalysis, analyse_gap, minimum_sufficient_sets
from .pipeline import MapRunner, VariantResult, default_registry
from .predicates import Leaf, RuleNode, evaluate_rule
from .pvs1 import compute_pvs1
from .spec import CriterionSpec, PVS1Config, VCEPSpec, load_spec

__all__ = [
    "MISSING",
    "CNVScore",
    "CNVScoringConfig",
    "Candidate",
    "CriterionSpec",
    "Evaluation",
    "EvidenceContext",
    "GapAnalysis",
    "Leaf",
    "LookupLog",
    "MapRunner",
    "PVS1Config",
    "RuleNode",
    "VCEPSpec",
    "VariantResult",
    "analyse_gap",
    "cnv_row",
    "compute_pvs1",
    "default_registry",
    "equivalence_class_id",
    "equivalence_key",
    "evaluate_rule",
    "evaluate_variant",
    "load_spec",
    "minimum_sufficient_sets",
    "render_evidence",
    "score_cnv",
]
