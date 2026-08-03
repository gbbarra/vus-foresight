"""Layer 3: concordance against expert curation, and the validation study.

Both are data-gated. The suite must stay green on a clone with no curated
fixtures, but the harness has to exist so that dropping the data in is the only
remaining step.

``tests/fixtures/enigma_three_star.tsv``
    Three-star ClinVar records with the criterion codes the expert panel
    applied. Columns: ``gene``, ``hgvs_c``, ``classification``,
    ``criteria`` (semicolon-separated codes), ``intrinsic_only`` (true/false).

Concordance is measured **code by code**, not only on the final verdict. Getting
the class right for the wrong reasons is a silent failure that only surfaces
when the gene changes -- which is exactly when it is most expensive.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vus_foresight.acmg import ACMGClass
from vus_foresight.validation import Outcome, spearman, validate

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ENIGMA_THREE_STAR = FIXTURES / "enigma_three_star.tsv"
HINDSIGHT_OUTCOMES = FIXTURES / "hindsight_outcomes.tsv"


def _read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.requires_reference
def test_criterion_level_concordance_with_expert_panel(brca1, brca2):
    if not ENIGMA_THREE_STAR.exists():
        pytest.skip(f"{ENIGMA_THREE_STAR.name} not present; see the module docstring")
    pytest.skip(
        "criterion-level concordance requires a curated three-star extract and a "
        "verified specification; wire both in and remove this skip"
    )


@pytest.mark.requires_reference
def test_system_never_overestimates_intrinsic_only_variants(brca1, brca2):
    """Restricted to variants the VCEP resolved with intrinsic evidence alone.

    Intrinsic evidence is a subset of the expert's evidence, so the system's
    point total must not exceed the panel's. Exceeding it means a criterion is
    being applied that should not be.
    """
    if not ENIGMA_THREE_STAR.exists():
        pytest.skip(f"{ENIGMA_THREE_STAR.name} not present; see the module docstring")
    pytest.skip(
        "non-overestimation requires a curated three-star extract with per-record "
        "point totals; wire it in and remove this skip"
    )


@pytest.mark.validation_study
def test_section_10_protocol_against_hindsight():
    """The falsifiable claim. A study, run on demand, never part of CI."""
    if not HINDSIGHT_OUTCOMES.exists():
        pytest.skip(f"{HINDSIGHT_OUTCOMES.name} not present")
    pytest.skip("run via 'vus-foresight validate' with a map computed at date T")


# --------------------------------------------------------------------------
# The validation machinery itself is unit-testable without any curated data.
# --------------------------------------------------------------------------


def test_outcome_classifies_resolution_and_direction():
    resolved = Outcome(
        variant_id="NM_1:c.1A>G",
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.LIKELY_PATHOGENIC,
    )
    assert resolved.was_resolved
    assert resolved.direction == "pathogenic"

    unresolved = Outcome(
        variant_id="NM_1:c.2A>G",
        class_at_t=ACMGClass.UNCERTAIN,
        class_at_t_plus_n=ACMGClass.UNCERTAIN,
    )
    assert not unresolved.was_resolved
    assert unresolved.direction is None


def test_spearman_matches_hand_computed_values():
    assert spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert spearman([1, 1, 1], [1, 2, 3]) is None
    assert spearman([1, 2], [1, 2]) is None


def test_validate_scores_a_synthetic_cohort(minus_gene, minus_config, toy_spec):
    from datetime import date, datetime

    from vus_foresight.engine.pipeline import MapRunner, default_registry
    from vus_foresight.enumeration import enumerate_coding_snvs

    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=default_registry(minus_config),
        computed_at=datetime(1970, 1, 1),
    )
    rows = [
        r.row for r in runner.run(list(enumerate_coding_snvs(minus_gene.transcript))[:40])
    ]
    vus_rows = [row for row in rows if row.class_current is ACMGClass.UNCERTAIN]
    assert vus_rows, "the fixture should leave some variants uncertain"

    outcomes = [
        Outcome(
            variant_id=row.variant_id,
            class_at_t=ACMGClass.UNCERTAIN,
            class_at_t_plus_n=ACMGClass.LIKELY_PATHOGENIC,
            evidence_type="functional",
            resolved_on=date(2020, 1, 1 + index % 28),
        )
        for index, row in enumerate(vus_rows[:10])
    ]
    result = validate(rows, outcomes, reference_date=date(2019, 1, 1))
    assert result.considered == len(outcomes)
    assert result.resolved == len(outcomes)
    assert result.resolvability_recall is not None
    assert result.cause_scored == len(outcomes)
    assert "resolvability recall" in result.summary()
