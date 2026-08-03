"""Layer 3: PS1 and PM5 against a dated ClinVar snapshot.

These are the only criteria whose answer depends on the state of a public
database on a date, and they are the two most likely to be implemented as a
tautology. The tests below are mostly about what must *not* count as evidence:

* a variant's own ClinVar record cannot satisfy PS1 -- PS1 asks whether a
  *different* nucleotide change reaching the same protein change is already
  established, and self-satisfaction would manufacture four points for every
  already-classified variant in the gene;
* the same protein change cannot satisfy PM5, which is PS1's question;
* a nonsense neighbour cannot satisfy PM5, which asks for a different *missense*
  at the residue.
"""

from __future__ import annotations

from datetime import date

import pytest

from vus_foresight.adapters.clinvar import (
    ClinVarRecord,
    ClinVarSnapshot,
    ClinVarSnapshotAdapter,
    normalise_classification,
    protein_change_kind,
    review_stars,
)
from vus_foresight.adapters.clinvar_import import (
    ParseStats,
    build_snapshot,
    parse_variant_summary,
    write_snapshot,
)
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.variant import Consequence


@pytest.fixture()
def missense(minus_gene):
    return next(
        v for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE and v.codon_index and v.codon_index > 3
    )


def _snapshot(records, when=date(2020, 1, 1)):
    return ClinVarSnapshot.from_records(records, when)


def test_a_variant_cannot_satisfy_ps1_from_its_own_record(missense):
    """The whole point. Self-satisfaction is circular, not evidence."""
    snapshot = _snapshot(
        [
            ClinVarRecord(
                hgvs_c=missense.hgvs_c,
                hgvs_p=missense.hgvs_p,
                classification="pathogenic",
                stars=3,
                codon=missense.codon_index,
            )
        ]
    )
    payload = snapshot.resolve(missense, min_stars=1)
    assert "same_protein_change" not in payload
    # It is still reported, for provenance and for the hindsight join.
    assert payload["self"]["classification"] == "pathogenic"


def test_ps1_fires_on_a_different_nucleotide_change_reaching_the_same_protein(missense):
    snapshot = _snapshot(
        [
            ClinVarRecord(
                hgvs_c="c.999999A>G",
                hgvs_p=missense.hgvs_p,
                classification="pathogenic",
                stars=3,
                codon=missense.codon_index,
            )
        ]
    )
    payload = snapshot.resolve(missense, min_stars=1)
    assert payload["same_protein_change"]["classification"] == "pathogenic"
    assert payload["same_protein_change"]["example"] == "c.999999A>G"
    assert payload["same_protein_change"]["count"] == 1


def test_pm5_ignores_the_same_protein_change(missense):
    """That is PS1's question; counting it twice would double the evidence."""
    snapshot = _snapshot(
        [
            ClinVarRecord(
                hgvs_c="c.999999A>G",
                hgvs_p=missense.hgvs_p,
                classification="pathogenic",
                stars=3,
                codon=missense.codon_index,
            )
        ]
    )
    payload = snapshot.resolve(missense, min_stars=1)
    assert "same_protein_change" in payload
    assert "codon" not in payload


def test_pm5_requires_the_neighbour_to_be_a_missense(missense):
    """A nonsense neighbour is a PVS1 observation, not a PM5 one."""
    codon = missense.codon_index
    nonsense_only = _snapshot(
        [
            ClinVarRecord(
                hgvs_c="c.999998C>T",
                hgvs_p=f"p.Gln{codon}Ter",
                classification="pathogenic",
                stars=3,
                codon=codon,
            )
        ]
    )
    assert "codon" not in nonsense_only.resolve(missense, min_stars=1)

    with_missense = _snapshot(
        [
            ClinVarRecord(
                hgvs_c="c.999997C>T",
                hgvs_p=f"p.Gln{codon}Trp",
                classification="likely_pathogenic",
                stars=2,
                codon=codon,
            )
        ]
    )
    payload = with_missense.resolve(missense, min_stars=1)
    assert payload["codon"]["other_change_classification"] == "likely_pathogenic"


def test_review_status_floor_excludes_uncriteried_submissions(missense):
    snapshot = _snapshot(
        [
            ClinVarRecord(
                hgvs_c="c.999999A>G",
                hgvs_p=missense.hgvs_p,
                classification="pathogenic",
                stars=0,
                codon=missense.codon_index,
            )
        ]
    )
    assert "same_protein_change" not in snapshot.resolve(missense, min_stars=1)
    assert "same_protein_change" in snapshot.resolve(missense, min_stars=0)


def test_the_strongest_pathogenic_precedent_wins(missense):
    """PS1 asks whether a pathogenic precedent exists, not what the consensus is."""
    snapshot = _snapshot(
        [
            ClinVarRecord("c.111A>G", missense.hgvs_p, "benign", 3, missense.codon_index),
            ClinVarRecord("c.222A>G", missense.hgvs_p, "pathogenic", 1, missense.codon_index),
            ClinVarRecord("c.333A>G", missense.hgvs_p, "uncertain", 3, missense.codon_index),
        ]
    )
    payload = snapshot.resolve(missense, min_stars=1)
    assert payload["same_protein_change"]["classification"] == "pathogenic"
    assert payload["same_protein_change"]["count"] == 3


def test_the_adapter_reaches_the_criteria_end_to_end(
    tmp_path, minus_gene, minus_config, toy_spec, missense
):
    from datetime import datetime

    from vus_foresight.adapters import AdapterRegistry
    from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
    from vus_foresight.engine.pipeline import MapRunner

    path = tmp_path / "clinvar.tsv"
    write_snapshot(
        [
            ClinVarRecord(
                "c.999999A>G", missense.hgvs_p, "pathogenic", 3, missense.codon_index
            )
        ],
        path,
    )
    registry = AdapterRegistry(
        [VariantAdapter(), TranscriptAdapter(), RegionAdapter(minus_config)]
    )
    registry.add(ClinVarSnapshotAdapter.from_path(path, snapshot_date=date(2020, 1, 1)))
    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=toy_spec,
        adapters=registry,
        computed_at=datetime(1970, 1, 1),
    )
    row = runner.evaluate(missense).row
    applied = {c.code: c for c in row.criteria_applied}
    assert "PS1" in applied
    assert applied["PS1"].points == 4
    assert "c.999999A>G" not in applied["PS1"].evidence  # template names the class, not the key
    assert applied["PS1"].source.startswith("clinvar@")


# --------------------------------------------------------------------------
# The snapshot builder.
# --------------------------------------------------------------------------

VARIANT_SUMMARY = """#AlleleID\tType\tName\tGeneSymbol\tClinicalSignificance\tLastEvaluated\tReviewStatus\tAssembly\tChromosome
1\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5074G>A (p.Asp1692Asn)\tBRCA1\tPathogenic\tMar 21, 2024\treviewed by expert panel\tGRCh38\t17
2\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5075A>G (p.Asp1692Gly)\tBRCA1\tLikely pathogenic\tJan 2, 2020\tcriteria provided, single submitter\tGRCh38\t17
3\tsingle nucleotide variant\tNM_007294.3(BRCA1):c.5074G>A (p.Asp1692Asn)\tBRCA1\tPathogenic\tMar 21, 2024\treviewed by expert panel\tGRCh38\t17
4\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5076C>T (p.Asp1692=)\tBRCA1\tBenign\t-\tcriteria provided, multiple submitters, no conflicts\tGRCh38\t17
5\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5077A>T (p.Ile1693Phe)\tBRCA1\tnot provided\t-\tno assertion criteria provided\tGRCh38\t17
6\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5074G>T (p.Asp1692Tyr)\tBRCA1\tPathogenic\t-\treviewed by expert panel\tGRCh37\t17
7\tDeletion\tGRCh38/hg38 17q21.31\tBRCA1\tPathogenic\t-\treviewed by expert panel\tGRCh38\t17
"""


def test_variant_summary_parser_keeps_only_the_named_transcript(tmp_path):
    source = tmp_path / "variant_summary.txt"
    source.write_text(VARIANT_SUMMARY, encoding="utf-8")
    stats = ParseStats()
    records = list(
        parse_variant_summary(source, transcript="NM_007294.4", stats=stats)
    )

    kept = {r.hgvs_c for r in records}
    assert kept == {"c.5074G>A", "c.5075A>G", "c.5076C>T"}
    # A classification made against NM_007294.3 is not evidence about a
    # coordinate in NM_007294.4 unless somebody checked that the two agree.
    assert stats.other_transcript >= 1
    assert stats.wrong_assembly == 1
    assert stats.unmapped_classification == 1  # "not provided"
    assert stats.unparsable_name == 1  # the cytogenetic deletion


def test_variant_summary_parser_derives_stars_codon_and_date(tmp_path):
    source = tmp_path / "variant_summary.txt"
    source.write_text(VARIANT_SUMMARY, encoding="utf-8")
    records = {r.hgvs_c: r for r in parse_variant_summary(source, transcript="NM_007294.4")}

    expert = records["c.5074G>A"]
    assert expert.stars == 3
    assert expert.codon == 1692
    assert expert.classification == "pathogenic"
    assert expert.last_evaluated == date(2024, 3, 21)
    assert expert.kind == "missense"

    assert records["c.5075A>G"].stars == 1
    assert records["c.5076C>T"].kind == "synonymous"
    assert records["c.5076C>T"].last_evaluated is None


def test_columns_are_read_by_name_so_an_added_column_does_not_shift_them(tmp_path):
    """ClinVar adds columns between releases; a positional parser reads garbage."""
    header, *rest = VARIANT_SUMMARY.strip().splitlines()
    columns = header.split("\t")
    reordered_header = "\t".join([columns[2], *columns[:2], *columns[3:]])
    reordered_rows = [
        "\t".join([(cells := line.split("\t"))[2], *cells[:2], *cells[3:]]) for line in rest
    ]
    source = tmp_path / "reordered.txt"
    source.write_text("\n".join([reordered_header, *reordered_rows]) + "\n", encoding="utf-8")

    records = {r.hgvs_c for r in parse_variant_summary(source, transcript="NM_007294.4")}
    assert records == {"c.5074G>A", "c.5075A>G", "c.5076C>T"}


def test_a_missing_required_column_is_an_error_not_an_empty_snapshot(tmp_path):
    source = tmp_path / "broken.txt"
    source.write_text("Name\tAssembly\nNM_1:c.1A>G\tGRCh38\n", encoding="utf-8")
    with pytest.raises(ValueError, match="ClinicalSignificance"):
        list(parse_variant_summary(source, transcript="NM_1"))


def test_snapshot_build_round_trips_and_is_byte_stable(tmp_path):
    source = tmp_path / "variant_summary.txt"
    source.write_text(VARIANT_SUMMARY, encoding="utf-8")
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    stats = build_snapshot(source, first, transcript="NM_007294.4")
    build_snapshot(source, second, transcript="NM_007294.4")

    assert stats.kept == 3
    assert first.read_bytes() == second.read_bytes()

    snapshot = ClinVarSnapshot.read(first, date(2024, 1, 1))
    assert len(snapshot) == 3
    assert snapshot.by_codon[1692]
    assert snapshot.snapshot_date == date(2024, 1, 1)


def test_vocabulary_helpers():
    assert normalise_classification("Pathogenic/Likely pathogenic") == "pathogenic"
    assert normalise_classification("Uncertain significance") == "uncertain"
    assert normalise_classification("something new") is None
    assert review_stars("reviewed by expert panel") == 3
    assert review_stars("a status nobody has seen") == 0
    assert protein_change_kind("p.Leu100Pro") == "missense"
    assert protein_change_kind("p.Leu100=") == "synonymous"
    assert protein_change_kind("p.Arg100Ter") == "nonsense"
    assert protein_change_kind("p.Arg100SerfsTer12") == "truncating"
    assert protein_change_kind("p.Met1?") == "unpredictable"
    assert protein_change_kind(None) == "unknown"


def test_the_map_never_reads_clinvars_own_call_as_a_criterion(toy_spec):
    """Classifying because ClinVar already did would make the map a mirror."""
    assert "clinvar.self.classification" not in toy_spec.field_footprint
    for criterion in toy_spec.criteria:
        assert not any(field.startswith("clinvar.self") for field in criterion.iter_fields())
