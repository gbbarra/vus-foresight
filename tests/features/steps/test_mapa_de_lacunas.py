"""Steps for mapa_de_lacunas.feature.

These scenarios are about the product itself, so they assert on a real map row
produced by the real pipeline. Nothing is constructed by hand: if the gap
engine stops naming what is missing, or the trace stops carrying its evidence,
these fail.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest
import yaml
from pytest_bdd import given, scenarios, then, when
from typer.testing import CliRunner

from vus_foresight.acmg import ACMGClass, BlockingReason, FeasibilityTag
from vus_foresight.adapters import AdapterRegistry
from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from vus_foresight.cli import app
from vus_foresight.engine.pipeline import MapRunner
from vus_foresight.enumeration import enumerate_coding_snvs, enumerate_exon_cnvs
from vus_foresight.variant import Consequence

scenarios("../mapa_de_lacunas.feature")

CLASSE_POR_NOME = {
    "Patogênica": ACMGClass.PATHOGENIC,
    "Provavelmente Patogênica": ACMGClass.LIKELY_PATHOGENIC,
    "Incerta": ACMGClass.UNCERTAIN,
    "Provavelmente Benigna": ACMGClass.LIKELY_BENIGN,
    "Benigna": ACMGClass.BENIGN,
}

SPEC_DIR = Path(__file__).resolve().parents[3] / "config" / "specs"


@given("um gene cujo mecanismo de doença é perda de função estabelecida")
def gene(mundo, minus_gene, minus_config, toy_spec):
    mundo["gene"] = minus_gene
    mundo["config"] = minus_config.model_copy(update={"lof_mechanism": "established"})
    mundo["spec"] = toy_spec


@given("uma variante de sentido trocado sobre a qual nada foi publicado")
def missense(mundo):
    mundo["variante"] = next(
        v
        for v in enumerate_coding_snvs(mundo["gene"].transcript)
        if v.consequence is Consequence.MISSENSE
    )


@given("uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado")
def nonsense_precoce(mundo):
    transcript = mundo["gene"].transcript
    mundo["variante"] = next(
        v
        for v in enumerate_coding_snvs(transcript)
        if v.consequence is Consequence.NONSENSE
        and v.codon_index
        and not transcript.ptc_escapes_nmd(v.codon_index)
    )


@given("uma deleção que remove o gene inteiro")
def delecao_total(mundo):
    from vus_foresight.variant import VariantKind

    mundo["cnv"] = max(
        (v for v in enumerate_exon_cnvs(mundo["gene"].transcript) if v.kind is VariantKind.CNV),
        key=lambda v: len(v.attributes.get("exons", "")),
    )


@given("uma especificação cujos limiares numéricos ainda não foram conferidos")
def spec_nao_conferida(mundo, tmp_path, minus_gene):
    transcript = minus_gene.transcript
    (tmp_path / "TOY1.fa").write_text(
        f">{transcript.transcript_id}\n{transcript.sequence}\n", encoding="ascii"
    )
    config = {
        "gene": "TOY1",
        "assembly": "TOY",
        "spec": "enigma_brca_v1.1.0",
        "lof_mechanism": "established",
        "transcript": {
            "id": transcript.transcript_id,
            "chrom": transcript.chrom,
            "strand": transcript.strand,
            "cds_length": transcript.cds_length,
            "protein_length": transcript.protein_length,
            "exon_count": len(transcript.exons),
            "exon_labels": [e.label for e in transcript.exons],
            "cds_start_tx": transcript.cds_start_tx,
            "cds_end_tx": transcript.cds_end_tx,
            "exons": [{"label": e.label, "start": e.start, "end": e.end} for e in transcript.exons],
        },
        "sequence": {"path": "TOY1.fa"},
    }
    caminho = tmp_path / "TOY1.yaml"
    caminho.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    mundo["config_path"] = caminho
    mundo["data_root"] = tmp_path


@when("a classificação for calculada")
def calcular(mundo):
    if "cnv" in mundo:
        from vus_foresight.engine.cnv_scoring import CNVScoringConfig, cnv_row

        cnv_config = CNVScoringConfig()
        mundo["cnv_config"] = cnv_config
        mundo["linha"] = cnv_row(
            mundo["cnv"],
            mundo["gene"].transcript,
            mundo["config"],
            cnv_config,
            computed_at=datetime(1970, 1, 1),
        )
        return

    runner = MapRunner(
        transcript=mundo["gene"].transcript,
        gene=mundo["config"],
        spec=mundo["spec"],
        adapters=AdapterRegistry(
            [VariantAdapter(), TranscriptAdapter(), RegionAdapter(mundo["config"])]
        ),
        computed_at=datetime(1970, 1, 1),
    )
    mundo["linha"] = next(iter(runner.run([mundo["variante"]]))).row


@when("alguém pedir o mapa desse gene")
def pedir_o_mapa(mundo, tmp_path):
    mundo["resultado_cli"] = CliRunner().invoke(
        app,
        [
            "map",
            "--gene",
            str(mundo["config_path"]),
            "--data-root",
            str(mundo["data_root"]),
            "--spec-dir",
            str(SPEC_DIR),
            "--out-dir",
            str(tmp_path / "saida"),
        ],
    )
    mundo["saida_esperada"] = tmp_path / "saida"


@then('o resultado deve ser "Incerta"')
def resultado_incerto(mundo):
    assert mundo["linha"].class_current is ACMGClass.UNCERTAIN


@then('o resultado não deve ser "Incerta"')
def resultado_nao_incerto(mundo):
    assert mundo["linha"].class_current is not ACMGClass.UNCERTAIN


@then("o mapa deve nomear qual evidência falta para resolvê-la")
def nomeia_o_que_falta(mundo):
    motivo = mundo["linha"].blocking_reason
    assert motivo is not BlockingReason.RESOLVED_NOT_BLOCKED
    assert motivo.value in {r.value for r in BlockingReason}


@then("deve informar quantos pontos ainda faltam para cada desfecho")
def informa_a_distancia(mundo):
    linha = mundo["linha"]
    assert linha.gap_to_LP is not None and linha.gap_to_LP > 0
    assert linha.gap_to_LB is not None and linha.gap_to_LB > 0


@then("cada conjunto de evidência suficiente deve declarar o quanto é viável obtê-lo")
def declara_viabilidade(mundo):
    conjuntos = mundo["linha"].minimum_sufficient_sets
    assert conjuntos, "nenhum conjunto suficiente foi proposto"
    assert all(c.feasibility in set(FeasibilityTag) for c in conjuntos)


@then(
    "deve haver pelo menos um conjunto cuja evidência já existe publicamente e só não foi ingerida"
)
def ha_conjunto_ja_publico(mundo):
    conjuntos = mundo["linha"].minimum_sufficient_sets
    assert any(c.feasibility is FeasibilityTag.AVAILABLE_UNINGESTED for c in conjuntos), (
        f"viabilidades propostas: {sorted({c.feasibility.value for c in conjuntos})}"
    )


@then("o mapa deve registrar que nada a está travando")
def nada_travando(mundo):
    assert mundo["linha"].blocking_reason is BlockingReason.RESOLVED_NOT_BLOCKED


@then("cada critério aplicado deve declarar a evidência e a fonte que o sustentam")
def criterios_declaram_evidencia(mundo):
    aplicados = mundo["linha"].criteria_applied
    assert aplicados, "nenhum critério aplicado nesta variante; o cenário não testa nada"
    for criterio in aplicados:
        assert criterio.evidence.strip() != ""
        assert criterio.source.strip() != ""


@then("cada critério avaliado e não aplicado deve declarar por que não se aplicou")
def criterios_declaram_o_porque(mundo):
    from vus_foresight.acmg import CriterionOutcome, SkipReason

    ignorados = mundo["linha"].criteria_evaluated_not_applied
    assert ignorados, "nenhum critério avaliado e não aplicado; o cenário não testa nada"
    for ignorado in ignorados:
        assert ignorado.outcome in set(CriterionOutcome)
        assert ignorado.reason in set(SkipReason)


@then("o sistema deve se recusar a escrevê-lo")
def recusa_escrever(mundo):
    assert mundo["resultado_cli"].exit_code != 0
    assert not mundo["saida_esperada"].exists()


@then("deve dizer o que é preciso fazer para prosseguir mesmo assim")
def diz_como_prosseguir(mundo):
    assert "--allow-unverified" in mundo["resultado_cli"].output


@then("o resultado deve ser reportado em uma escala própria")
def escala_propria(mundo):
    """Copy-number evidence is scored on its own scale, so its total must not
    be interpretable as Tavtigian points."""
    linha = mundo["linha"]
    assert linha.spec_version == mundo["cnv_config"].version_string
    assert linha.spec_version != mundo["spec"].spec_version


@then("deve declarar qual escala foi usada")
def declara_a_escala(mundo):
    assert mundo["cnv_config"].version_string in mundo["linha"].spec_version
    if not mundo["linha"].spec_version:  # pragma: no cover - guards the assertion above
        pytest.fail("a linha de número de cópias não declarou a sua escala")
