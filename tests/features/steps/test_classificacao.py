"""Steps for classificacao.feature.

The scenarios speak of "evidence adding up to N points", which is the curator's
language. Here that becomes applied criteria at the strengths that sum to N,
built from the shipped point system rather than by asserting on a number the
test itself chose.
"""

from __future__ import annotations

import pytest
from pytest_bdd import given, parsers, scenarios, then, when

from vus_foresight.acmg import (
    ACMGClass,
    BlockingReason,
    Direction,
    EvidenceClass,
    Strength,
)
from vus_foresight.engine.gap import derive_blocking_reason
from vus_foresight.gapmap import AppliedCriterion

scenarios("../classificacao.feature")

CLASSE_POR_NOME = {
    "Patogênica": ACMGClass.PATHOGENIC,
    "Provavelmente Patogênica": ACMGClass.LIKELY_PATHOGENIC,
    "Incerta": ACMGClass.UNCERTAIN,
    "Provavelmente Benigna": ACMGClass.LIKELY_BENIGN,
    "Benigna": ACMGClass.BENIGN,
}


def _forca_declarada(spec, codigo: str) -> tuple[Direction, Strength]:
    """Direction and strength as the *specification* declares them.

    Never as a table in this file. The first version of these steps carried its
    own, and put PM2 at moderate -- while both shipped specifications declare
    PM2_Supporting at supporting, following ClinGen's downgrade. PVS1 plus a
    moderate PM2 sums to 10 and reads Pathogenic; PVS1 plus a supporting one
    sums to 9 and reads Likely Pathogenic. A hand-written table in a test is
    how a scenario starts asserting something the system does not do.
    """
    criterio = spec.by_code(codigo)
    return criterio.direction, criterio.strength


def _criterio(spec, codigo: str) -> AppliedCriterion:
    direcao, forca = _forca_declarada(spec, codigo)
    return AppliedCriterion(
        code=codigo,
        direction=direcao,
        strength=forca,
        points=spec.point_system.points_for(direcao, forca),
        evidence_class=EvidenceClass.INTRINSIC,
        evidence=f"cenário de aceitação: {codigo}",
        source="acceptance",
    )


def _combinacao_que_soma(spec, alvo: int, direcao: Direction) -> list[str]:
    """Codes from the shipped specification whose points add up to ``alvo``.

    Chosen by search over what the specification actually declares, so that
    changing a criterion's strength changes what these scenarios exercise
    rather than silently leaving them testing a table nobody maintains.
    """
    sistema = spec.point_system
    disponiveis = [
        c.code
        for c in spec.criteria
        if c.direction is direcao
        and c.strength in sistema.pathogenic_points | sistema.benign_points
    ]
    escolhidos: list[str] = []
    restante = alvo
    for codigo in sorted(
        disponiveis, key=lambda c: -abs(sistema.points_for(*_forca_declarada(spec, c)))
    ):
        pontos = sistema.points_for(*_forca_declarada(spec, codigo))
        while abs(restante) >= abs(pontos) and (restante > 0) == (pontos > 0):
            escolhidos.append(codigo)
            restante -= pontos
    if restante != 0:
        pytest.fail(f"nenhuma combinação de critérios da spec soma exatamente {alvo}")
    return escolhidos


@given("um gene cujo mecanismo de doença é perda de função estabelecida")
def gene_com_lof_estabelecida(mundo):
    mundo["lof"] = "established"
    mundo.setdefault("aplicados", [])


@given(parsers.parse("uma variante com evidência somando {pontos:d} pontos"))
def variante_com_pontos(mundo, pontos, toy_spec):
    direcao = Direction.PATHOGENIC if pontos >= 0 else Direction.BENIGN
    mundo["aplicados"] = [
        _criterio(toy_spec, c) for c in _combinacao_que_soma(toy_spec, pontos, direcao)
    ]


@given(parsers.parse("uma variante com evidência de patogenicidade somando {pontos:d} pontos"))
def evidencia_a_favor(mundo, pontos, toy_spec):
    mundo["aplicados"] = [
        _criterio(toy_spec, c) for c in _combinacao_que_soma(toy_spec, pontos, Direction.PATHOGENIC)
    ]


@given(parsers.parse("evidência de benignidade somando {pontos:d} pontos"))
def evidencia_contra(mundo, pontos, toy_spec):
    mundo["aplicados"] += [
        _criterio(toy_spec, c) for c in _combinacao_que_soma(toy_spec, pontos, Direction.BENIGN)
    ]


@given(parsers.parse("uma variante com o critério {codigo} atribuído"))
@given(parsers.parse("o critério {codigo} atribuído"))
def criterio_atribuido(mundo, codigo, toy_spec):
    mundo["aplicados"] = [*mundo.get("aplicados", []), _criterio(toy_spec, codigo)]


@when("a classificação for calculada")
def calcular(mundo, toy_spec):
    """Everything here is the shipped code. Nothing is recomputed locally.

    An earlier version of this step summed the points and re-derived the
    conflict rule inside the test, which would have asserted on a
    reimplementation rather than on the system -- exactly the shape the
    standard forbids.
    """
    from vus_foresight.engine.evaluator import Evaluation

    avaliacao = Evaluation(applied=list(mundo["aplicados"]))
    avaliacao.points = sum(c.points for c in avaliacao.applied)
    avaliacao.acmg_class = toy_spec.point_system.classify(avaliacao.points)

    mundo["avaliacao"] = avaliacao
    mundo["pontos"] = avaliacao.points
    mundo["classe"] = avaliacao.acmg_class
    mundo["conflito"] = avaliacao.has_directional_conflict(
        toy_spec.blocking.conflict_min_points_each
    )
    mundo["motivo"] = derive_blocking_reason(avaliacao, toy_spec, [], ())


@then(parsers.parse('o resultado deve ser "{nome}"'))
def resultado_deve_ser(mundo, nome):
    assert mundo["classe"] is CLASSE_POR_NOME[nome], (
        f"{mundo['pontos']} pontos produziram {mundo['classe'].value}"
    )


@then("o mapa deve registrar que o que a trava é o conflito entre as evidências")
def registra_conflito(mundo):
    assert mundo["classe"] is ACMGClass.UNCERTAIN
    assert mundo["motivo"] is BlockingReason.CONFLICTING_INTRINSIC


@then("o mapa não deve registrar conflito")
def sem_conflito(mundo):
    assert mundo["conflito"] is False
    assert mundo["motivo"] is not BlockingReason.CONFLICTING_INTRINSIC


@then("deve nomear qual evidência falta para resolvê-la")
def nomeia_o_que_falta(mundo):
    assert mundo["motivo"] is not BlockingReason.RESOLVED_NOT_BLOCKED


@given("uma variante cuja frequência populacional satisfaz dois critérios de frequência")
def frequencia_satisfaz_dois(mundo, tmp_path, minus_gene, minus_config, toy_spec):
    """Runs the whole engine, because mutex resolution is the engine's job.

    The shipped frequency rules are written disjoint, so one is widened here to
    create the overlap the scenario describes -- otherwise the rule under test
    could never fire.
    """
    from datetime import datetime

    from vus_foresight.adapters import AdapterRegistry
    from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
    from vus_foresight.adapters.tabular import FrequencyAdapter
    from vus_foresight.engine.pipeline import MapRunner
    from vus_foresight.engine.predicates import Leaf, RuleNode
    from vus_foresight.enumeration import enumerate_coding_snvs
    from vus_foresight.variant import Consequence

    alargada = RuleNode(leaf=Leaf(field="frequency.gnomad.faf95_popmax", op="ge", value=0.0001))
    spec = toy_spec.model_copy(
        update={
            "criteria": tuple(
                c.model_copy(update={"rule": alargada}) if c.code == "BS1" else c
                for c in toy_spec.criteria
            )
        }
    )
    variante = next(
        v
        for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE
    )
    snapshot = tmp_path / "frequencia.tsv"
    snapshot.write_text(
        f"grch38_pos\tgnomad.faf95_popmax\n{variante.grch38_pos}\t0.005\n", encoding="utf-8"
    )
    runner = MapRunner(
        transcript=minus_gene.transcript,
        gene=minus_config,
        spec=spec,
        adapters=AdapterRegistry(
            [
                VariantAdapter(),
                TranscriptAdapter(),
                RegionAdapter(minus_config),
                FrequencyAdapter(snapshot, "acceptance"),
            ]
        ),
        computed_at=datetime(1970, 1, 1),
    )
    mundo["linha"] = next(iter(runner.run([variante]))).row


@then("apenas o critério de frequência mais forte deve contribuir com pontos")
def apenas_o_mais_forte(mundo):
    aplicados = {c.code for c in mundo["linha"].criteria_applied}
    assert "BA1" in aplicados
    assert "BS1" not in aplicados


@then("o critério descartado deve continuar registrado, com o motivo do descarte")
def descartado_fica_registrado(mundo):
    from vus_foresight.acmg import SkipReason

    descartado = next(s for s in mundo["linha"].criteria_evaluated_not_applied if s.code == "BS1")
    assert descartado.reason is SkipReason.SUPERSEDED_BY_MUTEX
