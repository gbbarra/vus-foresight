"""Steps for evidencia_de_vizinho.feature.

Everything here goes through the whole pipeline against dated snapshots of the
public database, because that is the only way the claim is testable: a variant
changing class with nothing new known about it is a statement about two runs,
not about one call.
"""

from __future__ import annotations

from datetime import date, datetime

import pytest
from pytest_bdd import given, parsers, scenarios, then, when

from vus_foresight.adapters import AdapterRegistry, ClinVarSnapshotAdapter
from vus_foresight.adapters.builtin import RegionAdapter, TranscriptAdapter, VariantAdapter
from vus_foresight.adapters.clinvar import ClinVarRecord
from vus_foresight.engine.pipeline import MapRunner
from vus_foresight.engine.timeline import TransitionCause, compare_maps
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.variant import Consequence

scenarios("../evidencia_de_vizinho.feature")


def _mapear(mundo, registros, gene, config, spec, quando: date):
    """One map over the same variants, against the snapshot given."""
    from vus_foresight.adapters.clinvar import ClinVarSnapshot

    registro = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(config)])
    registro.add(
        ClinVarSnapshotAdapter(
            ClinVarSnapshot.from_records(registros, quando), f"clinvar@{quando.isoformat()}"
        )
    )
    runner = MapRunner(
        transcript=gene.transcript,
        gene=config,
        spec=spec,
        adapters=registro,
        computed_at=datetime(1970, 1, 1),
        clinvar_snapshot=quando,
    )
    return [r.row for r in runner.run(mundo["variantes"])]


@given("um gene cujo mecanismo de doença é perda de função estabelecida")
def gene(mundo, minus_gene, minus_config, toy_spec):
    mundo["gene"] = minus_gene
    mundo["config"] = minus_config
    mundo["spec"] = toy_spec
    mundo["registros_antes"] = []


@given("uma variante de sentido trocado sobre a qual nada foi publicado")
def variante_missense(mundo):
    transcript = mundo["gene"].transcript
    candidatas = [
        v
        for v in enumerate_coding_snvs(transcript)
        if v.consequence is Consequence.MISSENSE and v.codon_index and v.codon_index > 3
    ]
    mundo["alvo"] = candidatas[0]
    # The whole codon plus a neighbouring one, so that "only this codon moved"
    # is something the scenario can actually check.
    mundo["variantes"] = candidatas[:40]


@given("que nenhuma outra variante do mesmo códon está classificada")
def nenhum_vizinho_classificado(mundo):
    mundo["registros_antes"] = []


@given("que a própria variante já está classificada como patogênica no banco público")
def a_propria_esta_classificada(mundo):
    alvo = mundo["alvo"]
    mundo["registros_antes"] = [
        ClinVarRecord(alvo.hgvs_c, alvo.hgvs_p, "pathogenic", 3, alvo.codon_index)
    ]
    mundo["registros_depois"] = list(mundo["registros_antes"])


@when("outra variante do mesmo códon for classificada como patogênica")
def vizinha_classificada(mundo):
    alvo = mundo["alvo"]
    mundo["registros_depois"] = [
        *mundo["registros_antes"],
        # A different nucleotide change reaching the same protein change, so the
        # neighbour is genuinely a neighbour and not the variant itself.
        ClinVarRecord("c.999999A>G", alvo.hgvs_p, "pathogenic", 3, alvo.codon_index),
    ]


@when(
    parsers.parse(
        "outra variante do mesmo códon for classificada como patogênica com {estrelas:d} estrela(s)"
    )
)
def vizinha_classificada_com_estrelas(mundo, estrelas):
    alvo = mundo["alvo"]
    mundo["registros_depois"] = [
        *mundo["registros_antes"],
        ClinVarRecord("c.999999A>G", alvo.hgvs_p, "pathogenic", estrelas, alvo.codon_index),
    ]


@when("a classificação for recalculada")
@when("a classificação for calculada")
def recalcular(mundo):
    gene, config, spec = mundo["gene"], mundo["config"], mundo["spec"]
    mundo["antes"] = _mapear(mundo, mundo["registros_antes"], gene, config, spec, date(2018, 1, 1))
    mundo["depois"] = _mapear(
        mundo,
        mundo.get("registros_depois", mundo["registros_antes"]),
        gene,
        config,
        spec,
        date(2024, 1, 1),
    )
    mundo["diff"] = compare_maps(mundo["antes"], mundo["depois"])
    mundo["transicao"] = next(
        (t for t in mundo["diff"].transitions if t.variant_id == mundo["alvo"].variant_id), None
    )


@then("a variante deve ganhar evidência")
def ganhou_evidencia(mundo):
    assert mundo["transicao"] is not None, "a variante alvo não se moveu"
    assert mundo["transicao"].criteria_gained != ()


@then("essa evidência deve ser atribuída à classificação da vizinha")
def atribuida_ao_vizinho(mundo):
    assert mundo["transicao"].cause is TransitionCause.NEIGHBOUR_EVIDENCE


@then("nenhuma evidência nova sobre a própria variante deve ter sido usada")
def nenhuma_evidencia_propria(mundo):
    """The published record about the variant itself is identical in both runs.

    Which is what makes the claim falsifiable rather than rhetorical: the only
    thing that changed between the two maps is a record about a *different*
    nucleotide change.
    """
    alvo = mundo["alvo"]
    antes = [r for r in mundo["registros_antes"] if r.hgvs_c == alvo.hgvs_c]
    depois = [r for r in mundo.get("registros_depois", []) if r.hgvs_c == alvo.hgvs_c]
    assert antes == depois


@then("todas as variantes daquele códon devem ter se movido")
def o_codon_inteiro_se_moveu(mundo):
    codon = mundo["alvo"].codon_index
    no_codon = {v.variant_id for v in mundo["variantes"] if v.codon_index == codon}
    moveram = {t.variant_id for t in mundo["diff"].transitions}
    assert no_codon <= moveram, f"{len(no_codon - moveram)} variantes do códon não se moveram"


@then("nenhuma variante de outro códon deve ter se movido")
def nenhum_outro_codon(mundo):
    codon = mundo["alvo"].codon_index
    por_variante = {v.variant_id: v.codon_index for v in mundo["variantes"]}
    outros = {
        t.variant_id for t in mundo["diff"].transitions if por_variante[t.variant_id] != codon
    }
    assert outros == set()


@then("a variante não deve ganhar evidência de vizinho")
def sem_evidencia_de_vizinho(mundo):
    """Its own pathogenic record must not become evidence about itself.

    Without the exclusion a classified variant would support PS1 for itself,
    fabricating four points out of its own conclusion.
    """
    linha = next(r for r in mundo["depois"] if r.variant_id == mundo["alvo"].variant_id)
    aplicados = {c.code for c in linha.criteria_applied}
    assert "PS1" not in aplicados
    assert "PM5" not in aplicados


@then("o mapa deve continuar dizendo que falta evidência para resolvê-la")
def continua_travada(mundo):
    from vus_foresight.acmg import BlockingReason

    linha = next(r for r in mundo["depois"] if r.variant_id == mundo["alvo"].variant_id)
    assert linha.blocking_reason is not BlockingReason.RESOLVED_NOT_BLOCKED


@then("nenhum critério deve ter consultado a classificação publicada da própria variante")
def ninguem_leu_a_propria_classificacao(mundo):
    """Audited, not assumed. The adapter publishes the field precisely so the
    validation study can use it as an oracle, which only holds if no criterion
    reads it."""
    gene, config, spec = mundo["gene"], mundo["config"], mundo["spec"]
    from vus_foresight.adapters.clinvar import ClinVarSnapshot

    registro = AdapterRegistry([VariantAdapter(), TranscriptAdapter(), RegionAdapter(config)])
    registro.add(
        ClinVarSnapshotAdapter(
            ClinVarSnapshot.from_records(mundo["registros_depois"], date(2024, 1, 1)),
            "clinvar@2024-01-01",
        )
    )
    runner = MapRunner(
        transcript=gene.transcript,
        gene=config,
        spec=spec,
        adapters=registro,
        computed_at=datetime(1970, 1, 1),
        audit_context_reads=True,
    )
    list(runner.run(mundo["variantes"]))

    lidos = runner.audited_paths
    assert lidos, "a auditoria não registrou leitura nenhuma"
    assert "clinvar.self.classification" not in lidos


@then(parsers.parse("a variante {veredito} ter ganhado evidência"))
def veredito_por_estrelas(mundo, veredito):
    if veredito == "deve":
        assert mundo["transicao"] is not None
        assert mundo["transicao"].criteria_gained != ()
    elif veredito == "não deve":
        assert mundo["transicao"] is None
    else:  # pragma: no cover - guards the scenario table itself
        pytest.fail(f"veredito não reconhecido no cenário: {veredito!r}")
