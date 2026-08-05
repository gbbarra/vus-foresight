"""Steps for perda_de_funcao.feature.

The scenarios describe where a protein is interrupted and what that is worth.
Here that becomes the shipped decision tree, called directly, with variants
taken from the real enumeration rather than hand-built -- so a scenario about
"interrupting the protein early enough for the transcript to be degraded" is
answered by the same NMD geometry the map uses.
"""

from __future__ import annotations

import pytest
from pytest_bdd import given, scenarios, then, when

from vus_foresight.acmg import Strength
from vus_foresight.engine.context import EvidenceContext
from vus_foresight.engine.pvs1 import compute_pvs1
from vus_foresight.enumeration import enumerate_coding_snvs
from vus_foresight.variant import Consequence, Variant, VariantKind

scenarios("../perda_de_funcao.feature")


def _nonsense_no_codon(transcript, codon: int) -> Variant | None:
    """A nonsense SNV at ``codon``, if one exists there."""
    return next(
        (
            v
            for v in enumerate_coding_snvs(transcript)
            if v.consequence is Consequence.NONSENSE and v.codon_index == codon
        ),
        None,
    )


@given("um gene cujo mecanismo de doença é perda de função estabelecida")
def lof_estabelecida(mundo, minus_config):
    mundo["gene"] = minus_config.model_copy(update={"lof_mechanism": "established"})


@given("um gene cujo mecanismo de doença não está estabelecido")
def lof_nao_estabelecida(mundo, minus_config):
    mundo["gene"] = minus_config.model_copy(update={"lof_mechanism": "not_established"})


@given("uma especificação que desliga a avaliação de perda de função")
def spec_desligada(mundo, minus_config, toy_spec):
    mundo.setdefault("gene", minus_config)
    mundo["pvs1_config"] = toy_spec.pvs1.model_copy(update={"enabled": False})


@given("uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado")
def interrompe_cedo(mundo, minus_gene):
    transcript = minus_gene.transcript
    variante = next(
        (
            v
            for codon in range(2, transcript.protein_length)
            if (v := _nonsense_no_codon(transcript, codon)) is not None
            and not transcript.ptc_escapes_nmd(codon)
        ),
        None,
    )
    if variante is None:
        pytest.fail("o gene sintético não produz nenhuma interrupção em zona de degradação")
    mundo["variante"] = variante


@given(
    "uma variante que interrompe a proteína tarde, escapando da degradação, "
    "mas removendo uma região crítica"
)
def interrompe_tarde_em_regiao_critica(mundo, plus_gene, plus_config):
    """A C-terminal critical domain, which is the case this branch is for.

    The synthetic gene's own critical region sits at residues 5-12 while NMD
    escape only begins at codon 29, so the two never meet and the branch is
    unreachable with the fixture as shipped. The configuration is therefore
    given a critical region over the last residues -- which is the real shape
    of the problem: BRCA1's BRCT repeats run 1646-1859 of 1863, entirely inside
    the NMD-escaping zone, and a stop there removes them while the transcript
    survives.
    """
    from vus_foresight.genome.reference import FunctionalRegion

    transcript = plus_gene.transcript
    escapam = [
        codon
        for codon in range(2, transcript.protein_length + 1)
        if transcript.ptc_escapes_nmd(codon) and _nonsense_no_codon(transcript, codon) is not None
    ]
    if not escapam:
        pytest.fail("o gene sintético não produz nenhuma interrupção que escape da degradação")
    alvo = escapam[0]
    gene = plus_config.model_copy(
        update={
            "lof_mechanism": "established",
            "functional_regions": (
                *plus_config.functional_regions,
                FunctionalRegion(
                    name="domínio crítico C-terminal",
                    start_aa=alvo,
                    end_aa=transcript.protein_length,
                    tags=("critical", "domain"),
                ),
            ),
        }
    )
    mundo["variante"] = _nonsense_no_codon(transcript, alvo)
    mundo["gene"] = gene
    mundo["transcrito"] = transcript


@given("uma variante truncante cuja posição de interrupção é desconhecida")
def truncante_sem_posicao(mundo, minus_gene):
    mundo["variante"] = Variant(
        transcript=minus_gene.transcript.transcript_id,
        gene=minus_gene.transcript.gene,
        hgvs_c="c.100del",
        hgvs_p=None,
        consequence=Consequence.FRAMESHIFT,
        kind=VariantKind.FRAMESHIFT_CLASS,
        codon_index=None,
        ptc_codon=None,
    )


@given("uma variante de sentido trocado sobre a qual nada foi publicado")
def missense_sem_dados(mundo, minus_gene):
    mundo["variante"] = next(
        v
        for v in enumerate_coding_snvs(minus_gene.transcript)
        if v.consequence is Consequence.MISSENSE
    )


@when("a evidência de perda de função for avaliada")
def avaliar_lof(mundo, minus_gene, toy_spec):
    transcrito = mundo.get("transcrito", minus_gene.transcript)
    config = mundo.get("pvs1_config", toy_spec.pvs1)
    mundo["resultado"] = compute_pvs1(
        mundo["variante"], transcrito, mundo["gene"], config, EvidenceContext()
    )


@then("ela deve ser aplicada com força muito forte")
def forca_muito_forte(mundo):
    assert mundo["resultado"]["applicable"] is True
    assert mundo["resultado"]["strength"] == Strength.VERY_STRONG.value


@then("ela deve ser aplicada com força menor que muito forte")
def forca_menor(mundo, toy_spec):
    resultado = mundo["resultado"]
    assert resultado["applicable"] is True

    pontos = toy_spec.point_system.pathogenic_points
    aplicada = Strength(resultado["strength"])
    assert pontos[aplicada] < pontos[Strength.VERY_STRONG], (
        f"{aplicada.value} não vale menos que muito forte"
    )


@then("ela deve ser aplicada")
def deve_ser_aplicada(mundo):
    assert mundo["resultado"]["applicable"] is True
    assert "strength" in mundo["resultado"]


@then("a justificativa deve nomear a região crítica removida")
def justificativa_nomeia_regiao(mundo):
    criticas = [r.name for r in mundo["gene"].functional_regions if "critical" in r.tags]
    justificativa = mundo["resultado"]["rationale"]
    assert any(nome in justificativa for nome in criticas), justificativa


@then("ela não deve ser aplicada")
def nao_deve_ser_aplicada(mundo):
    assert mundo["resultado"]["applicable"] is False
    # Absent, not empty: a criterion requiring pvs1.strength must report that it
    # could not be evaluated rather than quietly reading a zero.
    assert "strength" not in mundo["resultado"]


@then("o motivo registrado deve ser que a posição de interrupção é desconhecida")
def motivo_posicao_desconhecida(mundo):
    assert mundo["resultado"]["undetermined_reason"] == "ptc_position_unknown"


@then("o motivo registrado deve ser que a consequência não é truncante")
def motivo_nao_truncante(mundo):
    assert mundo["resultado"]["undetermined_reason"] == "consequence_not_truncating"


@then("o motivo registrado deve distinguir isso de uma regra que rodou e não se aplicou")
def motivo_desligada(mundo):
    assert mundo["resultado"]["undetermined_reason"] == "pvs1_disabled_in_spec"
