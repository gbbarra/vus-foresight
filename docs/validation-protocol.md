# Protocolo de validação (§10)

O mapa faz uma afirmação falsificável, e é isso que o transforma de taxonomia em ciência.

**Não é teste unitário.** As camadas 1–3 rodam em CI a cada commit; isto roda sob demanda, contra
os snapshots datados do `vus-hindsight`.

## Protocolo

1. Computar o mapa com o estado de evidência da data **T** — mesmos snapshots de gnomAD, dbNSFP,
   SpliceAI, MAVE e ClinVar que existiam em T. É o `--clinvar-date` e as versões declaradas nos
   adaptadores que tornam isso auditável.
2. Selecionar as variantes classificadas como VUS em T.
3. Extrair do ClinVar quais delas foram reclassificadas até **T+n**, com que classe, com que tipo de
   evidência no registro de submissão, e em que data.

```bash
vus-foresight map --gene config/genes/BRCA1.yaml \
  --clinvar data/snapshots/clinvar_2020-01-01.tsv --clinvar-date 2020-01-01 \
  --out-dir out/T2020

vus-foresight validate out/T2020/gene=BRCA1/gap_map.parquet \
  data/hindsight/brca1_outcomes_2020_2026.tsv \
  --reference-date 2020-01-01
```

### Formato da tabela de desfechos

TSV com cabeçalho:

| coluna | conteúdo |
|---|---|
| `variant_id` | `NM_007294.4:c.5074G>A` — transcrito mais HGVS `c.`, exatamente como o mapa emite |
| `class_at_t` | `P` `LP` `VUS` `LB` `B` |
| `class_at_t_plus_n` | idem |
| `evidence_type` | `functional`, `segregation`, `case_control`, `population_frequency`, `computational` |
| `resolved_on` | data ISO, opcional; necessária apenas para a métrica 4 |

## As quatro métricas

**1. Recall de resolubilidade.** Entre as que foram resolvidas, qual fração o mapa marcava como
resolvível? "Resolvível" aqui não é apenas ter lacuna finita: é ter pelo menos um conjunto mínimo
suficiente cuja evidência alguém consegue de fato obter. Uma variante cujos únicos conjuntos
suficientes são `intractable` conta como **não** resolvível — dizer "bastaria um evento de novo em
câncer hereditário de início adulto" não é um plano de pesquisa.

As que o mapa não antecipou são listadas nominalmente, nunca resumidas em um número.

**2. Acurácia da causa.** O `blocking_reason` previsto corresponde ao tipo de evidência que de fato
apareceu no registro de submissão? A saída inclui a matriz de confusão completa
previsto × observado, porque a estrutura do erro é mais informativa que a taxa.

É a métrica mais interessante e a mais publicável: ela testa se o sistema entende *por que* uma
variante está incerta, não apenas *que* está.

**3. Direção.** O mapa apontava para o lado certo? A direção prevista é o alvo do conjunto
suficiente mais barato — para onde o mapa diria a alguém para ir primeiro.

**4. Calibração temporal.** Variantes com lacuna menor foram resolvidas antes? Correlação de posto
de Spearman entre a magnitude da lacuna e os dias até a resolução; **negativa** significa
calibração correta. Implementada sem scipy, para não trazer a dependência por causa de uma métrica.

## Interpretando um resultado ruim

Cada métrica falha por um motivo diferente e o conserto é diferente:

- **recall baixo** — o mapa está declarando variantes irresolvíveis que o campo resolveu. Suspeitar
  primeiro das tags de viabilidade: `intractable` está sendo aplicado onde não deveria, ou
  `max_cardinality` está cortando conjuntos legítimos.
- **causa errada, direção certa** — o sistema acerta pelo motivo errado. Olhar a matriz de confusão:
  se `MISSING_FUNCTIONAL` está sendo previsto onde apareceu segregação, a prioridade em
  `blocking.priority` está mal ordenada para este gene, ou PS3/BS3 estão `NOT_EVALUABLE` cedo demais.
- **direção errada** — problema no lado intrínseco, não na lacuna. Provavelmente calibração de
  preditor ou limiar de frequência. Rodar antes a concordância em nível de critério da camada 3.
- **correlação temporal nula** — pode ser o mapa e pode ser o campo. O tempo até a reclassificação
  depende de quem estava olhando, não só de quanta evidência faltava; um `rho` próximo de zero é
  informativo mas não é, sozinho, uma refutação.

## A série temporal é o passo anterior

O protocolo acima compara o mapa com o que o *campo* fez. O `timeline` compara o mapa consigo mesmo
em datas distintas, e responde uma pergunta diferente: quantas variantes o próprio mapa moveu, e por
quê.

```bash
vus-foresight timeline \
  2018-01-01=out/T2018/gene=BRCA1/gap_map.parquet \
  2020-01-01=out/T2020/gene=BRCA1/gap_map.parquet \
  2024-01-01=out/T2024/gene=BRCA1/gap_map.parquet \
  --out out/transitions.tsv
```

Diferenças consecutivas, não todas contra a primeira: a pergunta da §4 é *quando* uma variante se
moveu e sobre o quê, e colapsar cinco anos num único antes/depois perde exatamente isso.

Cada transição é atribuída pelo tipo de evidência que a causou:

| causa | significado |
|---|---|
| `neighbour_evidence` | só critérios semi-intrínsecos se moveram — **nada novo se aprendeu sobre esta variante**; um vizinho foi classificado |
| `intrinsic_evidence` | frequência, preditor ou ensaio novos |
| `mixed_evidence` | os dois no mesmo intervalo; não são separáveis |
| `strength_only` | os mesmos critérios, em força diferente |
| `spec_change` | nenhum critério se moveu e a classe mudou — foram as regras, não a evidência |

`neighbour_evidence` cruzado com "saiu de VUS" é o achado da §4 em forma computável. Vale notar a
escala: uma classificação não move uma variante, move o códon inteiro — as que alcançam a mesma
alteração proteica ganham PS1, o resto do códon ganha PM5.

## Contrafactual barato

Para responder "quanto valeria ingerir este dataset", rode o mapa duas vezes, uma com o adaptador
e outra com ele desligado, e faça o diff dos `blocking_reason`. É para isso que `NullAdapter`
existe.
