# Proposta: as 18 funções acima do limiar de complexidade

CLAUDE.md §4.5 fixa complexidade ciclomática máxima de **10** por função, e manda
**propor o plano primeiro e esperar OK** antes de refatorar. Este documento é essa
proposta. Nenhuma linha de `src/` foi alterada para escrevê-lo.

## O panorama, medido

```
radon cc -s -a src/vus_foresight

420 blocos analisados        média: A (3.58)
  349  A  (1-5)
   53  B  (6-10)
   14  C  (11-20)   ← acima do limiar
    4  D  (21-30)   ← acima do limiar
    0  E/F

radon mi -s src/vus_foresight
  todos os arquivos: A
```

Vale começar pelo que isso **não** diz: 96% dos blocos estão dentro do limiar, a
média é 3,58, e o índice de manutenibilidade é A em todo arquivo. Não há aqui um
problema estrutural de complexidade — há 18 funções específicas, e elas se
dividem em três grupos com respostas diferentes.

## Grupo 1 — Refatorar: complexidade acidental em módulo crítico (4 funções)

Estas estão em módulos que o §7 declara críticos, e a complexidade delas é
acidental: são vários passos independentes num corpo só, não uma decisão
inerentemente ramificada.

| função | arquivo | CC | proposta |
|---|---|---:|---|
| `extract_transcript` | `genome/mane.py:135` | **28** | extrair `_parse_gtf_exons`, `_verify_strand_order`, `_verify_against_fasta`, `_apply_config_labels` |
| `validate` | `validation.py:394` | **22** | extrair uma função por métrica: `_score_resolvability`, `_score_cause`, `_score_direction`, `_collect_temporal` |
| `minimum_sufficient_sets` | `engine/gap.py:152` | **21** | separar a geração de candidatos do subset-sum e da rotulagem de viabilidade |
| `annotate_coding_edits` | `annotate/annotator.py:158` | **19** | separar o cálculo do códon editado da classificação da consequência |

**Por que estas primeiro.** `extract_transcript` é o portão pelo qual toda
coordenada entra no sistema; `validate` produz as métricas que dizem se o mapa
está certo; `minimum_sufficient_sets` produz o conjunto mínimo que é o produto;
`annotate_coding_edits` é o maior risco declarado do projeto. Em todas, o custo de
um erro é um número que parece clinicamente interpretável e está errado.

**Como, com segurança.** As quatro têm cobertura alta e teste de caracterização
(`extract_transcript` tem `test_mane_extract.py`; `minimum_sufficient_sets` e
`annotate_coding_edits` estão sob o golden snapshot). A refatoração é
comportamento-preservante por construção e verificável: suíte verde **mais**
golden byte-idêntico **mais** `selftest` produzindo parquet byte-idêntico. Se
qualquer um dos três mexer, a extração mudou semântica e volta atrás.

**Ordem sugerida:** `extract_transcript` → `validate` → `minimum_sufficient_sets`
→ `annotate_coding_edits`, uma por commit, cada uma com os três verificadores
colados no commit.

## Grupo 2 — Extrair do corpo, sem achatar a assinatura (3 funções)

| função | arquivo | CC | parâmetros |
|---|---|---:|---:|
| `validate_command` | `cli.py:435` | **23** | 11 |
| `map_command` | `cli.py:161` | **20** | 16 |
| `timeline_command` | `cli.py:719` | **17** | 7 |

O §4.5 já declara a exceção: *"A assinatura de um comando Typer **é** a interface
de linha de comando — cada parâmetro é uma flag. Contá-los como acoplamento mede o
tamanho da interface, não a complexidade do código. Se o corpo crescer, extraia
dele; não achate a assinatura."*

A exceção cobre a **contagem de parâmetros**, não a complexidade ciclomática. E a
própria exceção diz o que fazer: extrair do corpo. Proposta concreta:

- `map_command`: extrair `_build_extra_adapters(...)` (a cadeia de cinco `if x is
  not None`) e `_selected_classes(tiers)`. Estimativa: CC 20 → ~8.
- `validate_command`: extrair `_run_curated_mode` e `_run_time_series_mode`, que
  hoje são os dois ramos de um `if` de 60 linhas cada. Estimativa: 23 → ~7.
- `timeline_command`: extrair `_parse_snapshot_args(snapshots)` — que já tem
  validação própria e mereceria teste direto — e `_transitions_to_rows(diffs)`.
  Estimativa: 17 → ~6.

Ganho colateral que vale mais que o número: `_parse_snapshot_args` extraída é
testável sem `CliRunner`, e foi exatamente ali que a Fase 5 achou o defeito do
`label=path`.

## Grupo 3 — Não refatorar; a complexidade é a regra (11 funções)

| função | arquivo | CC |
|---|---|---:|
| `compare_maps` | `engine/timeline.py:205` | 16 |
| `derive_blocking_reason` | `engine/gap.py:213` | 16 |
| `build_synthetic_gene` | `testing.py:98` | 15 |
| `compute_pvs1` | `engine/pvs1.py:60` | 15 |
| `parse_variant_summary` | `adapters/clinvar_import.py:73` | 14 |
| `ClinVarSnapshot.resolve` | `adapters/clinvar.py:215` | 14 |
| `Transcript._check` | `genome/transcript.py:145` | 13 |
| `_evaluate_one` | `engine/evaluator.py:102` | 13 |
| `TranscriptConfig._consistent` | `genome/reference.py:113` | 11 |
| `_preflight` | `genome/importer.py:247` | 11 |
| `annotate_inframe_deletion` | `annotate/indel.py:201` | 11 |

**O argumento.** Em todas estas, os ramos **são** a especificação, e quebrá-las
esconderia a regra em vez de simplificá-la.

`compute_pvs1` é a árvore de decisão do PVS1 — cada ramo é um nó publicado da
árvore de Abou Tayoun, e a Fase 6 acabou de dar teste a cada um deles pelo nome do
motivo. Espalhá-los por quatro funções tornaria mais difícil, não mais fácil,
conferir a árvore contra o artigo.

`Transcript._check`, `TranscriptConfig._consistent` e `_preflight` são sequências
de validações independentes. Cada `if` é uma recusa distinta com sua própria
mensagem e seu próprio teste (Fase 7 escreveu 22 deles). Agrupá-las em
sub-funções trocaria 13 ramos legíveis por 4 funções e um despachante.

`derive_blocking_reason` ranqueia por prioridade declarada na spec; `_evaluate_one`
aplica um critério e produz o traço; `parse_variant_summary` e
`ClinVarSnapshot.resolve` lidam com as formas que um dado externo real assume.

**Contraproposta para este grupo:** em vez de refatorar, registrar a exceção como
dado — um `# radon: allow-complexity <motivo>` ou uma lista em `pyproject` — para
que a *próxima* função a passar de 10 seja um evento visível, e não mais uma numa
lista que ninguém lê. Isso transforma o limiar de aspiração em gate.

## O que eu recomendo

1. **Grupo 2 primeiro** — três funções, ganho de testabilidade imediato, risco
   baixo, e a própria §4.5 já autoriza a forma da mudança.
2. **Grupo 1 em seguida**, uma por commit, com golden e determinismo colados.
3. **Grupo 3 documentado, não tocado**, com a lista virando configuração para que
   o limiar passe a valer daí em diante.

Nada disso está feito. Aguardo o seu OK, e por qual grupo começar.
