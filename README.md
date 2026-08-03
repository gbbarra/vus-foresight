# vus-foresight

> Para **toda variante possível** em um gene — não apenas as observadas — computar qual
> classificação ACMG é alcançável somente com evidência independente de paciente, e **qual
> evidência específica faltaria** para resolver a incerteza.

O produto não é uma classificação. É um mapa de onde a incerteza mora e do que custaria eliminá-la.

Gêmeo prospectivo do `vus-hindsight`, cujos snapshots datados consome.
BRCA1/BRCA2 primeiro, arquitetura gene-agnóstica.

---

## Estado atual

| Fase | Escopo | Estado |
|---|---|---|
| 0 | Enumeração, anotação, schema, motor de pontos, testes | **implementada** |
| 1 | Critérios intrínsecos, PVS1, classes de equivalência, traço completo | **implementada**; concordância com o ENIGMA pendente de dados curados |
| 2 | Teto, conjuntos mínimos suficientes, `blocking_reason`, relatório `available_uningested` | **implementada** |
| 3 | PS1/PM5 contra snapshot datado do ClinVar | adaptador e critérios implementados; falta o snapshot |
| 4 | Protocolo de validação §10 | harness implementado, roda sob demanda |
| 5 | Segundo gene | testado com gene sintético; adicionar ATM é um YAML |

Duas coisas ainda **não** estão no repositório, e é deliberado:

1. **Coordenadas e sequência de BRCA1/BRCA2.** MANE Select é a única fonte de verdade, e este
   ambiente não tem acesso de rede ao Ensembl/NCBI. Inventar coordenadas produziria HGVS com
   aparência correta e conteúdo errado — a falha que a §11 nomeia como a de maior risco de todo o
   sistema. O carregador **falha explicitamente** em vez de improvisar. Veja
   [`docs/reference-data.md`](docs/reference-data.md).
2. **Limiares numéricos curados** da spec ENIGMA. A estrutura do YAML está completa e testada; os
   números são placeholders marcados `verified: false`, e a CLI se recusa a escrever um mapa a
   partir deles sem `--allow-unverified`.

Tudo o que não depende desses dois roda hoje, sobre genes sintéticos, com 136 testes.

```bash
pip install -e ".[dev]"
pytest                        # três camadas, sem nenhum dado de referência
vus-foresight selftest        # pipeline completo em um gene que não existe
```

---

## O modelo de pontos é a simplificação central

Sistema de Tavtigian (ClinGen SVI), não a tabela combinatória de 2015:

| Força | Patogênico | Benigno |
|---|---|---|
| Supporting | +1 | −1 |
| Moderate | +2 | — |
| Strong | +4 | −4 |
| Very Strong | +8 | −8 (stand-alone) |

**P ≥ 10 · LP 6–9 · VUS 0–5 · LB −1 a −6 · B ≤ −7**

Combinar critérios vira soma; a lacuna vira subtração. "Faltam 4 pontos para LP" é computável, e o
conjunto mínimo suficiente vira uma soma de subconjuntos trivialmente pequena
([`engine/gap.py`](src/vus_foresight/engine/gap.py)). Com a tabela de 2015 essa pergunta não tem
forma fechada.

---

## Arquitetura

```
config/genes/*.yaml      transcrito, domínios, mecanismo de LoF   <- conhecimento de domínio
config/specs/*.yaml      critérios, limiares, modulações de força <- conhecimento de domínio
src/vus_foresight/
  genome/        sequência, geometria do transcrito, importação da referência
  enumeration/   toda variante possível, por classe (Tier 1 e 2)
  annotate/      HGVS c./p., consequência, indels
  engine/        contexto, DSL de regras, PVS1, avaliador, lacuna, equivalência, CNV
  adapters/      fontes de dados versionadas, sempre offline
  output.py      Parquet determinístico e as agregações
  validation.py  o protocolo §10
```

O motor recebe `(variante, gene_config, spec)` e **não contém nenhuma constante de domínio**. Isso
não é uma promessa: [`test_domain_isolation.py`](tests/layer3_classification/test_domain_isolation.py)
faz o parse do código executável do motor, remove docstrings, e falha se `BRCA1`, `ENIGMA`,
`gnomAD` ou similares aparecerem.

### Escopo de enumeração

| Classe | Cardinalidade | BRCA1 | BRCA2 |
|---|---|---|---|
| SNV codificante | `3 × len(CDS)` | 16.776 | 30.771 |
| SNV intrônico (±50 nt de cada junção) | por posição, deduplicado | ~6.600 | ~7.800 |
| MNV intra-códon | `54 × n_códons` | 100.656 | 184.626 |
| Classe de frameshift (por posição do PTC) | ≤ `len(proteína)` | ≤ 1.863 | ≤ 3.418 |
| Deleção in-frame de 1–2 códons | por posição, 3'-normalizada | ~11k | ~20k |
| CNV exônico | `n(n+1)/2` por direção | 276 | 378 |

Fora de escopo, documentado e não silenciosamente omitido: inserções arbitrárias (crescem como
`4^k`), variantes profundamente intrônicas fora da janela, rearranjos complexos, e combinações em
*loci* distintos — estas dependem de fase, que é propriedade do paciente e não da variante.

**Nota sobre classes de frameshift.** A §3 pede uma classe por posição de códon. Algumas posições
não são o *primeiro* stop de nenhum indel mínimo de 1–2 nt; o enumerador reporta essas separadamente
em `FrameshiftEnumeration.unreachable` em vez de emitir classes vazias, e o invariante testado é
`len(classes) + len(unreachable) == len(proteína)`.

### Distinção intrínseco/extrínseco

É a espinha dorsal. Todo critério carrega `evidence_class`:

- `intrinsic` — computável sem paciente algum (PVS1, frequência, PP3/BP4, PS3/BS3, PM1, PM4, BP7)
- `semi_intrinsic` — depende do estado de um banco público numa data (PS1, PM5); sempre computado
  contra snapshot datado, e é o gancho direto para o `vus-hindsight`
- `extrinsic` — exige observação em paciente, família ou coorte; **nunca aplicado por este sistema**,
  apenas listado como o que faltaria

### Traço completo, nunca só o veredito

Toda avaliação retorna cada critério testado, aplicado ou não, com a razão. `MISSING` e `False` são
estados distintos e permanecem distintos: um critério cujo dado não foi carregado sai como
`NOT_EVALUABLE` com os campos ausentes nomeados, jamais como `NOT_MET`. É desse registro que sai o
relatório `available_uningested`.

---

## Uso

```bash
# quantas variantes existem
vus-foresight enumerate --gene config/genes/BRCA1.yaml --data-root data

# o mapa
vus-foresight map --gene config/genes/BRCA1.yaml --data-root data \
  --frequency data/snapshots/gnomad_v4_brca1.tsv \
  --predictor data/snapshots/dbnsfp_brca1.tsv \
  --splice    data/snapshots/spliceai_brca1.tsv \
  --functional data/snapshots/sge_brca1.tsv \
  --clinvar   data/snapshots/clinvar_2026-01-01.tsv --clinvar-date 2026-01-01 \
  --out-dir out

# as agregações
vus-foresight report out/gene=BRCA1/gap_map.parquet

# o estudo de validação (§10)
vus-foresight validate out/gene=BRCA1/gap_map.parquet outcomes.tsv --reference-date 2020-01-01
```

Consultável direto por DuckDB:

```sql
SELECT blocking_reason, count(*)
FROM 'out/gene=BRCA2/gap_map.parquet'
WHERE class_current = 'VUS'
GROUP BY 1 ORDER BY 2 DESC;
```

`blocking_reason` é enum exatamente para que isso seja um `GROUP BY`.

---

## Testes — três camadas, três oráculos

Cada camada falha de um jeito distinto. A camada 2 é a de maior risco: erro nela produz HGVS errado
com aparência correta.

**Camada 1 — enumeração.** Oráculo: contagem exata. Cardinalidade em forma fechada, completude
contra um produto cartesiano computado independentemente, ausência de duplicata e de no-op, 54 MNVs
por códon, `n(n+1)/2` intervalos de éxon, e — para frameshift — verificação de que os indels
listados de fato produzem aquele PTC **por tradução, não por metadado**. Property-based com
Hypothesis sobre posições aleatórias do CDS.

**Camada 2 — anotação.** Oráculos que não compartilham implementação com o código testado:

- *round-trip* genômico → `c.` → genômico é identidade;
- *tradução direta*: aplicar a variante ao CDS, traduzir, comparar com o `p.` predito;
- *simetria de fita*: os dois fixtures compartilham a sequência do transcrito, um na fita `+` e
  outro na `−`; todo `c.` e todo `p.` têm que bater caractere a caractere enquanto as coordenadas
  genômicas diferem. Bug de complemento reverso aparece em um e não no outro — o mesmo controle de
  graça que BRCA1(−)/BRCA2(+) dá;
- *fixtures de fronteira*: primeira e última base do CDS, códon de iniciação, terminador e
  stop-loss, as quatro posições canônicas de cada junção, limite de NMD, descontinuidade de
  numeração de éxons;
- *MNV*: códons onde duas substituições isoladas dão X e Y mas juntas dão um terceiro aminoácido Z.
  O sistema tem que retornar Z. Este teste pega a falha de compor consequências por posição em vez
  de traduzir o códon inteiro — o erro mais provável de toda a camada.

Os oráculos externos (ClinVar curado, VariantValidator/Mutalyzer) estão implementados em
[`test_external_oracle.py`](tests/layer2_annotation/test_external_oracle.py) e pulam limpo enquanto
os fixtures não existirem; o teste do validador **diagnostica** divergência sistemática (erro de
convenção) versus esparsa (erro de coordenada), porque o conserto é diferente.

**Camada 3 — classificação.** Monotonicidade por fuzzing, exclusividade mútua, aritmética da lacuna
consistente com a classificação, minimalidade por inclusão dos conjuntos suficientes, exclusão de
evidência intratável, determinismo (Parquet byte-idêntico entre duas execuções), isolamento de
domínio, e regressão por snapshot dourado versionado — qualquer diff exige justificativa no PR.

A concordância em nível de critério com registros três estrelas do ENIGMA e o invariante de
não-superestimação estão como harness e pulam sem os dados curados.

O protocolo da §10 é **estudo de validação**, não teste unitário, e não roda em CI.

---

## Não-objetivos

- **Não** classifica variantes de pacientes reais. Este sistema não vê dado de paciente.
- **Não** substitui o classificador VCEP. Compartilha o motor, responde outra pergunta.
- **Não** usa LLM. Em nenhuma etapa.
- **Não** emite recomendação clínica. A saída é insumo de curadoria e priorização de pesquisa.
