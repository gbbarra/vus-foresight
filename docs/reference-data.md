# Dados de referência

Nada neste diretório é baixado em tempo de execução do motor. Toda fonte entra por um adaptador com
interface única, declara versão, e é lida de um snapshot local. Reprodutibilidade não é alcançável
se uma fonte pode se mover debaixo de um job em andamento.

Por que a referência não está versionada aqui: coordenadas exônicas e a sequência do transcrito têm
que vir do MANE Select. Inventá-las produz HGVS com aparência correta e conteúdo errado, que a §11
nomeia como a falha de maior risco do sistema. `build_transcript` levanta `ReferenceUnavailable` em
vez de improvisar, e os testes que precisam de uma referência real são marcados
`requires_reference` e pulam limpo.

---

## 1. Transcrito e éxons

Para cada gene são necessários dois arquivos.

**a) FASTA do transcrito spliced**, um único registro, `NM_...` inteiro incluindo UTRs.

```
data/reference/BRCA1_NM_007294.4.fa
data/reference/BRCA2_NM_000059.4.fa
```

Fontes: RefSeq/MANE (`https://ftp.ncbi.nlm.nih.gov/refseq/MANE/`) ou Ensembl REST
(`/sequence/id/ENST...?type=cdna`).

**b) Tabela de éxons**, TSV `label<TAB>start<TAB>end` em **ordem de transcrito** (não genômica).
Para um gene na fita negativa as coordenadas descem; ordenar por coordenada silenciosamente
produziria um transcrito internamente consistente e biologicamente errado.

```
1	43125271	43125364
2	43124017	43124115
...
```

`label` é o nome clínico do éxon, como string. BRCA1 não tem éxon 4 — seus labels vão
`1,2,3,5,...,24`. Nada no código faz aritmética com label.

Então:

```bash
vus-foresight reference import \
  --gene config/genes/BRCA1.yaml \
  --sequence data/reference/BRCA1_NM_007294.4.fa \
  --exons    data/reference/BRCA1_exons.tsv
```

O importador:

- confere que o número e os labels dos éxons batem com o declarado no YAML;
- confere que a soma dos comprimentos exônicos é igual ao comprimento do FASTA;
- **deriva** `cds_start_tx`/`cds_end_tx` procurando a única ORF de exatamente `cds_length` bases que
  começa em ATG, termina em códon de parada e não tem stop interno. Se houver zero ou mais de uma,
  falha em vez de escolher;
- grava o bloco derivado de volta no YAML, junto com `sha256` e `expected_length` da sequência.

A partir daí, uma mudança silenciosa na referência falha a execução em vez de deslocar todas as
coordenadas em uma base.

## 2. Flancos intrônicos

TSV `posição_genômica<TAB>base_na_fita_plus`, cobrindo 50 nt de cada lado de cada junção.

```
data/reference/BRCA1_intronic_flanks.tsv
```

Opcional. Sem ele a enumeração intrônica ainda cobre **todas** as posições — a cardinalidade não
depende da base de referência — mas emite as quatro bases como alternativas e marca cada linha com
`reference_base_unknown`, em vez de adivinhar um alelo de referência.

## 3. Snapshots de evidência

Todos são TSV com cabeçalho, primeira coluna é a chave, e colunas com ponto viram estrutura
aninhada (`gnomad.faf95_popmax` → `frequency.gnomad.faf95_popmax`).

| Adaptador | Namespace | Chave | Colunas típicas |
|---|---|---|---|
| `FrequencyAdapter` | `frequency` | `grch38_pos` | `gnomad.faf95_popmax`, `gnomad.af_popmax`, `abraom.af` |
| `PredictorAdapter` | `predictor` | `grch38_pos` | `bayesdel`, `revel` |
| `SpliceAdapter` | `splice` | `grch38_pos` | `ds_max` |
| `FunctionalAdapter` | `functional` | `hgvs_p` | `classification`, `score`, `dataset` |

`grch38_pos` é `chrom-pos-ref-alt` na **fita plus do genoma**, independentemente da fita do
transcrito.

O ClinVar é a exceção e tem seção própria abaixo: PS1 e PM5 não são respondíveis por busca de chave
única.

### Ausente não é zero

`default_payload` é o que uma fonte afirma sobre uma variante que ela **não lista**, e a distinção
importa mais do que parece:

- o gnomAD não listar uma variante é uma afirmação positiva — procurou-se e não se viu, que é
  exatamente a evidência em que PM2 se apoia. Passe
  `default_payload={"gnomad": {"observed": false, "faf95_popmax": 0.0}}`;
- um ensaio funcional não listar uma variante é o oposto: nada foi medido. Não passe
  `default_payload`, e o critério sai `NOT_EVALUABLE` com o campo ausente nomeado.

### ABraOM

Frequência brasileira é o diferencial do projeto e é reportada em coluna própria
(`frequency.abraom.*`). **Nunca substitui o gnomAD silenciosamente.** Uma variante comum no Brasil e
ausente do gnomAD tem que ficar visível, não ser pontuada como se o mundo tivesse olhado e não
encontrado nada.

### Regiões ensaiadas

`FunctionalAdapter(assayed_regions=...)` recebe triplas `(primeiro_resíduo, último_resíduo, dataset)`.
É o que separa `available_uningested` de `assay_feasible` no relatório de lacunas: PS3 numa região
que uma tela de SGE já cobre é trabalho de engenharia; numa região que ninguém ensaiou é trabalho de
bancada. Colapsar os dois enterra o achado mais acionável do sistema.

BRCA1 e BRCA2 têm datasets de SGE **distintos e independentemente calibrados**. A proveniência é
guardada por variante para que um escore nunca seja atribuído ao ensaio errado.

## 4. Snapshots datados do ClinVar (PS1 e PM5)

O ClinVar não entra por `TableAdapter`. Duas exclusões separam um critério de uma tautologia, e
nenhuma é expressável numa tabela de chave única:

- **PS1 tem que excluir o registro da própria variante.** O critério pergunta se *outra* alteração
  de nucleotídeo que produz a mesma alteração proteica já é patogênica estabelecida. Uma variante
  que é ela mesma patogênica no ClinVar satisfazendo PS1 a partir do próprio registro é circular, e
  fabricaria quatro pontos para toda variante já classificada do gene.
- **PM5 tem que excluir a mesma alteração proteica** (isso é assunto do PS1) e exigir que o registro
  vizinho seja **missense** — "uma alteração *missense* diferente neste resíduo". Um nonsense no
  mesmo códon é observação de PVS1, não de PM5.

Por isso o snapshot é indexado por nucleotídeo, por proteína e por códon, e resolvido por variante.

### Construir o snapshot

```bash
vus-foresight clinvar build \
  --source data/raw/variant_summary_2024-01.txt.gz \
  --out    data/snapshots/clinvar_2024-01-01_BRCA1.tsv \
  --transcript NM_007294.4
```

O parser lê colunas **por nome de cabeçalho**, nunca por posição: o ClinVar adiciona colunas entre
releases, e um parser posicional passa a ler o campo errado silenciosamente na primeira vez que isso
acontece — o que parece mudança de dado, não bug.

O transcrito é casado **com a versão**. Uma classificação feita contra `NM_007294.3` não é evidência
sobre uma coordenada em `NM_007294.4` a menos que alguém tenha conferido que as duas coincidem;
casar frouxamente misturaria transcritos, que é exatamente o que a regra MANE-only existe para
impedir. O comando reporta quantas linhas foram descartadas e por quê, e falha se nada sobrar.

Formato normalizado (ordenado, de modo que dois builds batem byte a byte):

```
hgvs_c	hgvs_p	codon	classification	stars	last_evaluated
```

`stars` vem do `ReviewStatus` do ClinVar. `--clinvar-min-stars` (padrão 1) é o piso que um registro
vizinho precisa vencer para sustentar PS1 ou PM5: uma submissão sem critérios declarados não é um
precedente que alguém citaria.

### O que o mapa deliberadamente não lê

O adaptador reporta `clinvar.self.classification` — a própria opinião do ClinVar sobre a variante —
para proveniência e para o join com o `vus-hindsight`. **Nenhum critério lê esse campo**, e há um
teste que garante isso. Classificar uma variante porque o ClinVar já classificou tornaria o mapa um
espelho, não uma medida.

### Série temporal

Recompute o mapa contra snapshots de datas distintas e faça o diff:

```bash
vus-foresight timeline \
  2018-01-01=out/T2018/gene=BRCA1/gap_map.parquet \
  2024-01-01=out/T2024/gene=BRCA1/gap_map.parquet \
  --out out/transitions.tsv
```

A saída atribui cada transição pelo *tipo* de evidência que a causou —
`neighbour_evidence` quando só critérios semi-intrínsecos se moveram, ou seja, nada novo se
aprendeu sobre aquela variante. Ver [`validation-protocol.md`](validation-protocol.md).

## 5. Fixtures opcionais de teste

Ficam em `tests/fixtures/` e todos os testes que os usam pulam se ausentes.

| Arquivo | Para quê |
|---|---|
| `clinvar_hgvs_oracle.tsv` | oráculo externo primário da camada 2: `gene`, `hgvs_c`, `hgvs_p` curados no transcrito MANE, reproduzidos caractere a caractere |
| `reference_validator_sample.tsv` | oráculo secundário; o teste diagnostica divergência sistemática (convenção) vs esparsa (coordenada) |
| `enigma_three_star.tsv` | concordância em nível de critério e invariante de não-superestimação |
| `hindsight_outcomes.tsv` | protocolo §10 |

## 6. Convenções de HGVS fixadas neste projeto

Uma escolha de convenção errada aparece como divergência **sistemática** contra uma referência
externa, e por isso é diagnosticável. As escolhas são:

- aminoácidos de três letras (`p.Leu100Pro`), como o ClinVar exibe para transcritos MANE;
- sem parênteses em efeito proteico predito — todo efeito aqui é predito, então parentizar todos não
  carregaria informação e garantiria divergência caractere a caractere;
- sinônima como `p.Leu100=`;
- frameshift como `p.Arg100SerfsTer12`, com 12 contando do novo códon até o terminador inclusive;
- deleções e inserções normalizadas para a posição mais 3', com `dup` quando o inserido duplica as
  bases imediatamente anteriores;
- MNV como `delins` sobre o menor intervalo que contém todas as bases alteradas;
- `grch38_pos` estilo VCF na fita plus, com duas exceções documentadas: deleções usam `del` literal
  no campo alt (a base âncora que o VCF exige pode ser intrônica e portanto ausente da sequência do
  transcrito), e o campo é `None` quando a mudança atravessa uma junção — uma variante em espaço
  `c.`, dois registros em espaço genômico.
