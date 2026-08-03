# Como obter os dados externos

Nenhum comando deste documento foi executado neste repositório: o ambiente onde o código foi escrito
tem egresso bloqueado por política para todos os hosts abaixo. Trate as URLs como pontos de entrada
canônicos, não como caminhos de arquivo garantidos — nomes de release mudam a cada versão. O que
está garantido é o **formato de saída** que cada adaptador espera, e isso é verificável offline com
`vus-foresight data check`.

Ordem recomendada: 1 → 2 → 4 → 3 → 5 → 6. As fases 0–3 do projeto só precisam de 1, 2 e 4.

---

## 0. O caminho mais curto: deixar o GitHub buscar

Se a sessão do Claude Code tem egresso restrito — o caso comum — não é preciso liberar nada. Um
runner do GitHub Actions tem rede irrestrita, e os hosts do GitHub estão na allowlist de qualquer
política. Então a aquisição roda lá e a sessão lê o resultado de volta.

```bash
gh workflow run acquire-reference.yml -f clinvar_months=2018-01,2024-01
```

Ou, do celular: aba **Actions** → *acquire-reference* → **Run workflow**.

O workflow (`.github/workflows/acquire-reference.yml`) baixa o release do MANE e os
`variant_summary` do ClinVar, roda **a CLI deste próprio repositório** para extrair e normalizar, e
publica tudo como artifact. Rodar a nossa CLI e não `awk` no YAML é o ponto: o formato de saída é
garantido pela mesma implementação que o motor consome, e a etapa de verificação roda os invariantes
de cardinalidade da §11 contra os transcritos reais **antes** de publicar. Se as coordenadas
estiverem erradas, o workflow falha ali e não publica nada.

Isso também é melhor que um curador rodando `curl` à mão, independentemente de rede: a execução fica
logada, os inputs ficam registrados, e o resultado é reproduzível.

A referência derivada é commitada de volta pelo próprio workflow: são ~19 KB de FASTA para os dois
genes, e eles pertencem ao versionamento ao lado das coordenadas contra as quais estão conferidos.

### E o mapa em si: `gap-map`

Com a referência fixada, o segundo workflow roda a análise inteira no runner:

```bash
gh workflow run gap-map.yml -f snapshot_months=2020-01 -f include_current=true
```

Baixa os releases datados do ClinVar, computa um mapa por gene por data, diffa a série temporal,
roda o estudo da §10 — e commita **só as tabelas derivadas**, em `results/`. Os Parquet por variante
ficam no artifact: 4,6 MB por gene por data, e mudam por inteiro sempre que um snapshot muda, então
não pertencem ao versionamento. As agregações, as transições e as métricas pertencem: são texto,
pequenas, e um diff nelas é um diff em um achado.

O `results/manifest.tsv` registra o `sha256` e a contagem de linhas de cada Parquet de origem, mais
o `run_id`, o commit da CLI e se os limiares da spec estavam curados. É o que permite reproduzir uma
tabela committada em vez de acreditar nela.

Pré-requisito: `acquire-reference` já ter rodado. O workflow começa com `data check` e falha ali se
a referência não estiver fixada, em vez de produzir um mapa a partir de coordenadas ausentes.

O resto deste documento descreve os passos manuais, que continuam valendo quando há rede.

---

## Panorama

| # | Fonte | Para quê | Ordem de grandeza | Precisa de conta |
|---|---|---|---|---|
| 1 | MANE Select / RefSeq | transcrito, éxons | ~100 MB | não |
| 2 | ClinVar | PS1/PM5, série temporal, verdade-terreno da §10 | ~200 MB por snapshot | não |
| 3 | gnomAD v4 | BA1/BS1/PM2 | TB no total, ~MB depois de recortar | não |
| 4 | SpliceAI pré-computado | PP3/BP4 de splicing, PVS1 de sítio canônico | ~30 GB | sim (Illumina) |
| 5 | dbNSFP | BayesDel/REVEL para PP3/BP4 | dezenas de GB | não |
| 6 | MaveDB + suplementares de SGE | PS3/BS3 | KB a MB | não |
| 7 | ABraOM | frequência brasileira | ~GB | não |
| — | Spec do ENIGMA (cspec) | curar os limiares | documento | não |

A regra geral: **recorte por região antes de guardar.** BRCA1 é ~81 kb em chr17 e BRCA2 ~84 kb em
chr13. Baixar gnomAD ou dbNSFP inteiros para extrair dois genes é desperdício de ordens de grandeza;
todos os arquivos grandes são indexados por tabix.

---

## 1. MANE Select — transcrito e éxons

Sem isto nada roda: `build_transcript` levanta `ReferenceUnavailable` e todos os testes marcados
`requires_reference` pulam.

**Entrada:** `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/` → escolha um diretório
`release_X.Y/`. Os arquivos relevantes de cada release são o FASTA de RNA do RefSeq, o GTF genômico,
e um `summary.txt` que lista o par RefSeq/Ensembl de cada gene.

```bash
mkdir -p data/reference
BASE=https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4
curl -fSL -o data/reference/MANE.GRCh38.v1.4.refseq_rna.fna.gz \
  $BASE/MANE.GRCh38.v1.4.refseq_rna.fna.gz
curl -fSL -o data/reference/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz \
  $BASE/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz
```

**Não faça as duas coisas à mão.** Um comando faz o par inteiro, e faz as verificações que um
pipeline de `awk` não faz:

```bash
vus-foresight reference from-mane \
  --gene config/genes/BRCA1.yaml \
  --gtf   data/reference/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz \
  --fasta data/reference/MANE.GRCh38.v1.4.refseq_rna.fna.gz
```

Use o GTF **refseq**_genomic, não o ensembl_genomic: o primeiro é chaveado por acessos `NM_`, que é
o que os gene configs nomeiam; o segundo usa identificadores `ENST` e não acharia nada.

Três armadilhas que esse comando resolve, e as três produzem um transcrito internamente consistente
e biologicamente errado:

- **ordem de transcrito, não de coordenada.** O `exon_number` do GTF já está em ordem de transcrito;
  ordenar por `start` inverte um gene da fita negativa. O comando toma a ordem do `exon_number` e
  então **confere contra a fita** — descendente para `-`, ascendente para `+`.
- **os labels não são o `exon_number`.** BRCA1 não tem éxon 4: a numeração clínica vai
  `1,2,3,5,...,24` enquanto o `exon_number` do GTF vai `1..23`. Os labels vêm do gene config, e o
  comando avisa quando os dois divergem.
- **GTF e FASTA de releases diferentes.** Cada um plausível sozinho, juntos deslocam todas as
  coordenadas. O comando compara a soma dos comprimentos exônicos com o tamanho do registro FASTA e
  recusa se diferirem.

Ele já roda o `reference import` em seguida, que **deriva** `cds_start_tx`/`cds_end_tx` e grava o
`sha256`. Para importar de arquivos que você mesmo preparou, o passo continua disponível:

```bash
vus-foresight reference import \
  --gene config/genes/BRCA1.yaml \
  --sequence data/reference/BRCA1_NM_007294.4.fa \
  --exons    data/reference/BRCA1_exons.tsv
```

O importador **deriva** `cds_start_tx`/`cds_end_tx` procurando a única ORF de exatamente 5.592 bases
que começa em ATG, termina em códon de parada e não tem stop interno. Se houver zero ou duas, ele
falha em vez de escolher. Ele também grava `sha256` de volta no YAML, de modo que uma troca
silenciosa de release passa a falhar a execução em vez de deslocar tudo em uma base.

**Verificação:** `vus-foresight enumerate --gene config/genes/BRCA1.yaml` tem que imprimir
exatamente `16.776` SNVs codificantes e `100.656` MNVs. Qualquer outro número é erro de coordenada.

**Flancos intrônicos** (opcional) precisam do FASTA genômico do GRCh38, e saem com `samtools faidx`
nas 50 bases de cada lado de cada junção, gravadas como `posição<TAB>base_na_fita_plus`.

## 2. ClinVar — PS1/PM5 e a verdade-terreno

**Entrada:** `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`. O `variant_summary.txt.gz`
na raiz é o release corrente; o subdiretório `archive/` guarda releases mensais antigos, que é o que
a série temporal precisa.

```bash
mkdir -p data/raw data/snapshots
for month in 2018-01 2020-01 2022-01 2024-01; do
  curl -fSL -o data/raw/variant_summary_$month.txt.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_$month.txt.gz
  vus-foresight clinvar build \
    --source data/raw/variant_summary_$month.txt.gz \
    --out    data/snapshots/clinvar_${month}-01_BRCA1.tsv \
    --gene   config/genes/BRCA1.yaml
done
```

O comando reporta o que descartou e por quê, e falha se nada sobrar. Confira o número de linhas
mantidas: para BRCA1 num release recente espera-se dezenas de milhares de registros no transcrito
MANE. Se vier zero, quase sempre é a versão do transcrito — releases antigos do ClinVar usam
`NM_007294.3`, e o casamento é **com a versão** de propósito. Nesse caso, ou construa o snapshot
antigo com `--transcript NM_007294.3`, ou aceite que aquele intervalo não é comparável; misturar as
duas silenciosamente é o que a regra MANE-only existe para impedir.

**Uma nota sobre o que o ClinVar não carrega.** A métrica 2 da §10, na formulação original, quer o
tipo de evidência citado no registro de submissão. Isso não é campo estruturado no
`variant_summary`; está em prosa livre no `submission_summary.txt` e no XML dos SCVs. Como este
projeto não usa LLM em nenhuma etapa, extrair isso significaria um classificador de palavras-chave
que ninguém consegue calibrar. Por isso a métrica 2 aqui é computada de outro jeito — pela causa
que de fato se moveu nas fontes do próprio pipeline — e o `submission_summary.txt` fica como
curadoria opcional, alimentando a coluna `evidence_type` da tabela de desfechos, que tem
precedência quando presente.

## 3. gnomAD v4 — frequência

**Entrada:** `https://gnomad.broadinstitute.org/downloads`. Os VCFs de sítios estão em bucket
público do Google Cloud e são indexados por tabix, então recorte antes de guardar.

```bash
# coordenadas GRCh38 aproximadas; use as do gene YAML depois do import
tabix -h <URL_do_vcf_de_sitios_chr17> chr17:43,000,000-43,180,000 > brca1_gnomad.vcf
```

O que extrair: a **filtering allele frequency** (FAF95) do popmax, não a AF bruta. É o que os
critérios de frequência do VCEP usam, e é a diferença entre "raro no mundo" e "raro na amostra que
por acaso foi sequenciada". Formato de saída:

```
grch38_pos	gnomad.faf95_popmax	gnomad.af_popmax	abraom.af
chr17-43094500-C-T	0.0	0.0	
```

**Ausente não é zero, mas neste caso é uma afirmação.** O gnomAD não listar uma variante significa
que se procurou e não se viu — exatamente a evidência de PM2. Por isso o adaptador de frequência
recebe `default_payload={"gnomad": {"observed": false, "faf95_popmax": 0.0}}`. Ensaios funcionais
são o oposto e não recebem default: ausência ali significa que nada foi medido.

## 4. SpliceAI pré-computado

**Entrada:** distribuição da Illumina, via BaseSpace, com registro. Existem espelhos acadêmicos;
confira a licença antes de usar.

Cobre todos os SNVs possíveis do genoma, que é exatamente o que este projeto precisa — a enumeração
é de toda variante possível, não das observadas. Recorte por região como no gnomAD, e extraia o
máximo dos quatro delta scores:

```
grch38_pos	ds_max
chr17-43094500-C-T	0.02
```

Sem isto, PVS1 em sítio de splice canônico sai `NOT_EVALUABLE` com
`undetermined_reason=splice_prediction_missing`, e o relatório `available_uningested` mostra a
tabela faltante como lacuna acionável. Esse é o comportamento pretendido, não um erro.

## 5. dbNSFP — preditores

**Entrada:** a página de distribuição do dbNSFP (`sites.google.com/site/jpopgen/dbNSFP`) e seus
espelhos. Arquivo enorme, dividido por cromossomo, indexado por tabix.

Extraia **o preditor que a spec do VCEP determina**, com os limiares dela. Carregar REVEL e pontuar
com limiar de BayesDel é erro de especificação, e é pego pelo `requires:` do critério — o campo
declarado simplesmente não existe no contexto e PP3 sai `NOT_EVALUABLE` em vez de disparar errado.

```
grch38_pos	bayesdel	revel
chr17-43094500-C-T	0.41	0.62
```

## 6. MaveDB e os artigos de SGE — PS3/BS3

**Entrada:** `https://mavedb.org` para os datasets versionados, mais as tabelas suplementares dos
artigos de saturation genome editing. BRCA1 e BRCA2 têm datasets **distintos e independentemente
calibrados**; a proveniência é guardada por variante para que um escore nunca seja atribuído ao
ensaio errado.

```
hgvs_p	classification	score	dataset
p.Asp1692Asn	abnormal	-2.14	BRCA1-SGE-Findlay2018
```

E declare a **cobertura** no gene YAML:

```yaml
mave_datasets:
  - {name: BRCA1-SGE-Findlay2018, start_aa: 1, end_aa: 109}
  - {name: BRCA1-SGE-Findlay2018, start_aa: 1646, end_aa: 1859}
```

Isso é o que separa `available_uningested` de `assay_feasible` no relatório de lacunas: PS3 numa
região que uma tela já cobre é trabalho de engenharia; numa região que ninguém ensaiou é trabalho de
bancada. Uma lista vazia subdeclara quanto já é respondível.

## 7. ABraOM — frequência brasileira

**Entrada:** `https://abraom.ib.usp.br/`. É o diferencial do projeto e entra em coluna própria,
`frequency.abraom.*`. **Nunca substitui o gnomAD silenciosamente.** Uma variante comum no Brasil e
ausente do gnomAD tem que ficar visível, não ser pontuada como se o mundo tivesse olhado e não
encontrado nada.

Note que nenhum critério da spec lê `frequency.abraom.*` hoje — deliberadamente, porque o VCEP não
especifica limiar para ela. Para usá-la é preciso decidir e escrever esse limiar no YAML, o que é
uma decisão de curadoria, não de engenharia.

## Curar os limiares do ENIGMA

`https://cspec.genome.network/cspec/ui/svi/` publica as especificações dos VCEPs. O checklist está
no cabeçalho de `config/specs/enigma_brca_v1.1.0.yaml`; quando terminar, vire `verified: true` e a
CLI passa a escrever o mapa sem `--allow-unverified`.

Enquanto isso, o `toy_v0.1.0.yaml` é internamente consistente e é o que a suíte usa, então a
curadoria não bloqueia nenhum teste.

## Verificar o que foi materializado

```bash
vus-foresight data check --gene config/genes/BRCA1.yaml 
```

Reporta, por fonte: presente ou não, versão declarada, número de registros, e quanto da região do
transcrito está coberta. É o comando para rodar depois de cada passo acima, antes de gastar tempo
computando um mapa sobre dados incompletos.
