# CLAUDE.md

Instruções permanentes para o Claude neste repositório.
Leia este arquivo antes de qualquer alteração de código.

---

## 1. Regras invioláveis

Estas regras têm precedência sobre qualquer pedido meu feito no calor do momento.
Se um pedido meu conflitar com elas, **pare e me avise** em vez de obedecer.

1. **Nunca altere um teste existente para fazê-lo passar.**
   Se um teste falha, o código está errado — ou o teste está errado e você deve me
   perguntar antes de tocar nele. Editar a asserção para casar com o output é proibido.

2. **Nunca remova, comente, pule (`skip`/`xfail`) ou afrouxe um teste** sem
   autorização explícita minha, mensagem por mensagem.

3. **Nunca reduza um limiar de qualidade** (cobertura mínima, complexidade máxima,
   regras de lint) para fazer o pipeline passar.

4. **Não declare "pronto" ou "funcionando" sem ter executado os testes** e colado a
   saída real. Se você não rodou, diga que não rodou.

5. **Teste antes de código.** Para qualquer comportamento novo: escreva o teste,
   mostre-o falhando, só então implemente. Sem exceção para "mudança pequena".

6. **Nenhum `assert` decorativo.** São proibidos como asserção única:
   `assert result is not None`, `assert result`, `assert len(x) > 0`,
   `assert isinstance(x, dict)`. Toda asserção verifica um valor concreto esperado.

7. **Nada de mock do que está sendo testado.** Mock é para I/O externo (rede, banco,
   sistema de arquivos, API paga). Se você precisou mockar a lógica de negócio para o
   teste passar, o desenho está errado — me avise.

---

## 2. Fluxo de trabalho padrão

Para cada tarefa, siga esta ordem e **pare para confirmação** entre as etapas 2 e 3:

1. **Entender** — releia o código existente antes de escrever. Não presuma a API.
2. **Especificar** — descreva em uma frase o comportamento esperado e os casos de
   borda que vai cobrir. Espere meu OK.
3. **Teste vermelho** — escreva o(s) teste(s) e rode. Cole a saída mostrando a falha.
4. **Implementar** — o mínimo de código para passar. Nada de funcionalidade extra
   "que pode ser útil depois".
5. **Verde** — rode a suíte inteira, não só o teste novo. Cole a saída.
6. **Refatorar** — só com a suíte verde, e rodando novamente ao final.
7. **Gate** — rode lint, tipos e cobertura. Reporte os números.

### Caracterização de código já existente

Escrever teste para código que já existe **não** viola a regra 5 — é teste de
caracterização, e a ordem correta ali é: escreva o teste contra o comportamento que
você **espera**, rode, e se ele falhar você achou um defeito, não um teste errado.
Nunca ajuste a expectativa ao que o código faz sem antes decidir qual dos dois está
errado. Diga explicitamente quando estiver caracterizando em vez de guiando.

---

## 3. Comandos do projeto

```bash
# a suíte inteira, como o CI roda
pytest -m "not validation_study"

# suíte + cobertura de branches
pytest --cov --cov-report=term-missing -m "not validation_study"

# sem nenhum dado de referência (clone limpo): tudo o que não precisa de MANE
pytest -m "not requires_reference and not validation_study"

# só o que exige a referência real fixada em config/genes/*.yaml
pytest -m requires_reference

# um teste específico durante o desenvolvimento
pytest tests/layer3_classification/test_timeline.py::test_a_series_is_diffed_consecutively -x -vv

# lint e formatação
ruff check src tests && ruff format --check src tests

# tipagem estática
mypy

# testes de aceitação (BDD)
pytest tests/features/ -v

# métricas de complexidade e manutenibilidade
radon cc -s -a src/vus_foresight && radon mi src/vus_foresight

# segurança e segredos (o CI reprova em medium+, mais estrito que o §4.6)
bandit -q -r src/vus_foresight -c pyproject.toml --severity-level medium
detect-secrets scan --baseline .secrets.baseline \
  --exclude-files '\.git/' --exclude-files 'reference/.*\.fa$'

# o piso de 95% dos módulos críticos, que o coverage.py não sabe impor
pytest -m "not validation_study" --cov --cov-report=json
python tools/check_module_coverage.py

# os mesmos gates rápidos, antes do push
pre-commit install
pre-commit run --all-files

# teste de mutação (lento — sob demanda; escopo em [tool.mutmut] do pyproject)
pip install -e ".[mutation]"
mutmut run
mutmut results
mutmut show <id>              # o diff de um mutante sobrevivente

# o pipeline em um gene sintético, sem nenhum dado externo
vus-foresight selftest
```

Marcadores registrados em `pyproject.toml`:

- `requires_reference` — precisa do transcrito MANE real, que não é vendorizado no
  repositório mas **está** fixado por `sha256` em `config/genes/*.yaml`.
- `validation_study` — o protocolo da §10. Roda sob demanda, nunca no CI: é um
  estudo de validação, não um teste.

---

## 4. Padrões por tipo de teste

### 4.1 Testes unitários

- Nome descreve o comportamento, não a função:
  `test_retorna_erro_quando_campo_info_ausente`, não `test_parse_1`.
- Estrutura **Arrange / Act / Assert**, com linha em branco separando os blocos.
- Um comportamento por teste. Vários `assert` só se verificarem o mesmo comportamento.
- Sem rede, sem banco, sem I/O real. Use `tmp_path` do pytest para arquivos.
- Para toda função nova, cubra no mínimo: caso feliz, entrada vazia/nula, entrada
  malformada, e o limite (`boundary`) de qualquer comparação numérica.
- Erros esperados se testam com `pytest.raises(TipoDoErro, match="trecho da mensagem")`.
- Use `@pytest.mark.parametrize` em vez de copiar e colar variações do mesmo teste.

#### Organização: por oráculo, não por módulo

Os testes **não** espelham `src/`. Eles são organizados pela pergunta "como este
teste consegue estar errado", que é o que decide se você pode confiar nele:

| diretório | oráculo | falha típica que ele pega |
|---|---|---|
| `tests/layer1_enumeration/` | contagem exata, em forma fechada | variante faltando ou contada duas vezes |
| `tests/layer2_annotation/` | oráculo externo e geometria conhecida | HGVS com aparência correta e coordenada errada |
| `tests/layer3_classification/` | invariantes e snapshot golden | regra aplicada onde não devia |
| `tests/unit/` | o módulo isolado | contrato local: entrada vazia, malformada, limite |

Espelhar `src/` diria onde o código mora; isto diz de que forma o teste pode mentir.
Um teste novo vai para `tests/unit/` só quando é genuinamente local a um módulo e
não pertence a nenhuma das três camadas.

### 4.2 Testes Gherkin / BDD

- Reservados para **regras de negócio que um revisor não-programador precisa validar**.
  Neste projeto isso significa: um geneticista clínico ou curador de VCEP. Não use BDD
  para função utilitária — isso é teste unitário.
- `.feature` em `tests/features/`, escrito em português, com `# language: pt`.
- Steps em `tests/features/steps/`, usando `pytest-bdd`.
- O `.feature` descreve **o quê**, nunca **como**. Proibido citar nome de função,
  classe, endpoint ou estrutura de dados no cenário.
- Escreva o `.feature` primeiro, mostre para mim, e só depois implemente os steps.

```gherkin
# language: pt
Funcionalidade: Classificação de variantes
  Cenário: Critérios de patogenicidade forte e moderado
    Dado uma variante com o critério PVS1 atribuído
    E o critério PM2 atribuído
    Quando a classificação for calculada
    Então o resultado deve ser "Provavelmente Patogênica"
```

### 4.3 Cobertura

- Sempre com `--cov-branch`. Cobertura de linha sozinha esconde `if` não testado.
- Mínimo do projeto: **80%**. Módulos de lógica de decisão: **95%** (lista na §7).
- Cobertura é métrica de **ausência** — mostra o que não foi testado, não valida o que
  foi. Nunca use "100% de cobertura" como argumento de que está correto.
- Ao priorizar, ataque branches de decisão. Ignore boilerplate, `__init__.py`,
  `__repr__` e blocos `if TYPE_CHECKING` — as exclusões estão em `[tool.coverage.report]`.

### 4.4 Teste de mutação

- Rode sob demanda (antes de release ou ao fechar um módulo), nunca no loop de
  desenvolvimento — é lento.
- Para cada mutante sobrevivente: mostre o diff do mutante e escreva o teste que o mata.
- Mutante sobrevivente em código de decisão clínica ou de cálculo é **bloqueante**.
- Se um mutante for genuinamente equivalente (não altera comportamento observável),
  explique por quê em vez de escrever teste artificial.

### 4.5 Métricas de qualidade

- Complexidade ciclomática máxima por função: **10**.
- Função acima de 50 linhas ou com mais de 5 parâmetros: proponha refatoração.
- Ao encontrar violação, **proponha o plano primeiro** e espere meu OK. Não refatore
  código que eu não pedi para tocar.
- Exceção declarada: as funções de comando da CLI. A assinatura de um comando Typer
  **é** a interface de linha de comando — cada parâmetro é uma flag. Contá-los como
  acoplamento mede o tamanho da interface, não a complexidade do código. Se o corpo
  crescer, extraia dele; não achate a assinatura.

### 4.6 Quality gates (CI)

O pipeline deve falhar (exit code ≠ 0) se qualquer etapa não passar:

- [ ] `ruff check` sem erros
- [ ] `ruff format --check` sem diferenças
- [ ] `mypy` sem erros
- [ ] `pytest` com toda a suíte verde
- [ ] cobertura de branches ≥ 80%
- [ ] `bandit` sem achado de severidade alta — **o CI reprova em `medium`**, mais
      estrito que este item e nunca mais frouxo. Os cinco achados `LOW` existentes
      estão revisados e registrados em `[tool.bandit]` do `pyproject`, com o motivo
      de cada um, em vez de silenciados
- [ ] nenhum segredo commitado (`detect-secrets`)
- [ ] cobertura ≥95% em cada módulo crítico da §7 (`tools/check_module_coverage.py`;
      o `coverage.py` só sabe impor um número global, e sem isto um módulo crítico
      apodrece enquanto o total se sustenta nas costas do resto)

Espelhe os mesmos gates em `pre-commit` para pegar antes do push — mas só os
rápidos. A suíte inteira e os pisos de cobertura ficam no CI, onde levar um minuto
não custa nada: um `pre-commit` que as pessoas aprendem a contornar com
`--no-verify` é pior que nenhum, porque transforma a falha do CI em surpresa em vez
de confirmação.

---

## 5. O que NÃO fazer

- ❌ Escrever teste depois do código "para bater a cobertura".
- ❌ Ajustar o teste até passar.
- ❌ `try/except: pass` para silenciar erro que o teste expôs.
- ❌ Adicionar dependência nova sem me perguntar.
- ❌ Alterar arquivo de configuração (`pyproject.toml`, CI, `.pre-commit-config.yaml`)
  no meio de uma tarefa de código, sem avisar.
- ❌ Reescrever módulo inteiro quando eu pedi uma correção pontual.
- ❌ Resumir a saída dos testes. Cole a saída real, inclusive as falhas.

---

## 6. Comunicação

- Responda em português.
- Se um requisito estiver ambíguo, **pergunte antes de implementar** — não escolha por mim.
- Se você discordar de uma decisão minha, diga, com o motivo técnico. Não obedeça em
  silêncio a algo que você acha errado.
- Ao terminar, reporte: testes que passaram/falharam, cobertura antes e depois, e o
  que ficou de fora.

---

## 7. Contexto do projeto

- **Domínio:** genômica clínica. Para **toda variante possível** de um gene — não só as
  observadas — o sistema computa qual classificação ACMG é alcançável apenas com
  evidência independente de paciente, e qual evidência específica faltaria para
  resolver a incerteza. O produto não é uma classificação: é um mapa de onde a
  incerteza mora e do que custaria eliminá-la.

- **Stack:** Python ≥3.11 · Polars (dataframes) · Pydantic v2 (schemas) · DuckDB
  (consulta) · Typer (CLI) · PyYAML (regra-como-dado) · pyarrow (Parquet) ·
  pytest + Hypothesis.

- **Estrutura de diretórios:**

  ```
  src/vus_foresight/
    acmg.py            vocabulário ACMG e o sistema de pontos de Tavtigian
    variant.py         a variante e suas consequências
    gapmap.py          o schema de saída (§6 da spec)
    enumeration/       toda variante possível, por classe
    annotate/          HGVS c./p. e consequência
    genome/            transcrito, coordenadas, MANE, importador da referência
    engine/            avaliador de regras, PVS1, lacuna, equivalência, série temporal
    adapters/          fontes externas (ClinVar, frequência, preditor, splice, MAVE)
    validation.py      o protocolo de validação da §10
    output.py          Parquet determinístico e as agregações
    cli.py             a interface
  tests/               três camadas por oráculo (ver §4.1) + tests/unit + tests/features
  config/genes/        um YAML por gene: coordenadas MANE fixadas por sha256
  config/specs/        um YAML por especificação de VCEP: regra como dado
  results/             tabelas derivadas commitadas pelo workflow gap-map
  docs/                aquisição de dados, referência, protocolo de validação
  ```

- **Módulos críticos** (cobertura ≥95% e mutação limpa). O critério é: um erro aqui
  produz um número que **parece clinicamente interpretável e está errado**.

  | módulo | por quê |
  |---|---|
  | `acmg.py` | a soma de pontos vira a classe. Um limiar errado reclassifica o gene inteiro |
  | `engine/gap.py` | teto, lacuna, conjuntos mínimos suficientes, `blocking_reason` |
  | `engine/pvs1.py` | a árvore de decisão do PVS1; erro aqui dá LP a variante benigna |
  | `engine/evaluator.py` | aplica critério e produz o traço |
  | `engine/predicates.py` | o DSL fechado que lê a spec |
  | `engine/equivalence.py`, `engine/signature.py` | se colapsarem variantes diferentes, a saída inteira mente |
  | `engine/cnv_scoring.py` | escala separada; misturar com a de sequência é erro silencioso |
  | `annotate/*` | **o maior risco do projeto**: HGVS errado com aparência correta |
  | `genome/transcript.py` | fita, éxon, coordenada. Reverter uma fita é biologicamente fatal |
  | `validation.py` | as métricas que dizem se o mapa está certo |

- **Restrições regulatórias / de dado sensível:**
  - **Nenhum dado de paciente entra neste repositório**, em nenhuma forma, nem em
    fixture, nem em teste. Toda evidência é agregada e pública.
  - **Nenhum LLM em nenhuma etapa** do pipeline (§13 da spec). Nem para parsear texto
    de submissão, nem para decidir critério. Se uma tarefa parecer exigir isso, o
    desenho está errado — me avise.
  - **Nenhuma recomendação clínica** é emitida. A saída é um mapa de lacunas.
  - **Determinismo byte-a-byte**: mesma entrada + mesma versão de spec + mesmo snapshot
    = mesma saída. Há teste que compara bytes; qualquer coisa que quebre isso (ordem de
    dicionário, paralelismo, relógio) falha ali.
  - Limiares numéricos ainda não curados vivem sob `verified: false`, e a CLI se recusa
    a escrever um mapa a partir deles sem `--allow-unverified`.
