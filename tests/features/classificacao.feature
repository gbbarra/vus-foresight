# language: pt
Funcionalidade: Como a evidência acumulada vira uma classificação

  Um curador precisa poder conferir, sem ler código, que a soma da evidência
  produz a classe correta e que a mesma evidência nunca é contada duas vezes.

  Contexto:
    Dado um gene cujo mecanismo de doença é perda de função estabelecida

  Esquema do Cenário: A soma da evidência determina a classe
    Dado uma variante com evidência somando <pontos> pontos
    Quando a classificação for calculada
    Então o resultado deve ser "<classe>"

    Exemplos: as faixas publicadas e as suas fronteiras
      | pontos | classe                   |
      | 10     | Patogênica               |
      | 9      | Provavelmente Patogênica |
      | 6      | Provavelmente Patogênica |
      | 5      | Incerta                  |
      | 0      | Incerta                  |
      | -1     | Provavelmente Benigna    |
      | -6     | Provavelmente Benigna    |
      | -7     | Benigna                  |

  Cenário: Critérios de patogenicidade forte e moderado
    Dado uma variante com o critério PVS1 atribuído
    E o critério PM2_Supporting atribuído
    Quando a classificação for calculada
    Então o resultado deve ser "Provavelmente Patogênica"

  Cenário: A mesma evidência não é contada duas vezes
    Dado uma variante cuja frequência populacional satisfaz dois critérios de frequência
    Quando a classificação for calculada
    Então apenas o critério de frequência mais forte deve contribuir com pontos
    E o critério descartado deve continuar registrado, com o motivo do descarte

  # Atenção, curador: aqui este sistema difere das regras de combinação de 2015,
  # nas quais evidência benigna e patogênica presentes juntas tornavam a variante
  # incerta por definição. No sistema de pontos adotado aqui, as direções se
  # somam. A consequência é que uma variante pode sair patogênica apesar de haver
  # evidência benigna aplicada, e isso é intencional, não um descuido.

  Esquema do Cenário: Direções opostas se somam, em vez de anular o veredito
    Dado uma variante com evidência de patogenicidade somando <a favor> pontos
    E evidência de benignidade somando <contra> pontos
    Quando a classificação for calculada
    Então o resultado deve ser "<classe>"

    Exemplos:
      | a favor | contra | classe                   |
      | 8       | -8     | Incerta                  |
      | 8       | -4     | Incerta                  |
      | 8       | -1     | Provavelmente Patogênica |
      | 12      | -1     | Patogênica               |

  Cenário: Quando a variante fica incerta, o conflito é nomeado como a causa
    Dado uma variante com evidência de patogenicidade somando 8 pontos
    E evidência de benignidade somando -8 pontos
    Quando a classificação for calculada
    Então o mapa deve registrar que o que a trava é o conflito entre as evidências

  Cenário: Um único critério de apoio de cada lado não é um conflito
    Dado uma variante com evidência de patogenicidade somando 1 pontos
    E evidência de benignidade somando -1 pontos
    Quando a classificação for calculada
    Então o mapa não deve registrar conflito
    E deve nomear qual evidência falta para resolvê-la
