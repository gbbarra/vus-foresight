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
    E o critério PM2 atribuído
    Quando a classificação for calculada
    Então o resultado deve ser "Provavelmente Patogênica"

  Cenário: A mesma evidência não é contada duas vezes
    Dado uma variante cuja frequência populacional satisfaz dois critérios de frequência
    Quando a classificação for calculada
    Então apenas o critério de frequência mais forte deve contribuir com pontos
    E o critério descartado deve continuar registrado, com o motivo do descarte

  Cenário: Evidência a favor e contra ao mesmo tempo não produz um veredito
    Dado uma variante com evidência forte de patogenicidade
    E evidência forte de benignidade
    Quando a classificação for calculada
    Então o resultado deve ser "Incerta"
    E o mapa deve registrar que há conflito de direção
