# language: pt
Funcionalidade: Uma variante pode mudar de classe sem evidência sobre ela mesma

  Esta é a observação que justifica separar a evidência que depende do estado
  de um banco público daquela que depende só da variante. Um curador precisa
  validá-la diretamente: se ela não for verdadeira, a distinção é decorativa.

  Contexto:
    Dado um gene cujo mecanismo de doença é perda de função estabelecida
    E uma variante de sentido trocado sobre a qual nada foi publicado

  Cenário: Classificar uma vizinha move a variante, sem nenhum dado novo sobre ela
    Dado que nenhuma outra variante do mesmo códon está classificada
    Quando outra variante do mesmo códon for classificada como patogênica
    E a classificação for recalculada
    Então a variante deve ganhar evidência
    E essa evidência deve ser atribuída à classificação da vizinha
    E nenhuma evidência nova sobre a própria variante deve ter sido usada

  Cenário: Uma classificação move o códon inteiro, não uma variante isolada
    Quando outra variante do mesmo códon for classificada como patogênica
    E a classificação for recalculada
    Então todas as variantes daquele códon devem ter se movido
    E nenhuma variante de outro códon deve ter se movido

  Cenário: A própria variante não é evidência de si mesma
    Dado que a própria variante já está classificada como patogênica no banco público
    E que nenhuma outra variante do mesmo códon está classificada
    Quando a classificação for calculada
    Então a variante não deve ganhar evidência de vizinho
    E o mapa deve continuar dizendo que falta evidência para resolvê-la

  Cenário: A classificação publicada da própria variante nunca entra no cálculo
    Dado que a própria variante já está classificada como patogênica no banco público
    Quando a classificação for calculada
    Então nenhum critério deve ter consultado a classificação publicada da própria variante

  Esquema do Cenário: Uma submissão sem revisão suficiente não move ninguém
    Dado que nenhuma outra variante do mesmo códon está classificada
    Quando outra variante do mesmo códon for classificada como patogênica com <estrelas> estrela(s)
    E a classificação for recalculada
    Então a variante <deve ou não> ter ganhado evidência

    Exemplos:
      | estrelas | deve ou não |
      | 0        | não deve    |
      | 1        | deve        |
      | 3        | deve        |
