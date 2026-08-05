# language: pt
Funcionalidade: A força da evidência de perda de função depende de onde a proteína é interrompida

  A árvore de decisão de perda de função é a regra que mais pesa no sistema, e
  aplicá-la onde não cabe classifica como patogênica uma variante que não é.
  Um curador precisa validar as suas três decisões: se o transcrito é
  degradado, se a região perdida é crítica, e se o mecanismo de doença do gene
  é mesmo perda de função.

  Cenário: Interrupção precoce, com transcrito degradado e mecanismo estabelecido
    Dado um gene cujo mecanismo de doença é perda de função estabelecida
    E uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado
    Quando a evidência de perda de função for avaliada
    Então ela deve ser aplicada com força muito forte

  Cenário: O mesmo achado vale menos quando o mecanismo não está estabelecido
    Dado um gene cujo mecanismo de doença não está estabelecido
    E uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado
    Quando a evidência de perda de função for avaliada
    Então ela deve ser aplicada com força menor que muito forte

  Cenário: Interrupção tardia que ainda remove uma região crítica
    Dado um gene cujo mecanismo de doença é perda de função estabelecida
    E uma variante que interrompe a proteína tarde, escapando da degradação, mas removendo uma região crítica
    Quando a evidência de perda de função for avaliada
    Então ela deve ser aplicada
    E a justificativa deve nomear a região crítica removida

  Cenário: Sem saber onde a proteína é interrompida, a regra não é aplicada
    Dado um gene cujo mecanismo de doença é perda de função estabelecida
    E uma variante truncante cuja posição de interrupção é desconhecida
    Quando a evidência de perda de função for avaliada
    Então ela não deve ser aplicada
    E o motivo registrado deve ser que a posição de interrupção é desconhecida

  Cenário: Uma troca de aminoácido não é evidência de perda de função
    Dado um gene cujo mecanismo de doença é perda de função estabelecida
    E uma variante de sentido trocado sobre a qual nada foi publicado
    Quando a evidência de perda de função for avaliada
    Então ela não deve ser aplicada
    E o motivo registrado deve ser que a consequência não é truncante

  Cenário: Uma especificação pode desligar a regra, e isso fica registrado
    Dado uma especificação que desliga a avaliação de perda de função
    E uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado
    Quando a evidência de perda de função for avaliada
    Então ela não deve ser aplicada
    E o motivo registrado deve distinguir isso de uma regra que rodou e não se aplicou
