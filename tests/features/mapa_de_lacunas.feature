# language: pt
Funcionalidade: O mapa diz o que falta, não apenas o que se sabe

  O produto não é a classificação: é o mapa de onde a incerteza mora e do que
  custaria removê-la. Um curador precisa poder conferir que, para cada variante
  incerta, o sistema nomeia a evidência que falta e diz se ela já existe
  publicamente ou se ainda precisa ser gerada na bancada.

  Contexto:
    Dado um gene cujo mecanismo de doença é perda de função estabelecida

  Cenário: Toda variante incerta nomeia a evidência que a está travando
    Dado uma variante de sentido trocado sobre a qual nada foi publicado
    Quando a classificação for calculada
    Então o resultado deve ser "Incerta"
    E o mapa deve nomear qual evidência falta para resolvê-la
    E deve informar quantos pontos ainda faltam para cada desfecho

  Cenário: O mapa distingue o que já é público do que precisa ser gerado
    Dado uma variante de sentido trocado sobre a qual nada foi publicado
    Quando a classificação for calculada
    Então cada conjunto de evidência suficiente deve declarar o quanto é viável obtê-lo
    E deve haver pelo menos um conjunto cuja evidência já existe publicamente e só não foi ingerida

  Cenário: Uma variante já resolvida não é reportada como travada
    Dado uma variante que interrompe a proteína cedo o suficiente para o transcrito ser degradado
    Quando a classificação for calculada
    Então o resultado não deve ser "Incerta"
    E o mapa deve registrar que nada a está travando

  Cenário: Nenhum veredito é emitido sem o traço que o sustenta
    Dado uma variante de sentido trocado sobre a qual nada foi publicado
    Quando a classificação for calculada
    Então cada critério aplicado deve declarar a evidência e a fonte que o sustentam
    E cada critério avaliado e não aplicado deve declarar por que não se aplicou

  Cenário: O sistema se recusa a publicar números a partir de limiares não curados
    Dado uma especificação cujos limiares numéricos ainda não foram conferidos
    Quando alguém pedir o mapa desse gene
    Então o sistema deve se recusar a escrevê-lo
    E deve dizer o que é preciso fazer para prosseguir mesmo assim

  Cenário: Variantes de número de cópias não são somadas na mesma escala
    Dado uma deleção que remove o gene inteiro
    Quando a classificação for calculada
    Então o resultado deve ser reportado em uma escala própria
    E deve declarar qual escala foi usada
