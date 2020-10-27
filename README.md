# Statistic-Functions-R

Neste respositório estão presentes funções em linguagem R que construí durante o Doutorado. Estas funções foram desenvolvidas em grande maioria para resolver pequenos problemas mas que também podem ajudar outras pessoas. Antes de seguirmos a uma breve explicação do que cada função é capaz de realizar, observe as sugestões dispostas neste arquivo.

## Sugestões

-> Utilizo o RStudio para realizar análises em R. Contudo, isso não o(a) impede de utilizar a própria interface gráfica do R.

-> Costumo organizar meus dados inicialmente no Excel e posteriormente importo ao R. Portanto, talvez seja necessário carregar um pacote que permite a importação da tabela para o R, sugiro o uso do pacote "readxl".

## Função multiple_mann-whitney

Essa função permite que seja calculado o teste de Mann-Whitney para cada variável disposta em linhas de uma tabela enquanto as amostras estejam nas colunas. 

Para o uso correto da função, organize seu dataframe como mostra a imagem abaixo. Na primeira coluna coloque as variáveis e nas demais coloque as amostras. Primeiro insira todas amostras de um mesmo grupo e depois insira as de outro grupo. Assim, será garantido que amostras de um mesmo grupo fiquem próximas.

![img1](https://user-images.githubusercontent.com/32198100/97354724-9af5b900-1874-11eb-85aa-5e2b44c088b0.png)
