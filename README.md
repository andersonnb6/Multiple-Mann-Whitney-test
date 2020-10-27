# Statistic-Functions-R

Neste respositório estão presentes funções em linguagem R que construí durante o Doutorado. Estas funções foram desenvolvidas em grande maioria para resolver pequenos problemas mas que também podem ajudar outras pessoas. Antes de seguirmos a uma breve explicação do que cada função é capaz de realizar, observe as sugestões dispostas neste arquivo.

## Função: multiple_mann_whitney()

Essa função permite que seja calculado o teste de Mann-Whitney para cada variável disposta em linhas de uma tabela enquanto as amostras estejam nas colunas. 

### Parâmetros

multiple_mann_whitney(**x, y, z, w, p**)  

**x** => dataframe  
**y** => coluna que inicia as amostras do grupo 1  
**z** => coluna que termina as amostras do grupo 1  
**w** => coluna que inicia as amostras do grupo 2  
**p** => coluna que termina as amostras do grupo 2  

### Como utilizar

Para iniciarmos o uso da função, copie todo o texto do arquivo "multiple_mann_whitney.R" e cole em seu ambiente de trabalho (RStudio por exemplo). Em seguida, execute a função para que ela seja armazenada na memória.

Organize seu dataframe como mostra a imagem abaixo. Na primeira coluna coloque as variáveis e nas demais coloque as amostras. Primeiro insira todas amostras de um mesmo grupo e depois insira as de outro grupo. Assim, será garantido que amostras de um mesmo grupo fiquem próximas. Isso será muito importante.

![img1](https://user-images.githubusercontent.com/32198100/97354724-9af5b900-1874-11eb-85aa-5e2b44c088b0.png)

Em nosso exemplo, "dataset" será o nome da variável que conterá nosso dataframe criado acima. Estando o dataframe organizado, basta chamar a função multiple_mann_whitney() com os parâmetros devidamente preenchidos. Procure salvar o resultado em uma nova variável. Em nosso exemplo, executamos a função e salvamos o resultado em uma variável chamada "result" de acordo com o código abaixo:

`result <- multiple_mann_whitney(dataset,2,6,7,11)`

Em caso de dúvida, comparare o código acima com o tópico "Parâmetros".

### Resultado

O resultado pode ser armazenado em uma nova variável
