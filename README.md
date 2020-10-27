# multiple_mann_whitney()

Essa função permite que seja calculado o teste de Mann-Whitney para cada variável disposta em linhas de uma tabela enquanto as amostras estejam nas colunas. 

## Parâmetros

multiple_mann_whitney(**x, y, z, w, p**)  

**x** => dataframe  
**y** => coluna que inicia as amostras do grupo 1  
**z** => coluna que termina as amostras do grupo 1  
**w** => coluna que inicia as amostras do grupo 2  
**p** => coluna que termina as amostras do grupo 2  

## Como utilizar

Para iniciarmos o uso da função, copie todo o texto do arquivo "multiple_mann_whitney.R" e cole em seu ambiente de trabalho (RStudio por exemplo). Em seguida, execute a função para que ela seja armazenada na memória.

Organize seu dataframe como mostra a imagem abaixo. Na primeira coluna coloque as variáveis e nas demais coloque as amostras. Primeiro insira todas amostras de um mesmo grupo e depois insira as de outro grupo. Assim, será garantido que amostras de um mesmo grupo fiquem próximas. Isso será muito importante.

![img1](https://user-images.githubusercontent.com/32198100/97354724-9af5b900-1874-11eb-85aa-5e2b44c088b0.png)

Em nosso exemplo, "dataset" será o nome da variável que conterá nosso dataframe criado acima. Estando o dataframe organizado, basta chamar a função multiple_mann_whitney() com os parâmetros devidamente preenchidos. Procure salvar o resultado em uma nova variável. Em nosso exemplo, executamos a função e salvamos o resultado em uma variável chamada "result" de acordo com o código abaixo:

`result <- multiple_mann_whitney(dataset,2,6,7,11)`

Em caso de dúvida, comparare o código acima com o tópico "Parâmetros".

## Resultado

Será obtido como resultado um dataframe contendo inúmeras colunas. Abaixo, listarei o que cada coluna apresenta.

| Coluna 	| Rótulo 	| Conteúdo 	|
|-	|-	|-	|
| Coluna 1 	| variable 	| lista de variáveis 	|
| Coluna 2 	| w_value 	| valor de W 	|
| Coluna 3 	| p_value 	| valor de P 	|
| Coluna 4 	| p_value_corrected_FDR 	| valor de p corrigido pelo método FDR 	|
| Coluna 5 	| p_value_corrected_Bonferroni 	| valor de p corrigido pelo método Bonferroni 	|
| Coluna 6 	| p_value_corrected_Holm 	| valor de p corrigido pelo método Holm 	|
| Coluna 7 	| p_value_corrected_Hommel 	| valor de p corrigido pelo método Hommel 	|
| Coluna 8 	| p_value_corrected_Hochberg 	| valor de p corrigido pelo método Hochberg 	|
| Coluna 9 	| group1_mean 	| média do grupo 1 	|
| Coluna 10 	| group1_sd 	| desvio-padrão do grupo 1 	|
| Coluna 11 	| group1_median 	| madiana do grupo 1 	|
| Coluna 12 	| group1_IQR 	| Amplitude interquartil do grupo 1 	|
| Coluna 13 	| group1_1st_quantile 	| Primeiro Quartil do grupo 1 	|
| Coluna 14 	| group1_3st_quantile 	| Terceiro Quartil do grupo 1 	|
| Coluna 15 	| group1_max_value 	| Maior valor do grupo 1 	|
| Coluna 16 	| group1_min_value 	| Menor do grupo 1 	|
| Coluna 17 	| group2_mean 	| média do grupo 2 	|
| Coluna 18 	| group2_sd 	| desvio-padrão do grupo 2 	|
| Coluna 19 	| group2_median 	| madiana do grupo 2 	|
| Coluna 20 	| group2_IQR 	| Amplitude interquartil do grupo 2 	|
| Coluna 21 	| group2_1st_quantile 	| Primeiro Quartil do grupo 2 	|
| Coluna 22 	| group2_3st_quantile 	| Terceiro Quartil do grupo 2 	|
| Coluna 23 	| group2_max_value 	| Maior valor do grupo 2 	|
| Coluna 24 	| group2_min_value 	| Menor valor do grupo 2 	|
| Coluna 25 	| group_higt_mean 	| Grupo com maior média 	|
| Coluna 26 	| group_higt_median 	| Grupo com menor média 	|
