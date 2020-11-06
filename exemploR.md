# R Markdown
## Testando Markdown com R


```R
mann_whitney_analysis <- function(x,y,z,w,p) {
  # Obs.:  usa coluna 1 para distinguir as linhas no arquivo de saida
  # x: tabela (dataframe)
  # y: coluna que inicia as amostras com periodontitis
  # z: coluna que termina as amostras com periodontitis
  # w: coluna que inicia as amostras saudaveis
  # p: coluna que termina as amostras saudaveis
  
  dataset <- x  #tabela (dataframe)
  col_ini_periodontits <- y    # coluna que inicia as amostras com periodontitis
  col_fim_periodontits <- z   # coluna que termina as amostras com periodontitis
  col_ini_healthy <- w    # coluna que inicia as amostras saudaveis
  col_fim_healthy <- p    # coluna que termina as amostras saudaveis
  
  taxon <- c()
  w_value <- c()
  p_value <- c()
  p_value_corrected <- c()
  periodontitis_mean <- c()
  periodontitis_sd <- c()
  group_higt_median <- c()
  periodontitis_median <- c()
  periodontitis_IQR <- c()
  periodontitis_1st_quantile <- c()
  periodontitis_3st_quantile <- c()
  periodontitis_max_value <- c()
  periodontitis_min_value <- c()
  healthy_mean <- c()
  healthy_sd <- c()
  healthy_median <- c()
  healthy_IQR <- c()
  healthy_1st_quantile <- c()
  healthy_3st_quantile <- c()
  healthy_max_value <- c()
  healthy_min_value <- c()
  
  # loop executado em cada linha do conjunto de dados
  for (i in 1:nrow(dataset)) {
    # intervalo de colunas dos doentes
    periodontits <- as.numeric(paste(c(dataset[i,col_ini_periodontits:col_fim_periodontits]))) 
    
    # intervalo de linhas dos saudaveis
    healthy <- as.numeric(paste(c(dataset[i,col_ini_healthy:col_fim_healthy])))  # intervalo de linhas dos saudaveis
    
    # media de cada grupo
    periodontitis_mean[i] <- mean(periodontits)
    healthy_mean[i] <- mean(healthy)
    
    # desvio padrao (sd) de cada grupo
    periodontitis_sd[i] <- sd(periodontits)
    healthy_sd[i] <- sd(healthy)
    
    # mediana de cada grupo
    periodontitis_median[i] <- median(periodontits)
    healthy_median[i] <- median(healthy)
    
    # 1st e 3st quartil
    periodontitis_1st_quantile[i] <- as.numeric(round(quantile(periodontits, c(0.25)),2))
    periodontitis_3st_quantile[i] <- as.numeric(round(quantile(periodontits, c(0.75)),2))
    healthy_1st_quantile[i] <- as.numeric(round(quantile(healthy, c(0.25)),2))
    healthy_3st_quantile[i] <- as.numeric(round(quantile(healthy, c(0.75)),2))
    
    # Amplitude Inter Quartil do ingles: Interquartile Range (IQR)
    periodontitis_IQR[i] <- IQR(periodontits)
    healthy_IQR[i] <- IQR(healthy)
    
    # verificando qual grupo tem maior mediana
    if (periodontitis_median[i] == healthy_median[i]){
      group_higt_median[i] <- "Healthy=Periodontitis"
    }
    if (periodontitis_median[i] > healthy_median[i]){
      group_higt_median[i] <- "Periodontitis"
    }
    if (periodontitis_median[i] < healthy_median[i]){
      group_higt_median[i] <- "Healthy"
    }

    # valores max e min de cada grupos
    periodontitis_max_value[i] <- max(periodontits)
    periodontitis_min_value[i] <- min(periodontits)
    healthy_max_value[i] <- max(healthy)
    healthy_min_value[i] <- min(healthy)
    
    # executando teste de mann-whitney
    result_test <- wilcox.test(healthy,periodontits)  # teste
    taxon[i] <- paste(dataset[i,1]) # taxon (salvando)
    p_value[i] <- result_test[["p.value"]] # p-value (salvando)
    w_value[i] <- result_test[["statistic"]][["W"]] # w-value (salvando)
    
  }
  
  #correção de p pelo metodo Benjamini & Hochberg
  p_value_corrected <- p.adjust(p_value, method = "fdr")
  
  # tabulando dados
  result_table <- data.frame(taxon,w_value,p_value,p_value_corrected,
                             group_higt_median,periodontitis_mean,periodontitis_sd, periodontitis_median,
                             periodontitis_IQR,periodontitis_1st_quantile,
                             periodontitis_3st_quantile,periodontitis_max_value,
                             periodontitis_min_value,healthy_mean,healthy_sd,healthy_median,healthy_IQR,
                             healthy_1st_quantile,healthy_3st_quantile,healthy_max_value,
                             healthy_min_value)
  return(result_table)
}
```
