mann_whitney_analysis <- function(x,y,z,w,p) {
  
  # x: dataframe
  # y: coluna que inicia as amostras do grupo 1
  # z: coluna que termina as amostras do grupo 1
  # w: coluna que inicia as amostras do grupo 2
  # p: coluna que termina as amostras do grupo 2
  
  # Armazenando informações fornecidas pelo usuário
  dataset <- x
  col_ini_group1 <- y
  col_fim_group1 <- z
  col_ini_group2 <- w
  col_fim_group2 <- p

  # Criando variáveis
  variable <- c()
  w_value <- c()
  p_value <- c()
  p_value_corrected_FDR <- c()
  p_value_corrected_Bonferroni <- c()
  p_value_corrected_Holm <- c()
  p_value_corrected_Hommel <- c()
  p_value_corrected_Hochberg <- c()
  group2_mean <- c()
  group2_sd <- c()
  group2_median <- c()
  group2_IQR <- c()
  group2_1st_quantile <- c()
  group2_3st_quantile <- c()
  group2_max_value <- c()
  group2_min_value <- c()
  group1_mean <- c()
  group1_sd <- c()
  group1_median <- c()
  group1_IQR <- c()
  group1_1st_quantile <- c()
  group1_3st_quantile <- c()
  group1_max_value <- c()
  group1_min_value <- c()
  group_higt_mean <- c()
  group_higt_median <- c()

  
  # Loop executado para cada linha do dataset
  for (i in 1:nrow(dataset)) {
    
    # Capturando dados do grupo 1
    group1 <- as.numeric(paste(c(dataset[i,col_ini_group1:col_fim_group1])))  
    
    # Capturando dados do grupo 2
    group2 <- as.numeric(paste(c(dataset[i,col_ini_group2:col_fim_group2]))) 
    
    # Calculando média
    group1_mean[i] <- mean(group1)
    group2_mean[i] <- mean(group2)
    
    # Calculando desvio-padrão
    group1_sd[i] <- sd(group1)
    group2_sd[i] <- sd(group2)
    
    # Calculando mediana
    group1_median[i] <- median(group1)
    group2_median[i] <- median(group2)
    
    # Calculando 1st e 3st quartil
    group1_1st_quantile[i] <- as.numeric(round(quantile(group1, c(0.25)),2))
    group1_3st_quantile[i] <- as.numeric(round(quantile(group1, c(0.75)),2))
    group2_1st_quantile[i] <- as.numeric(round(quantile(group2, c(0.25)),2))
    group2_3st_quantile[i] <- as.numeric(round(quantile(group2, c(0.75)),2))
    
    # Calculando Interquartile Range (IQR)
    group1_IQR[i] <- IQR(group1)
    group2_IQR[i] <- IQR(group2)
    
    # Calculando qual grupo tem maior média
    if (group2_mean[i] == group1_mean[i]){
      group_higt_mean[i] <- "group1=group2"
    }
    if (group2_mean[i] < group1_mean[i]){
      group_higt_mean[i] <- "group1"
    }
    if (group2_mean[i] > group1_mean[i]){
      group_higt_mean[i] <- "group2"
    }  
    
    # Calculando qual grupo tem maior mediana
    if (group2_median[i] == group1_median[i]){
      group_higt_median[i] <- "group1=group2"
    }
    if (group2_median[i] < group1_median[i]){
      group_higt_median[i] <- "group1"
    }
    if (group2_median[i] > group1_median[i]){
      group_higt_median[i] <- "group2"
    }

    # Valores max e min de cada grupos
    group1_max_value[i] <- max(group1)
    group1_min_value[i] <- min(group1)
    group2_max_value[i] <- max(group2)
    group2_min_value[i] <- min(group2)

    # Executando teste de mann-whitney
    result_test <- wilcox.test(group1,group2)
    
    # Salvando nome da variável
    variable[i] <- paste(dataset[i,1])
    
    # Salvando p-value 
    p_value[i] <- result_test[["p.value"]]
    
    # Salvando w-value
    w_value[i] <- result_test[["statistic"]][["W"]]
    
  } # fim do loop
  
  
  # Corrigindo p-value (Benjamini & Hochberg)
  p_value_corrected_FDR <- p.adjust(p_value, method = "fdr")
  
  # Corrigindo p-value (Bonferroni)
  p_value_corrected_Bonferroni <- p.adjust(p_value, method = "bonferroni")
  
  # Corrigindo p-value (Holm)
  p_value_corrected_Holm <- p.adjust(p_value, method = "holm")
  
  # Corrigindo p-value (Hommel)
  p_value_corrected_Hommel <- p.adjust(p_value, method = "hommel")
  
  # Corrigindo p-value (Hochberg)
  p_value_corrected_Hochberg <- p.adjust(p_value, method = "hochberg")

  # Gerando dataframe com resultados
  result_table <- data.frame(variable,w_value,p_value, p_value_corrected_FDR,
                             p_value_corrected_Bonferroni, p_value_corrected_Holm,
                             p_value_corrected_Hommel, p_value_corrected_Hochberg,
                             group1_mean,group1_sd,group1_median,group1_IQR,
                             group1_1st_quantile,group1_3st_quantile,group1_max_value,
                             group1_min_value,
                             group2_mean,group2_sd, group2_median,
                             group2_IQR,group2_1st_quantile,
                             group2_3st_quantile,group2_max_value,
                             group2_min_value, group_higt_mean, group_higt_median)
  return(result_table)
}

library(readxl)
tabela <- read_excel("input.xlsx")
resultado <- mann_whitney_analysis(tabela,2,6,7,11)
write.table(resultado, file = "output.txt", sep = "\t", row.names = FALSE)













