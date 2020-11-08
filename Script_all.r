#############################################################################
#
#
# Todos os pacotes necessários
#
#
#############################################################################

library(dplyr)
library(tidyverse)
library(stringr)
library(tibble)
library(gplots)
library(ggplot2)
library(gridExtra)
library(lattice)
library(psych)
library(multtest)
library(permute)
library(vegan)
library(pheatmap)
library(ggConvexHull) # font: https://github.com/cmartin/ggConvexHull

#############################################################################
#
# Formatando arquivo de saída da Análise pelo Metaphlan 2.0
#    input = AbundanceTable_to_R.csv
#    output = GeneralDataset.csv
#
#  The "GeneralDataset.csv" file will be used in the following analyzes
#
#############################################################################

# Carregando dados
dataset <- read.csv('AbundanceTable_to_R.csv', header = TRUE, sep = ',', dec = ".")

# Declarando variaveis
col_ini_healthy <- 2          # coluna que inicia as amostras saudaveis
col_fim_healthy <- 14         # coluna que termina as amostras saudaveis
col_ini_periodontits <- 15    # coluna que inicia as amostras com periodontitis
col_fim_periodontits <- 28    # coluna que termina as amostras com periodontitis
id <- 1:nrow(dataset)
taxon <- c()
level <- c()
w_value <- c()
p_value <- c()
p_value_corrected <- c()
periodontitis_mean <- c()
periodontitis_sd <- c()
group_high_median <- c()
group_high_mean <- c()
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
abundances_values <- dataset[col_ini_healthy:col_fim_periodontits]

# Loop executado em cada linha 
for (i in 1:nrow(dataset)) {
  # Restringindo intervalo de dados (Periodontitis)
  periodontits <- as.numeric(paste(c(dataset[i,col_ini_periodontits:col_fim_periodontits])))
  
  # Restringindo intervalo de dados (Healthy)
  healthy <- as.numeric(paste(c(dataset[i,col_ini_healthy:col_fim_healthy])))
  
  # Calculando Média
  periodontitis_mean[i] <- mean(periodontits)
  healthy_mean[i] <- mean(healthy)
  
  # Calculando Desvio-Padrão
  periodontitis_sd[i] <- sd(periodontits)
  healthy_sd[i] <- sd(healthy)
  
  # Calculando Mediana
  periodontitis_median[i] <- median(periodontits)
  healthy_median[i] <- median(healthy)
  
  # Calculando 1st e 3st quartil
  periodontitis_1st_quantile[i] <- as.numeric(round(quantile(periodontits, c(0.25)),2))
  periodontitis_3st_quantile[i] <- as.numeric(round(quantile(periodontits, c(0.75)),2))
  healthy_1st_quantile[i] <- as.numeric(round(quantile(healthy, c(0.25)),2))
  healthy_3st_quantile[i] <- as.numeric(round(quantile(healthy, c(0.75)),2))
  
  # Calculando Interquartile Range (IQR)
  periodontitis_IQR[i] <- IQR(periodontits)
  healthy_IQR[i] <- IQR(healthy)
  
  # Comparando qual grupo apresenta maior Mediana
  if (periodontitis_median[i] == healthy_median[i]){
    group_high_median[i] <- "Healthy=Periodontitis"
  }
  if (periodontitis_median[i] > healthy_median[i]){
    group_high_median[i] <- "Periodontitis"
  }
  if (periodontitis_median[i] < healthy_median[i]){
    group_high_median[i] <- "Healthy"
  }
  
  # Comparando qual grupo apresenta maior Media
  if (periodontitis_mean[i] == healthy_mean[i]){
    group_high_mean[i] <- "Healthy=Periodontitis"
  }
  if (periodontitis_mean[i] > healthy_mean[i]){
    group_high_mean[i] <- "Periodontitis"
  }
  if (periodontitis_mean[i] < healthy_mean[i]){
    group_high_mean[i] <- "Healthy"
  }
  
  # Obtendo maior valor
  periodontitis_max_value[i] <- max(periodontits)
  periodontitis_min_value[i] <- min(periodontits)
  
  # Obtendo menor valor
  healthy_max_value[i] <- max(healthy)
  healthy_min_value[i] <- min(healthy)
  
  # Salvando nivel taxonomico para criar uma coluna especifica
  simbol_tax = str_sub(paste(dataset[i,1]), start = 1, end = 3)
  if (str_detect("k__",simbol_tax) == TRUE){level[i] <- "domain"} else {
    if (str_detect("p__",simbol_tax) == TRUE){level[i] <- "phylum"} else {
      if (str_detect("c__",simbol_tax) == TRUE){level[i] <- "class"} else{
        if (str_detect("o__",simbol_tax) == TRUE){level[i] <- "order"} else{
          if (str_detect("f__",simbol_tax) == TRUE){level[i] <- "family"} else{
            if (str_detect("g__",simbol_tax) == TRUE){level[i] <- "genus"} else{
              if (str_detect("s__",simbol_tax) == TRUE){level[i] <- "specie"} else{
                level[i] <- "NaN"
              }
            }
          }
        }
      }
    }
  }
  
  # Removendo "_" e outros caracteres dos Taxons
  taxon[i] <- str_replace_all(str_sub(paste(dataset[i,1]), start = 4), "_", " ")
  
  # Aplicando o Mann-Whitney Test
  result_test <- wilcox.test(healthy,periodontits)
  p_value[i] <- result_test[["p.value"]]
  w_value[i] <- result_test[["statistic"]][["W"]]
}
remove(i)

# Corrigindo p-value pelo Método Benjamini & Hochberg
p_corrected_FDR <- p.adjust(p_value, method = "fdr")

# Tabulando Resultados
result_table <- data.frame(id, level, taxon,w_value, p_value, p_corrected_FDR,
                           group_high_median, group_high_mean, periodontitis_mean, periodontitis_sd,
                           periodontitis_median,
                           periodontitis_IQR, periodontitis_1st_quantile,
                           periodontitis_3st_quantile, periodontitis_max_value,
                           periodontitis_min_value, healthy_mean,healthy_sd, healthy_median, healthy_IQR,
                           healthy_1st_quantile, healthy_3st_quantile, healthy_max_value,
                           healthy_min_value, abundances_values)

# Salvando resultados tabulados
write.csv(result_table, 
          file = "GeneralDataset.csv",
          sep = ",",
          quote = FALSE,
          row.names = FALSE, 
          col.names = NA)



#############################################################################
#
#
# Shannon Diversity Analysis
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando dataset
dataset <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Especie" para analise
dados_especies <- filter(dataset, level == "specie")

# Removendo 'unclassified' e 'noname' dos dados
dados_especies <- dados_especies %>% filter(!str_detect(taxon, 'unclassified'))
dados_especies <- dados_especies %>% filter(!str_detect(taxon, 'noname'))

# Preparando tabela (Healthy group)
healthy_table <- select(dados_especies, H01, H02, H03, H04, H05, H06, H07, H08, H09, H10, H11, H12, H13)
row.names(healthy_table) <- dados_especies$taxon
healthy_table <- t(healthy_table)

# Preparando tabela (Periodontitis group)
periodontitis_table <- select(dados_especies, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14)
row.names(periodontitis_table) <- dados_especies$taxon
periodontitis_table <- t(periodontitis_table)

# Calculando diversidade de Shannon
# Periodontitis group
shannon_periodontitis <- diversity(periodontitis_table, index = "shannon")
shannon_periodontitis <- as.numeric(shannon_periodontitis)
# Healthy group
shannon_healthy <- diversity(healthy_table, index = "shannon")
shannon_healthy <- as.numeric(shannon_healthy)


# Aplicando Shapiro test sobre os dados de Diversidade de Shannon
# Periodontitis group
shapiro_periodontitis <- shapiro.test(shannon_periodontitis)
# Healthy group
shapiro_healthy <- shapiro.test(shannon_healthy)

# Mann-Whitney test
p_value <- wilcox.test(shannon_healthy,shannon_periodontitis)$p.value

# Definindo intervalo do p-value
if (p_value > 0.05) {
  p_simbol <- "p>0.05"
} else {
  if (p_value < 0.001) {
    p_simbol <- "p<0.001"
  } else{
    if (p_value < 0.01) {
      p_simbol <- "p<0.01"
    } else {
      if (p_value < 0.05) {
        p_simbol <- "p<0.05"
      }
    }
  }
}

# Criando Dataframe para boxplot
groups <- c(rep("Healthy",length(shannon_healthy)),rep("Periodontitis",length(shannon_periodontitis)))
samples <- c(row.names(healthy_table),row.names(periodontitis_table))
shannon_diversity <- c(shannon_healthy,shannon_periodontitis)
shannon_index_table <- data.frame(samples,shannon_diversity)

# Criando boxplot sem jitter
boxplot_ShannonDiversity <- ggplot(shannon_index_table, aes(x=groups, y=shannon_diversity, fill=groups)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=3) +
  annotate("text", y=4.1, x=0.5, label=paste("Mann Whitney,",p_simbol), size = 3, hjust = 0) +
  scale_fill_manual(values=c("#008000", "#ff0000")) +
  labs(title = "",
       x = "",
       y = "Shannon's Diversity Index",
       color = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = rel(1.2), angle = 90, color = "black"),
        axis.text.y = element_text(size = rel(1.2),color = "black"),
        axis.text.x = element_text(size = rel(1.2),color = "black"),
        legend.position = "none")
boxplot_ShannonDiversity

# Criando semente para o parâmetro jitter do gráfico
set.seed(123)

# Boxplot com jitter
boxplot_ShannonDiversity_Jitter <- ggplot(shannon_index_table, aes(x=groups, y=shannon_diversity, fill=groups)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2)) +
  annotate("text", y=4.1, x=0.5, label="Mann Whitney, p<0.001", size = 3, hjust = 0) +
  scale_fill_manual(values=c("#008000", "#ff0000")) +
  labs(title = "",
       x = "",
       y = "Shannon's Diversity Index",
       color = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = rel(1.2), angle = 90, color = "black"),
        axis.text.y = element_text(size = rel(1.2),color = "black"),
        axis.text.x = element_text(size = rel(1.2),color = "black"),
        legend.position = "none")
boxplot_ShannonDiversity_Jitter

# Criando diretório para guardar os resultados
dir.create("ShannonDiversity_result")

# Salvando gráfico "boxplot_ShannonDiversity.png"
png(filename = "ShannonDiversity_result/boxplot_ShannonDiversity.png", width = 10, height = 10, units = "cm", pointsize = 25, res=300, type = c("cairo"))
boxplot_ShannonDiversity
dev.off()

# Salvando gráfico "boxplot_ShannonDiversity_Jitter.png"
png(filename = "ShannonDiversity_result/boxplot_ShannonDiversity_Jitter.png", width = 10, height = 10, units = "cm", pointsize = 25, res=300, type = c("cairo"))
boxplot_ShannonDiversity_Jitter
dev.off()

# Salvando tabela com indíces de shannon
write.csv(shannon_index_table, file = "ShannonDiversity_result/shannon_index_table.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = NA)



#############################################################################
#
#
# Non-metric multidimensional scaling (NMDS)
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando dataset
tabela_geral <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Especie" para analise
tabela_filtrada <- filter(tabela_geral, level == "specie")

# Removendo 'unclassified' e 'noname' dos dados
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'unclassified'))
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'noname'))

# Preparando tabela e convertendo em matriz
dados_especies <- select(tabela_filtrada, H01, H02, H03, H04, H05, H06, H07, H08, H09, H10, H11, H12, H13, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14)
dados_especies <- t(dados_especies)
colnames(dados_especies) <- tabela_filtrada$taxon
dados_especies <- as.matrix(dados_especies)

# Criando vetor com grupos para análise ANOSIM e construção do gráfico NMDS
groups <- c(rep("Healthy", 13), rep("Periodontitis", 14))

# Análise NMDS
nmds <- metaMDS(dados_especies, distance = "bray")
nmds_estress <- nmds$stress
plot(nmds)

# Vetor scores do NMDS
NMDS_scores_table <- as.data.frame(scores(nmds)) 

# Análise ANOSIM
anosim_result <- anosim(dados_especies, groups, distance = "bray", permutations = 9999, parallel = 4)
# Valor de R da ANOSIM
anosim_R <- paste("R=", round(anosim_result$statistic,3))
# Definindo intervalo de R da ANOSIM
if (anosim_result$signif < 0.001) {anosim_p <- "p<0.001"
} else { 
  if (anosim_result$signif < 0.01) {
    anosim_p <- "p<0.01"
  } else {
    if (anosim_result$signif < 0.05) {
      anosim_p <- "p<0.05"
    } else {
      if (anosim_result$signif > 0.05) {
        anosim_p <- "p>0.05"
      }
    }
  }
}

# Vetores necessários para criar tabela utilizada para construção do gráfico NMDS
samples <- rownames(dados_especies)
groups <- groups
NMDS1 <- NMDS_scores_table$NMDS1
NMDS2 <- NMDS_scores_table$NMDS2

# Criando tabela que utlizada para criação do gráfico NMDS
NMDS_table <- data.frame(samples, groups, NMDS1, NMDS2)
colnames(NMDS_table) <- c("Sample","Groups","NMDS1","NMDS2")

# Gráfico NMDS
NMDS_graph <- ggplot(NMDS_table, aes(x = NMDS1, y = NMDS2, color = Groups)) +
  geom_point() +
  geom_convexhull(alpha = 0.3,aes(fill = Groups)) +
  geom_text(aes(label=Sample),hjust=-0.1, vjust=-0.3, show.legend = FALSE) +
  annotate("text", y=1, x=-0.7, label=paste("Stress: ",round(nmds_estress,3),"\n","ANOSIM: ",anosim_p,", ",anosim_R, sep = ""), size = 3, hjust = 0) +
  scale_fill_manual(values = c("#008000", "#ff0000")) +
  scale_colour_manual(values = c("#008000", "#ff0000")) +
  ylim(-0.7,1) +
  xlim(-0.7,1) +
  labs(x = 'NMDS1',
       y = 'NMDS2') + 
  theme_minimal() +
  theme(axis.title.y = element_text(size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),
        legend.title = element_text(color = "white"),
        legend.text = element_text(size = rel(1.2)),
        legend.position = "top")
NMDS_graph

# Criando diretório para guardar os resultados
dir.create("NMDS_result")

# Salvando gráfico "NMDS_species.png"
png(filename = "NMDS_result/NMDS_species.png", width = 15, height = 10, units = "cm", res=300, type = c("cairo"))
NMDS_graph
dev.off()

# Salvando tabela com dados do NMDS
write.csv(NMDS_table, file = "NMDS_result/NMDS_table.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = NA)



#############################################################################
#
#
# Heatmap Analysis
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando tabelas
table <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Especie" para analise
table <- filter(table, level == "specie")

# Removendo 'unclassified' e 'noname' dos dados
table <- table %>% filter(!str_detect(taxon, 'unclassified'))
table <- table %>% filter(!str_detect(taxon, 'noname'))

# Definindo nível de significância de p<0.05 sobre os valores de p corrigidos por FDR
table <- filter(table,as.numeric(p_corrected_FDR) < 0.05)

# Armazenando nomes dos táxons em vetor
taxons <-  table$taxon

# Selecionando amostras e renomeando as linhas com os táxons armazerenados no vetor 'taxons'
table <- select(table, H01, H02, H03, H04, H05, H06, H07, H08, H09, H10, H11, H12, H13, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14)
rownames(table) <- taxons

# Mudando nome da tabela e removendo variaveis que não são mais necessárias
abundance_table <- table
remove(table, taxons)

# Coletando nomes das espécies e amostras
species_list <- rownames(abundance_table)
sample_list <-  colnames(abundance_table)

# Criando metadados das amostras com base na presença de 'H' ou 'P' no nome
# Definindo Grupos
Condition <- c()
for (i in 1:length(sample_list)) {
  if (str_count(sample_list[i], "H") > 0) {
    Condition[i] <- "Healthy"
  } else {
    Condition[i] <- "Periodontitis"
  }
}; remove(i)

# Criando Tabela de Metadados para ser usada durante a contrução do heatmap
metadados <- data.frame(Condition)
rownames(metadados) <- sample_list

# Criando lista de anotacoes com base nos metadados
my_colour = list(
  Condition = c(Periodontitis = "#ff0000", Healthy = "#008000")
)

# Tabela de distancia de bray-curtis (Linhas / Species)
dist_bray_row <- vegdist(abundance_table, "bray", diag = TRUE)

# Tabela de distancia de bray-curtis (Colunas / Samples)
dist_bray_col <- vegdist(t(abundance_table), "bray", diag = TRUE)

# Criando Heatmap padronizado por amostra em Z-score
heatmap_Zscore <- pheatmap(abundance_table,
                           color = colorRampPalette(c("navy", "white", "firebrick3"))(1000),
                           scale = "column",
                           clustering_distance_rows = dist_bray_row, # Species
                           clustering_distance_cols = dist_bray_col, # Samples
                           clustering_method="ward.D2",
                           border_color=FALSE,
                           annotation_colors = my_colour,
                           annotation_col=metadados,
                           drop_levels = FALSE,
                           angle_col = "90")

# Criando diretório para guardar os resultados
dir.create("Heatmap_result")

# Salvando gráfico "heatmap_Z-score_species.png"
png(filename = "Heatmap_result/heatmap_Z-score_species.png", width = 25, height = 20, units = "cm", res=300, type = c("cairo"))
heatmap_Zscore
dev.off()



#############################################################################
#
#
# Stacked Bar chart Analysis
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando dataset
tabela_geral <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Phylum" para analise
tabela_filtrada <- filter(tabela_geral, level == "phylum") # filtrando nivel taxonomico especificado

# Removendo 'unclassified' e 'noname' dos dados
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'unclassified'))
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'noname'))

# Preparando tabela
dados_phylum <- select(tabela_filtrada, H01, H02, H03, H04, H05, H06, H07, H08, H09, H10, H11, H12, H13, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14)
dados_phylum <- add_column(dados_phylum, tabela_filtrada$taxon, .before = "H01")
names(dados_phylum)[1] <- "phylum"

# Ajustando para o formato apropriado ao ggplot2
dados_phylum <- pivot_longer(dados_phylum, cols = c(
  "H01","H02","H03","H04","H05","H06",
  "H07","H08","H09","H10","H11","H12",
  "H13","P01","P02","P03","P04","P05",
  "P06","P07","P08","P09","P10","P11",
  "P12","P13","P14"),
  names_to="amostras",
  values_to="abundancia")

# Stacked bar chart 
Stacked_Bar_Chart <- ggplot(data=dados_phylum, aes(x=amostras, y=abundancia, fill=phylum)) +
  geom_bar(stat="identity", color="black", size=0.3) +
  labs(y = 'Relative Abundance (%)',
       x = "",
       fill = "Phylum:") +
  theme_minimal() + 
  theme() +
  theme(axis.title.y = element_text(size = rel(1.2), angle = 90, color = "black"),
        axis.title.x = element_text(size = rel(1.2), color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, vjust=0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  scale_fill_brewer(palette="Set2")
Stacked_Bar_Chart

# Criando diretório para guardar os resultados
dir.create("StackedBar_result")

# Salvando gráfico "StackedBar_phylum.png"
png(filename = "StackedBar_result/StackedBar_phylum.png", width = 20, height = 10, units = "cm", res=300, type = c("cairo"))
Stacked_Bar_Chart
dev.off()



#############################################################################
#
#
# Scatterplot for Taxonomic Level = 'Specie'
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando tabela_filtrada
dataset <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Especie" para analise
tabela_filtrada <- filter(dataset, level == "specie")

# Removendo 'unclassified' e 'noname' dos dados
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'noname'))
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'unclassified'))

# Definindo nível de significância de p<0.05 sobre os valores de p corrigidos por FDR
tabela_filtrada <- filter(tabela_filtrada,as.numeric(p_corrected_FDR) < 0.05)
tabela_filtrada <- select(tabela_filtrada, taxon, periodontitis_mean, healthy_mean, periodontitis_sd, healthy_sd, group_high_mean)

# Calculando R2 (Coeficiente de determinacao)
r2 <- cor(tabela_filtrada$periodontitis_mean, tabela_filtrada$healthy_mean)^2
r2 <- round(r2,3)
r2_label <- paste("R^2 == ", r2)

# Scatterplot nivel = "Specie"
scatterplot_specie <- ggplot(data = tabela_filtrada) +
  aes(x = healthy_mean, y = periodontitis_mean, colour = group_high_mean) +
  geom_point(size = 3) +
  geom_text(aes(label=ifelse(
    (healthy_mean > 2.0), taxon, ifelse(periodontitis_mean > 2.8, taxon, ""))), hjust=0.13, vjust=-0.5, angle = 7, show.legend = FALSE, size = 5) +
  geom_abline(intercept = c(0), slope = 1, color = c("black"), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 15, by=2.5), limits=c(-1, 12)) +
  scale_x_continuous(breaks = seq(-10, 30, by=5), limits=c(-1, 30)) +
  labs(title = "",
       x = "Healthy (%)",
       y = "Periodontitis (%)",
       color = "") +
  annotate("text", y=11.75, x=0.5, label= r2_label, parse =TRUE, size = 5) +
  scale_color_manual(values=c("#008000","#ff0000")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, color = "black"),
        axis.title.x = element_text(size = rel(1.7), color = "black"),
        axis.text.x = element_text(size = rel(1.6),color = "black"),
        axis.text.y = element_text(size = rel(1.6),color = "black"),
        legend.text = element_text(size = rel(1.7)),
        legend.position = "top")
scatterplot_specie

# Criando diretório para guardar os resultados
dir.create("Scatterplot_result")

# Salvando gráfico "scatterplot_specie.png"
png(filename = "Scatterplot_result/scatterplot_specie.png", width = 20, height = 15, units = "cm", pointsize = 25, res=300, type = c("cairo"))
scatterplot_specie
dev.off()



#############################################################################
#
#
# Scatterplot for Taxonomic Level = 'Genus'
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando tabela_filtrada
dataset <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "Genus" para analise
tabela_filtrada <- filter(dataset, level == "genus")

# Removendo 'unclassified' e 'noname' dos dados
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'noname'))
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'unclassified'))

# Definindo nível de significância de p<0.05 sobre os valores de p corrigidos por FDR
tabela_filtrada <- filter(tabela_filtrada,as.numeric(p_corrected_FDR) < 0.05)
tabela_filtrada <- select(tabela_filtrada, taxon, periodontitis_mean, healthy_mean, periodontitis_sd, healthy_sd, group_high_mean)

# Calculando R2 (Coeficiente de determinacao)
r2 <- cor(tabela_filtrada$periodontitis_mean, tabela_filtrada$healthy_mean)^2
r2 <- round(r2,3)
r2_label <- paste("R^2 == ", r2)

# Scatterplot nivel = "Genus"
scatterplot_genus <- ggplot(data = tabela_filtrada) +
  aes(x = healthy_mean, y = periodontitis_mean, colour = group_high_mean) +
  geom_point(size = 3) +
  geom_text(aes(label=ifelse(
    (healthy_mean > 2.0), taxon, ifelse(periodontitis_mean > 2.6, taxon, ""))), hjust=0.13, vjust=-0.5, show.legend = FALSE, size = 5) +
  geom_abline(intercept = c(0), slope = 1, color = c("black"), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 14, by=2), limits=c(-1, 14)) +
  scale_x_continuous(breaks = seq(-10, 30, by=5), limits=c(-1, 30)) +
  labs(title = "",
       x = "Healthy (%)",
       y = "Periodontitis (%)",
       color = "") +
  annotate("text", y=14, x=0.5, label= r2_label, parse =TRUE, size = 5) +
  scale_color_manual(values=c("#008000","#ff0000")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, color = "black"),
        axis.title.x = element_text(size = rel(1.7), color = "black"),
        axis.text.x = element_text(size = rel(1.6),color = "black"),
        axis.text.y = element_text(size = rel(1.6),color = "black"),
        legend.text = element_text(size = rel(1.7)),
        legend.position = "top")
scatterplot_genus

# Criando diretório para guardar os resultados
dir.create("Scatterplot_result")

# Salvando gráfico "scatterplot_specie.png"
png(filename = "Scatterplot_result/scatterplot_genus.png", width = 20, height = 15, units = "cm", pointsize = 25, res=300, type = c("cairo"))
scatterplot_genus
dev.off()



#############################################################################
#
#
# Scatterplot for Taxonomic Level = 'Phylum'
#
#
#############################################################################

# Limpando vetores e dataframes do ambiente
rm(list=ls())

# Carregando tabela_filtrada
dataset <- read.csv('GeneralDataset.csv', header = TRUE, sep = ',', dec = ".")

# Filtrando nivel taxonomico "phylum" para analise
tabela_filtrada <- filter(dataset, level == "phylum")

# Removendo 'unclassified' e 'noname' dos dados
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'noname'))
tabela_filtrada <- tabela_filtrada %>% filter(!str_detect(taxon, 'unclassified'))

# Definindo nível de significância de p<0.05 sobre os valores de p corrigidos por FDR
tabela_filtrada <- filter(tabela_filtrada,as.numeric(p_corrected_FDR) < 0.05)
tabela_filtrada <- select(tabela_filtrada, taxon, periodontitis_mean, healthy_mean, periodontitis_sd, healthy_sd, group_high_mean)

# Calculando R2 (Coeficiente de determinacao)
r2 <- cor(tabela_filtrada$periodontitis_mean, tabela_filtrada$healthy_mean)^2
r2 <- round(r2,3)
r2_label <- paste("R^2 == ", r2)

# Scatterplot nivel = "phylum"
scatterplot_phylum <- ggplot(data = tabela_filtrada) +
  aes(x = healthy_mean, y = periodontitis_mean, colour = group_high_mean) +
  geom_point(size = 3) +
  geom_text(aes(label=ifelse(
    (healthy_mean > 2.0), taxon, ifelse(periodontitis_mean > 0, taxon, ""))), hjust=-0.1, vjust=-0.3, show.legend = FALSE, size = 4) +
  geom_abline(intercept = c(0), slope = 1, color = c("black"), linetype = "dashed") +
  scale_y_continuous(breaks = seq(-10, 30, by=5), limits=c(-1, 30)) +
  scale_x_continuous(breaks = seq(-10, 35, by=5), limits=c(-1, 35)) +
  labs(title = "",
       x = "Healthy (%)",
       y = "Periodontitis (%)",
       color = "") +
  annotate("text", y=30, x=0.5, label= r2_label, parse =TRUE, size = 5) +
  scale_color_manual(values=c("#008000","#ff0000")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, color = "black"),
        axis.title.x = element_text(size = rel(1.7), color = "black"),
        axis.text.x = element_text(size = rel(1.6),color = "black"),
        axis.text.y = element_text(size = rel(1.6),color = "black"),
        legend.text = element_text(size = rel(1.7)),
        legend.position = "top")
scatterplot_phylum

# Criando diretório para guardar os resultados
dir.create("Scatterplot_result")

# Salvando gráfico "scatterplot_specie.png"
png(filename = "Scatterplot_result/scatterplot_phylum.png", width = 20, height = 15, units = "cm", pointsize = 25, res=300, type = c("cairo"))
scatterplot_phylum
dev.off()