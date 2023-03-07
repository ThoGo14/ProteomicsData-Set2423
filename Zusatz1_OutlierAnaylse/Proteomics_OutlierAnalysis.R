library(tidyverse)
library(matrixStats)
library(openxlsx)
library(stats)
library(plotly)
library(manipulateWidget)
library(htmlwidgets)

theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("../ProteomicsNRE2__Expression.RData")
ProteomicsData_expression <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]


## Table ----
MyoID = ProteomicsDataNRE2[["ProteomicsData__Metadata"]]$MyoID

Protein_MeanSD_List <- list()

for (TP in c("pre_resting", "pre_acuteexercise", "post_acuteexercise", "post_resting")) {
  Protein_MeanSD <- data.frame(Proteins = rownames(ProteomicsData_expression))
  
  Proteomics_short <- ProteomicsData_expression[, grep(paste0("__", TP), colnames(ProteomicsData_expression))]
  
  ## 4209 == M102
  Mean <- t(t(rowMeans(Proteomics_short)))
  colnames(Mean) = paste0(TP, "_Mean")
  
  SD <- t(t(rowSds(as.matrix(Proteomics_short))))
  colnames(SD) = paste0(TP, "_SD")
  
  
  Protein_MeanSD <- cbind.data.frame(Protein_MeanSD, Mean, SD, Proteomics_short)
  Protein_MeanSD_List[[paste0("TP_", TP)]] <- Protein_MeanSD
}

write.xlsx(Protein_MeanSD_List, file = "Proteomics_OutlierAnalysis.xlsx")


## PCS Analysis ----
MyoID = sapply(colnames(ProteomicsData_expression), function (x) unlist(str_split(x, "__"))[1])
TP = sapply(colnames(ProteomicsData_expression), function (x) unlist(str_split(x, "__"))[2])
samples = colnames(ProteomicsData_expression)

pca = prcomp(t(ProteomicsData_expression), scale. = FALSE, retx = TRUE)
pca.res_summary <- summary(pca)[["importance"]]
summary(pca)
pca$x

screeplot(pca)
Screeplot <- pca.res_summary %>% as_tibble(rownames = "Datatype") %>% 
  pivot_longer(cols = -1, names_to = "PC", values_to = "values") %>% 
  filter(Datatype == "Proportion of Variance") %>% 
  mutate(PC = factor(PC, levels = paste0("PC", 1:ncol(pca.res_summary)))) %>% 
  filter(PC %in% paste0("PC", 1:10)) %>% 
  ggplot(aes(x = PC,
             y = values)) +
  geom_col() +
  # geom_line(aes(group = 1)) +
  geom_label(aes(x = PC, y = 0.25, label = round(values, 4)*100)) +
  ylim(0, 1) +
  labs(title="Screeplot") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_blank()
  )

pc.var<-pca$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

pca.df <-as_tibble(pca$x)
p1 <- pca.df %>% 
  ggplot(aes(x = PC1,
             y = PC2,
             label = samples,
             color = MyoID,
             )) +
  geom_point(size = 3, aes(shape = TP)) +
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) + 
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  labs(title = "PCA plot") +
  theme(
    legend.position = "right"
  )
ggplotly(p1)  


save(Screeplot, p1, file = "OutlierAnylsis.RData")

combined_plots <- combineWidgets(ggplotly(Screeplot), 
                                 ggplotly(p1), 
                                 ncol = 2,
                                 title = "PCA Analyse von Olink qPCR Daten")
saveWidget(combined_plots, selfcontained = TRUE, file = "OutlierAnalysis_PCA.html")

## Select outliers by SD and mean (>4x)
apply(pca.df, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) ))

## Select outliers by MAD (median absolute deviation) and median (>9x)
apply(pca.df, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) # %>% reduce(union)
row.names(pca$x)[c(4, 76, 83)]


dist <- apply(pca.df, 2, function(x) abs(x - median(x)) / mad(x)) %>%
  apply(1, max)

qplot(y = sort(dist, decreasing = TRUE)) +
  geom_hline(yintercept = 6, color = "red")
