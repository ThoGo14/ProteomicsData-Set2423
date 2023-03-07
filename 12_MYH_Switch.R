library(openxlsx)
library(tidyverse)
library(ggsignif)
library(gridExtra)
theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")
ProteomicsData__Expression_log2 <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]
Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]


ProteomicsData__Expression_log2 <- ProteomicsData__Expression_log2[sort(row.names(ProteomicsData__Expression_log2)), ]


MYH_IDs = Proteomics_Protein_Info[grep(paste0("^", c("MYH1", "MYH2", "MYH4", "MYH7"), "$", collapse = "|"),
                                       Proteomics_Protein_Info$Gene_name), "Uniprot_ID"]
names(MYH_IDs) = c("MYH1", "MYH2", "MYH4", "MYH7")


# MYH_IDs = Proteomics_Protein_Info[grep("MYH", Proteomics_Protein_Info$Gene_name), "Uniprot_ID"]
# names(MYH_IDs) = Proteomics_Protein_Info[grep("MYH", Proteomics_Protein_Info$Gene_name), "Gene_name"]

ProteomicsData__Expression_MyH <- ProteomicsData__Expression_log2[MYH_IDs, ]


## Prepare Data as long table for Plots -----
ProteomicsData__Expression_long <- ProteomicsData__Expression_MyH %>% 
  mutate(Gene = names(MYH_IDs)) %>% 
  relocate(Gene, .before = names(.)[1]) %>% 
  gather(., 
         key = "Colname", value = "value",
         2:ncol(.)) %>% 
  separate(Colname, c("MyoID", "TP"), sep = "__", remove = FALSE) %>% 
  mutate(Timepoint = case_when(TP == "pre_resting" ~ "T0",
                               TP == "pre_acuteexercise" ~ "T2",
                               TP == "post_acuteexercise" ~ "T4",
                               TP == "post_resting" ~ "T5",
                               TRUE ~ "UNDFINED"),
         Timepoint = factor(Timepoint, levels = c("T0", "T2", "T4", "T5")))



Mean_Values = ProteomicsData__Expression_long %>% 
  group_by(Timepoint, Gene) %>% 
  summarise(
    Mean = mean(value),
    Sum = sum(value)
  ) %>% 
  arrange(Timepoint)

for(TP in unique(Mean_Values$Timepoint)){
  Sum_all = colSums(Mean_Values[grep(TP, Mean_Values$Timepoint), "Sum"])
  Mean_Values[grep(TP, Mean_Values$Timepoint), "Freq"] <- 
    Mean_Values[grep(TP, Mean_Values$Timepoint), "Sum"] / Sum_all
}

Mean_Values %>% 
  ggplot(aes(x = Timepoint,
             y = Freq,
             fill = Gene,
             label = (round(Freq*100, 2)))) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, vjust = 2),
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 15, face = "bold"),
    legend.background=element_rect(fill = alpha("white", 0)),
    plot.title = element_text(face = "bold", size = 25),
    axis.ticks.length=unit(.25, "cm"),
  ) 



res.aov2 <- aov(value ~ Gene + Timepoint, data = ProteomicsData__Expression_long)
summary(res.aov2)
TukeyHSD(res.aov2)


ProteomicsData__Expression_long$unique = paste(ProteomicsData__Expression_long$Gene,
                                               ProteomicsData__Expression_long$Timepoint,
                                               sep = "_")

res.aov2 <- aov(value ~ unique, data = ProteomicsData__Expression_long)
summary(res.aov2)
TukeyConmparission <- as.data.frame(TukeyHSD(res.aov2)$unique)
TukeyConmparission$x <- sapply(row.names(TukeyConmparission), function (x) unlist(str_split(x, "-"))[1])
TukeyConmparission$y <- sapply(row.names(TukeyConmparission), function (x) unlist(str_split(x, "-"))[2])

len <- length(unique(TukeyConmparission$y))
TukeyConmparissionMatrix <- as.data.frame(matrix(rep(NA, times = len*len), ncol = len))
colnames(TukeyConmparissionMatrix) <- unique(TukeyConmparission$y)
row.names(TukeyConmparissionMatrix) <- unique(TukeyConmparission$x)

for(i in 1:nrow(TukeyConmparission)){
  ForY = TukeyConmparission[i, "y"]
  ForX = TukeyConmparission[i, "x"]
  
  pVal = TukeyConmparission[i, "p adj"]
  
  TukeyConmparissionMatrix[ForX, ForY] <- pVal
}

write.xlsx(TukeyConmparissionMatrix, "12_MYH_Switch_MYH1-2-4-7.xlsx",
           row.names = TRUE, overwrite = TRUE)

# write.xlsx(TukeyConmparissionMatrix, "12_MYH_Switch.xlsx",
#            rowNames = TRUE, overwrite = TRUE)
