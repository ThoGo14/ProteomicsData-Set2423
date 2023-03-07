library(openxlsx)
library(tidyverse)
theme_set(theme_classic())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")

# ProteomicsData__Expression_Head <- as.matrix(ProteomicsDataNRE2[["ProteomicsData__Expression_Head"]])
# ProteomicsData__Expression_log2 <- log2(ProteomicsDataNRE2[["ProteomicsData__Expression"]])
# 
# ProteomicsData__Expression <- ProteomicsDataNRE2[["ProteomicsData__Expression"]]
# 
Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]
Proteomics_Protein_Info_IDs <- ProteomicsDataNRE2[["Proteomics_Protein_Info_IDs"]]

# row.names(Proteomics_Protein_Info) <- Proteomics_Protein_Info$Gene_name

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/20210315_Klinische_Daten.xlsx")

load("05_Korrelation_ProteomicsNRE2_mit_KlinischenDaten.RData")

#### Selecting all Mitochondria Proteins ----
## Select Mitochondrion (GO:0005739) 
Proteomics_Protein_Info_GO <- ProteomicsDataNRE2[["Proteomics_Protein_Info_GO"]]
Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]

ProteomicsDataNRE2_Mitochondrion_GO <- Proteomics_Protein_Info_GO %>% 
  filter(GO.Domain == "cellular_component", GO.ID == "GO:0005739") %>% 
  select(Gene_name, GO.Term) %>% 
  distinct(Gene_name, .keep_all = TRUE)

ProteomicsDataNRE2_Mitochondrion_GO <- rename(ProteomicsDataNRE2_Mitochondrion_GO,
                                              Protein = Gene_name,
                                              !!paste0("Mitochondrion (GO:0005739) (",
                                                       nrow(ProteomicsDataNRE2_Mitochondrion_GO), ")") := GO.Term)


## Select Mitochondria from MitoCarta

MitoCarta <- read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/MitoCarta/Human.MitoCarta3.0.xlsx",
                       sheet = "A Human MitoCarta3.0")

MitoCarta_short <- MitoCarta %>% 
  select(UniProt, MitoCarta3.0_MitoPathways, MitoCarta3.0_SubMitoLocalization) %>% 
  filter(UniProt %in% Proteomics_Protein_Info$Uniprot_ID) %>% 
  mutate(MitoCarta3.0_MitoPathways = ifelse(MitoCarta3.0_MitoPathways == "0",
                                            NA, MitoCarta3.0_MitoPathways))

MitoCarta_short <- rename(MitoCarta_short, 
                          Uniprot_ID = UniProt,
                          !!paste0("MitoCarta3.0_MitoPathways (", 
                                   table(!is.na(MitoCarta_short$MitoCarta3.0_MitoPathways))["TRUE"], ")") := MitoCarta3.0_MitoPathways,
                          !!paste0("MitoCarta3.0_SubMitoLocalization (", 
                                   nrow(MitoCarta_short), ")") := MitoCarta3.0_SubMitoLocalization)
MitoCarta_short$Protein <- Proteomics_Protein_Info[grep(paste0(MitoCarta_short$Uniprot_ID, collapse = "|"), 
                                                        Proteomics_Protein_Info$Uniprot_ID),
                                              "Gene_name"]


## Combine both Tables

Mito_Proteins <- unique(c(ProteomicsDataNRE2_Mitochondrion_GO$Protein,
                          MitoCarta_short$Protein))


## Korrelationsanalyse ----
Timepoints <- names(KorrelationTable_Results_List)
Welche_KD <- colnames(KorrelationTable_Results_List[[1]])[(grep("Correlation__Analysis", 
                                                               colnames(KorrelationTable_Results_List[[1]]))+1L):
                                                            (grep("Gene Annotations",
                                                                 colnames(KorrelationTable_Results_List[[1]]))-1L)]
Welche_KD <- unique(sapply(Welche_KD, function (x) unlist(str_split(x, "__"))[1]))
  

All_Proteins = Proteomics_Protein_Info$Uniprot_ID
names(All_Proteins) = Proteomics_Protein_Info$Gene_name

for(TP in Timepoints){
  
  ProteinTable <- ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]]
  colnames(ProteinTable) <- sapply(colnames(ProteinTable), function(x){
    x = unlist(str_split(x, "_"))[1]
    return(x)
  })
  ProteinTable <- data.frame("MyoID" = colnames(ProteinTable), t(ProteinTable), check.names = F)
  
  Table_For_Korrelation <- right_join(KlinischeDaten, ProteinTable, by = "MyoID")
  
  for(KD in Welche_KD){
    for(Protein in All_Proteins){
      
      Gene_Name = names(which(Protein == All_Proteins))
      
      Correlation <- cor.test(Table_For_Korrelation[, Protein], 
                              Table_For_Korrelation[, KD])
      
      if (Correlation[["p.value"]] < 0.05) {
        if (Correlation[["p.value"]] > 0.0001) {
          
          Annot_Text <- paste("r = ", format(round(Correlation[["estimate"]][["cor"]], digits = 4), nsmall = 4), 
                              "\np-Value = ", format(round(Correlation[["p.value"]], digits = 4), nsmall = 4), 
                              sep = "")
        } else {
          
          Annot_Text <- paste("r  = ", format(round(Correlation[["estimate"]][["cor"]], digits = 4), nsmall = 4),
                              "\np-Value = ", formatC(Correlation[["p.value"]], format = "e",  digits = 2), 
                              sep = "")
        }
        
        if(Correlation[["estimate"]][["cor"]] > 0) {
          annotations <- data.frame(
            xpos = Inf,
            ypos =  -Inf,
            hjustvar = 1,
            vjustvar = -0.25) #<- adjust
        } else {
          annotations <- data.frame(
            xpos = Inf,
            ypos =  Inf,
            hjustvar = 1,
            vjustvar = 1) #<- adjust
        }
        
        Table_For_Korrelation %>% 
          ggplot(aes(y = .data[[Protein]],
                     x = .data[[KD]])) +
          geom_smooth(method = "lm", formula = "y ~ x", se = T, size = 1.5, color = "red", fill = "red", alpha = .2) +
          geom_point(size = 4, stroke = 2) +
          geom_text(data = annotations,aes(x = xpos,y = ypos, hjust = hjustvar, vjust = vjustvar, label = Annot_Text), size = 7) +
          labs(title = TP,
               y = paste0(Gene_Name, " (", Protein, ")"),
               x = KD) +
          theme(
            axis.text = element_text(size = 22, color = "black"),
            axis.title.x = element_text(size = 22, vjust = 0),
            axis.title.y = element_text(size = 22, vjust = 2),
            legend.position = "top",
            legend.text = element_text(size = 16),
            # legend.title = element_text(size = 15, face = "bold"),
            plot.title = element_text(face = "bold", size = 25),
            axis.ticks.length = unit(.25, "cm"),
          ) 
        
        ggsave(filename = paste("16_Korrelation mit Klinischen Daten/", TP, "/",
                                "Korrelation__", TP, "__", KD, "_", Gene_Name, 
                                ifelse(Gene_Name %in% Mito_Proteins,
                                       "", "__NonMito"),
                                ".png", sep = ""), 
               width = 21, height = 15, dpi = 200, units = "cm")
      }
   
    }
  }
  
  print(paste0(format(Sys.time(), "%x %T"), "  Finished:  ", TP))
}



