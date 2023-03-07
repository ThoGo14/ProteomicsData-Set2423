library(openxlsx)
library(tidyverse)
theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Select Protein, Timepoint and Clinical Data

Welche_KD_forPlot <- c("ISI_FC")
Welche_Prot_forPlot <- c("ANKRD2")

Welche_TP_forPlot <- c("T2T0")
Welche_BLKD_forPlot <- NA



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")

Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/20210315_Klinische_Daten.xlsx")

## Plots with FC Data ----
if(all(!is.na(Welche_TP_forPlot)) & all(!is.na(Welche_KD_forPlot))){  
  for(TP in Welche_TP_forPlot){
    
    ProteinTable <- ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]]
    colnames(ProteinTable) <- sapply(colnames(ProteinTable), function(x){
      x = unlist(str_split(x, "_"))[1]
      return(x)
    })
    ProteinTable <- data.frame("MyoID" = colnames(ProteinTable), t(ProteinTable), check.names = F)
    
    Table_For_Korrelation <- right_join(KlinischeDaten, ProteinTable, by = "MyoID")
    
    for(KD in Welche_KD_forPlot){
      for(Protein in Welche_Prot_forPlot){
        
        Correlation <- cor.test(Table_For_Korrelation[, Protein], 
                                Table_For_Korrelation[, KD])
        
        if (Correlation[["p.value"]] > 0.0001) {
          
          Annot_Text <- paste("r = ", round(Correlation[["estimate"]][["cor"]], digits = 4), 
                              "\np-Value = ", round(Correlation[["p.value"]],  digits = 4), 
                              sep = "")
        } else {
          
          Annot_Text <- paste("r = ", round(Correlation[["estimate"]][["cor"]], digits = 4),
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
        
        
        
        p <- Table_For_Korrelation %>% 
          select(all_of(Protein), all_of(KD)) %>% 
          na.omit %>% 
          ggplot(aes(y = .data[[Protein]], x = .data[[KD]])) +
          geom_smooth(method = "lm", formula = "y ~ x", 
                      se = T, size = 1.5, color = "red", fill = "red", alpha = .2) +
          geom_point(size = 4) +
          geom_text(data = annotations,aes(x = xpos,y = ypos, hjust = hjustvar, vjust = vjustvar, label = Annot_Text), size = 7) +
          # scale_color_manual("", breaks = c("resp", "nonr"), 
          #                    values = c("#0072B2", "#D55E00")) +
          ggtitle(TP) +
          theme(
            axis.text = element_text(size = 22, color = "black"),
            axis.title.x = element_text(size = 22, vjust = 0),
            axis.title.y = element_text(size = 22, vjust = 2),
            legend.position = "bottom",
            legend.text = element_text(size = 20),
            plot.title = element_text(face = "bold", size = 25),
            axis.ticks.length = unit(.25, "cm"),
          ) 
        
        ggsave(plot = p, filename = paste0("08_SelectedPictures Korrelation ProteomicsNRE2 mit KlinischenDaten/", TP, "__",
                                           KD, "__", Protein, ".png"), 
               width = 21, height = 15, dpi = 200, units = "cm")
        
        print(paste0(TP, " -- ", KD, " -- ", Protein))
        
        
      }
      
    }
  }
}

## Plots with Baseline Data ----
TP_Table <- data.frame(TP_old = c("T0", "T2", "T4", "T5"),
                       TP_new = c("pre_resting", "pre_acuteexercise",
                                  "post_acuteexercise", "post_resting"))


if(all(!is.na(Welche_BLKD_forPlot))){
  for(KD in Welche_BLKD_forPlot){
    
    TP_old = unlist(str_split(KD, "__"))[2]
    TP = TP_Table[grep(TP_old, TP_Table$TP_old), "TP_new"]
    
    TP_Cols <- grep(TP, colnames(ProteomicsDataNRE2[["ProteomicsData__Expression_log2"]]))
    
    ProteinTable <- ProteomicsDataNRE2[["ProteomicsData__Expression_log2"]][, TP_Cols]
    colnames(ProteinTable) <- sapply(colnames(ProteinTable), function(x){
      x = unlist(str_split(x, "_"))[1]
      return(x)
    })
    ProteinTable <- data.frame("MyoID" = colnames(ProteinTable), t(ProteinTable), check.names = F)
    
    Table_For_Korrelation <- right_join(KlinischeDaten, ProteinTable, by = "MyoID")
    
    KD_new <- unlist(str_split(KD, "__"))[1]
    
    for(Protein in Welche_Prot_forPlot){
      
      Correlation <- cor.test(Table_For_Korrelation[, Protein], 
                              Table_For_Korrelation[, KD_new])
      
      if (Correlation[["p.value"]] > 0.0001) {
        
        Annot_Text <- paste("r = ", round(Correlation[["estimate"]][["cor"]], digits = 4), 
                            "\np-Value = ", round(Correlation[["p.value"]],  digits = 4), 
                            sep = "")
      } else {
        
        Annot_Text <- paste("r = ", round(Correlation[["estimate"]][["cor"]], digits = 4),
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
      
      
      
      p <- Table_For_Korrelation %>% 
        select(all_of(Protein), all_of(KD_new)) %>% 
        na.omit %>% 
        ggplot(aes(y = .data[[Protein]], x = .data[[KD_new]])) +
        geom_smooth(method = "lm", formula = "y ~ x", 
                    se = T, size = 1.5, color = "red", fill = "red", alpha = .2) +
        geom_point(size = 4) +
        geom_text(data = annotations,aes(x = xpos,y = ypos, hjust = hjustvar, vjust = vjustvar, label = Annot_Text), size = 7) +
        # scale_color_manual("", breaks = c("resp", "nonr"), 
        #                    values = c("#0072B2", "#D55E00")) +
        ggtitle(paste0(TP, " (", TP_old, ")")) +
        theme(
          axis.text = element_text(size = 22, color = "black"),
          axis.title.x = element_text(size = 22, vjust = 0),
          axis.title.y = element_text(size = 22, vjust = 2),
          legend.position = "bottom",
          legend.text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 25),
          axis.ticks.length = unit(.25, "cm"),
        ) 
      
      ggsave(plot = p, filename = paste0("08_SelectedPictures Korrelation ProteomicsNRE2 mit KlinischenDaten/", TP_old, "_", TP, "__",
                                         KD_new, "__", Protein, ".png"), 
             width = 21, height = 15, dpi = 200, units = "cm")
      
      print(paste0(TP, " -- ", KD_new, " -- ", Protein))
      
    }
  }
}


