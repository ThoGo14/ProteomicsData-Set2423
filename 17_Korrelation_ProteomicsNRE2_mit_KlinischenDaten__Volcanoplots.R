library(openxlsx)
library(tidyverse)
library(rlang)
library(scales)
library(ggrepel)
library(ggforce)
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

## Load Limma-analysis table for FC
load("01_ProteomicsNRE2_Limma_Analyse.RData")
Limma_FC <- final_data_complete[, c(3, grep("__FoldChange$", colnames(final_data_complete)))] %>% 
  # rename(Protein = Gene_name) %>% 
  mutate_if(is.numeric, round, 2)


Which_FC <- 1.2
Limma_FC_Cols <- colSums(abs(Limma_FC[, -1]) >= Which_FC)
colnames(Limma_FC)[-1] <- paste0(colnames(Limma_FC)[-1], "__(FC", Which_FC, "=", Limma_FC_Cols, ")")

Limma_FC_Rows <- lapply(Limma_FC[, -1], function(x) which(abs(x) >= Which_FC))
names(Limma_FC_Rows) <- sapply(names(Limma_FC_Rows), function (x){
  x <- str_replace_all(x, "__", " ")
})




#### Select Mitochondrion (GO:0005739) -----
Proteomics_Protein_Info_GO <- ProteomicsDataNRE2[["Proteomics_Protein_Info_GO"]]

ProteomicsDataNRE2_Mitochondrion_GO <- Proteomics_Protein_Info_GO %>% 
  filter(GO.Domain == "cellular_component", GO.ID == "GO:0005739") %>% 
  select(Gene_name, GO.Term) %>% 
  distinct(Gene_name, .keep_all = TRUE)
  
ProteomicsDataNRE2_Mitochondrion_GO <- rename(ProteomicsDataNRE2_Mitochondrion_GO,
                                              Protein = Gene_name,
                                              !!paste0("Mitochondrion (GO:0005739) (",
                                                     nrow(ProteomicsDataNRE2_Mitochondrion_GO), ")") := GO.Term)


#### Select Mitochondria from MitoCarta -----

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
  



#### Select Mitochondria from MitoEvidenceIMPI -----

MitoEvidenceIMPI <- read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/MitoEvidenceIMPI/impi-2021-q4pre-20211001-dist_0.xlsx",
                       sheet = "IMPI-2021Q4pre")

MitoEvidenceIMPI$HGNC.Symbol <- ifelse(is.na(MitoEvidenceIMPI$HGNC.Symbol),
                                       MitoEvidenceIMPI$Symbol, MitoEvidenceIMPI$HGNC.Symbol)

MitoEvidenceIMPI_short <- MitoEvidenceIMPI %>% 
  select(HGNC.Symbol, IMPI.Class) %>% 
  filter(HGNC.Symbol %in% Proteomics_Protein_Info$Gene_name) %>% 
  distinct(HGNC.Symbol, .keep_all = T)

MitoEvidenceIMPI_short <- rename(MitoEvidenceIMPI_short, 
                                 Protein = HGNC.Symbol,
                                 !!paste0("MitoEvidence_IMPI.Class (", 
                                          nrow(MitoEvidenceIMPI_short), ")")  := IMPI.Class)

Mito_Proteins <- unique(c(ProteomicsDataNRE2_Mitochondrion_GO$Protein,
                  final_data_complete[final_data_complete$Uniprot_ID %in% MitoCarta_short$Uniprot_ID, "Gene_name"]))

## Korrelationsanalyse ----
Timepoints <- c("T2T0", "T4T5", "T5T0", "T4T2", "T5T2", "T4T0")
Welche_KD <- c("ISI_FC", "SM_MODP.S_FC", "IAS_ergo_BW_FC")
All_Proteins <- Proteomics_Protein_Info$Uniprot_ID

KorrelationTable_Results <- Proteomics_Protein_Info
colnames(KorrelationTable_Results)[2] <- "Protein"

KorrelationTable_Results[, "Limma__Analysis"] <- NA
KorrelationTable_Results <- left_join(KorrelationTable_Results, 
                                      Limma_FC,
                                      by = "Uniprot_ID")

KorrelationTable_Results[, "Correlation__Analysis"] <- NA

for (KD in Welche_KD){
  KorrelationTable_Results[, paste0(KD, "__slope")] <- NA
  KorrelationTable_Results[, paste0(KD, "__intercept")] <- NA
  KorrelationTable_Results[, paste0(KD, "__pValue")] <- NA
  KorrelationTable_Results[, paste0(KD, "__pearsonR")] <- NA
  
}

KorrelationTable_Results_List <- list()

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
      RowNR <- grep(paste0("^", Protein, "$"), KorrelationTable_Results$Uniprot_ID)
      ColNR_slope <- grep(paste0(KD, "__slope"), colnames(KorrelationTable_Results))
      ColNR_intercept <- grep(paste0(KD, "__intercept"), colnames(KorrelationTable_Results))
      ColNR_pValue <- grep(paste0(KD, "__pValue"), colnames(KorrelationTable_Results))
      ColNR_pearson <- grep(paste0(KD, "__pearsonR"), colnames(KorrelationTable_Results))
      
      
      
      Correlation <- cor.test(Table_For_Korrelation[, Protein],
                              Table_For_Korrelation[, KD])
      
      LM1 <- lm(Table_For_Korrelation[, Protein] ~ Table_For_Korrelation[, KD])
      
      
      
      KorrelationTable_Results[RowNR, ColNR_slope] <- LM1[["coefficients"]][2]
      KorrelationTable_Results[RowNR, ColNR_intercept] <- LM1[["coefficients"]][1]
      KorrelationTable_Results[RowNR, ColNR_pValue] <- Correlation[["p.value"]]
      KorrelationTable_Results[RowNR, ColNR_pearson] <- Correlation[["estimate"]][["cor"]]
    }
  }
  
  KorrelationTable_Results_List[[TP]] <- KorrelationTable_Results
  print(paste0(format(Sys.time(), "%x %T"), "  Finished:  ", TP))
}




## Save Data ----

save(KorrelationTable_Results_List, 
     file = "17_Korrelation_ProteomicsNRE2_mit_KlinischenDaten__Volcanoplots.RData")

## functions ----
# function for reversing and Log Transform x/y-axis
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
# function for filtering data inside ggplot
pick <- function(condition){
  function(d) d %>% filter_(condition)
}

#### PlotData -----
# Add color
colors <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 


for (TP in Timepoints) {
  
  Table_for_plot <- KorrelationTable_Results_List[[TP]]
  
  for (KD in Welche_KD) {
    
    Slope_Col = paste0(KD, "__slope")
    Interc_Col = paste0(KD, "__intercept")
    pval_Col = paste0(KD, "__pValue")
    pearson_Col = paste0(KD, "__pearsonR")
    
    RegulatedGenes_Count <- Table_for_plot %>%
      select(Protein, all_of(Slope_Col), all_of(Interc_Col), all_of(pval_Col), all_of(pearson_Col)) %>% 
      mutate(Regulation = case_when((!!parse_expr(Slope_Col) >= 0.2 & !!parse_expr(pval_Col) <= 0.05) ~ "up",
                                     (!!parse_expr(Slope_Col) <= -0.2 & !!parse_expr(pval_Col) <= 0.05) ~ "down",
                                     TRUE ~ "ns")) %>% 
      count(Regulation) %>% 
      mutate(Subtitles = paste0(Regulation, " (", n, ")"))
    
    # Subtitle <- RegulatedGenes_Count$Subtitles
    Subtitle <- RegulatedGenes_Count$n
    names(Subtitle) = RegulatedGenes_Count$Regulation
    
    Ratios = Table_for_plot %>% select(Protein, all_of(Slope_Col), all_of(pval_Col)) %>% 
      arrange(-!!parse_expr(Slope_Col)) %>% filter(!!parse_expr(pval_Col) <= 0.05) %>%  head(5)
    
    MaxVal = Table_for_plot %>% select(all_of(Slope_Col)) %>% max() %>% ceiling()
    
    
   Table_for_plot %>% 
      select(Protein, all_of(Slope_Col), all_of(Interc_Col), all_of(pval_Col), all_of(pearson_Col)) %>% 
      mutate(Regulation = case_when(!!parse_expr(Slope_Col) > 0.2 & !!parse_expr(pval_Col) <= 0.05 ~ "up",
                                    !!parse_expr(Slope_Col) < -0.2 & !!parse_expr(pval_Col) <= 0.05 ~ "down",
                                    TRUE ~ "ns")) %>%
      ggplot(aes(x = .data[[Slope_Col]],
                 y = .data[[pval_Col]])) +
      geom_point(aes(fill = Regulation), shape = 21, color = "black", size = 2) +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      geom_vline(xintercept = 0.2, linetype = "dashed") +
      geom_vline(xintercept = -0.2, linetype = "dashed") +
      geom_label_repel(data = Ratios, 
                       aes(alpha = "Red", label = Protein),
                       box.padding   = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'black',
                       min.segment.length = 0,
                       max.overlaps = 20,
                       force = 20,
                       seed = 1,
                       color = "black",
                       direction = "y",
                       nudge_x = MaxVal - Ratios[, Slope_Col],
                       alpha = 1) +
      scale_y_continuous(trans = reverselog_trans(10)) +
      scale_fill_manual("", values = colors, breaks = RegulatedGenes_Count$Regulation, labels = Subtitle) +
      labs(x ="Ratio (slope)",
           y = "pValue",
           title = paste(TP, " -- ", KD)) +
      theme(
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 20, vjust = 2),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(linetype = "solid", size = 0.5, color = "black"),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(face = "bold", size = 13),
        axis.ticks.length=unit(.25, "cm"),
      ) +
      guides(
        fill = guide_legend(override.aes = list(size = 3))
      )
      
    print(paste(TP, " -- ", KD))
    
    ggsave(filename = paste0("17_Volcanoplots from Korrelation mit KlinischeDaten/", "Volcanoplot",
                             TP, "__", KD, ".png"), 
           width = 20, height = 15, units = "cm", dpi = 200)
      
  }
}


#### PlotData with Mitos-----
# Add color
colors <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
shape <- c("Mito" = 22, "Non Mito" = 21)
size <- c("Mito" = 2.5, "Non Mito" = 2)


for (TP in Timepoints) {
  
  Table_for_plot <- KorrelationTable_Results_List[[TP]]
  
  for (KD in Welche_KD) {
    
    Slope_Col = paste0(KD, "__slope")
    Interc_Col = paste0(KD, "__intercept")
    pval_Col = paste0(KD, "__pValue")
    pearson_Col = paste0(KD, "__pearsonR")
    
    
    Table_for_plot_mutated <-  Table_for_plot %>%
      select(Protein, all_of(Slope_Col), all_of(Interc_Col), all_of(pval_Col), all_of(pearson_Col)) %>% 
      mutate(Regulation = case_when((!!parse_expr(Slope_Col) >= 0.2 & !!parse_expr(pval_Col) <= 0.05) ~ "up",
                                    (!!parse_expr(Slope_Col) <= -0.2 & !!parse_expr(pval_Col) <= 0.05) ~ "down",
                                    TRUE ~ "ns"),
             Mitochondrion = case_when(Protein %in% Mito_Proteins ~ "Mito",
                                       TRUE ~ "Non Mito"),
             Non_Mito_up = case_when(Regulation == "up" & Mitochondrion == "Non Mito" ~ "NonMitoUp",
                                     TRUE ~ "Not"))
    
    RegulatedGenes_Count <- Table_for_plot_mutated %>% 
      count(Regulation) %>% 
      mutate(Subtitles = paste0(Regulation, " (", n, ")"))
    
    # Subtitle <- RegulatedGenes_Count$Subtitles
    Subtitle <- RegulatedGenes_Count$n
    names(Subtitle) = RegulatedGenes_Count$Regulation
    
    Non_Mito_up <- Table_for_plot_mutated %>%
      filter(Regulation == "up", Mitochondrion == "Non Mito") %>%
      nrow()
    
    
    if (Non_Mito_up > 0) {
      if (Non_Mito_up > 10) {
        
        Proteins = Table_for_plot_mutated %>% 
          filter(., Non_Mito_up == "NonMitoUp") %>% 
          mutate(MostImportant = !!parse_expr(Slope_Col) - !!parse_expr(pval_Col)) %>% 
          arrange(-MostImportant) %>% 
          mutate(MostImportant_Select = c(rep(TRUE, times = 10), rep(NA, times = nrow(.)-10))) %>% 
          filter(MostImportant_Select == TRUE) %>% pull(Protein)
        
        Ratios_df = Table_for_plot_mutated %>% filter(Protein %in% Proteins) %>% 
          select(Protein, all_of(Slope_Col))
        
        Ratios = Ratios_df[, 2]
        names(Ratios) = Ratios_df[, 1]
        
      } else {
        Ratios = Table_for_plot_mutated %>% 
          filter(., Non_Mito_up == "NonMitoUp") %>% pull(all_of(Slope_Col))
        
      }
      MaxVal = Table_for_plot %>% select(all_of(Slope_Col)) %>% max() %>% ceiling()
    }
    
    
    Table_for_plot_mutated %>%
      ggplot(aes(x = .data[[Slope_Col]],
                 y = .data[[pval_Col]])) +
      geom_point(aes(fill = Regulation, shape = Mitochondrion, size = Mitochondrion), color = "black") + {
        if (Non_Mito_up > 0) 
          geom_point(data = . %>% filter(Non_Mito_up == "NonMitoUp"), 
                     aes(alpha = "Red"), size = 2, shape = 21, fill = "red", color = "black")
      } + {
        if (Non_Mito_up > 0) 
          scale_alpha_manual("", breaks = c("Red"), values = 1, labels = "Non Mito & up")
      } +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      geom_vline(xintercept = 0.2, linetype = "dashed") +
      geom_vline(xintercept = -0.2, linetype = "dashed") + {
        if (Non_Mito_up > 10)
          geom_label_repel(data = . %>% filter(Protein %in% Proteins), 
                           aes(alpha = "Red", label = Protein),
                           box.padding   = 0.35, 
                           point.padding = 0.5,
                           segment.color = 'black',
                           min.segment.length = 0,
                           max.overlaps = 20,
                           force = 20,
                           seed = 1,
                           color = "black",
                           direction = "y",
                           nudge_x = MaxVal - Ratios,
                           alpha = 1)
        
      } + {
        if (Non_Mito_up > 0 & Non_Mito_up <= 10)
          geom_label_repel(data = . %>% filter(Non_Mito_up == "NonMitoUp"), 
                           aes(alpha = "Red", label = Protein),
                           box.padding   = 0.35, 
                           point.padding = 0.5,
                           segment.color = 'black',
                           min.segment.length = 0,
                           max.overlaps = 20,
                           force = 20,
                           seed = 1,
                           color = "black",
                           direction = "y",
                           nudge_x = MaxVal - Ratios,
                           alpha = 1)
        
      } +
      scale_y_continuous(trans = reverselog_trans(10)) +
      scale_fill_manual("", values = colors, breaks = RegulatedGenes_Count$Regulation, labels = Subtitle) +
      scale_shape_manual("", values = shape) +
      scale_size_manual("", values = size) +
      labs(x ="Slope linear model",
           y = "pValue",
           title = paste(TP, " -- ", 
                         ifelse(KD == "SM_MODP.S_FC",
                                "Maximal coupled respiration", KD))) +
      theme(
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 20, vjust = 2),
        axis.title.y = element_text(size = 20, vjust = 2),
        legend.position = "top",
        legend.justification = "left",
        legend.direction = "horizontal", 
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"),
        legend.background = element_rect(linetype = "solid", size = 0.5, color = "black"),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(face = "bold", size = 13),
        axis.ticks.length=unit(.25, "cm"),
      ) +
      guides(
        # fill = guide_legend(override.aes = list(shape = 21, size = 3), 
        #                     order = 1,
        #                     ),
        fill = guide_none(),
        shape = guide_legend(override.aes = list(shape = c("Mito" = 0, "Non Mito" = 1),
                                                 size = 3, stroke = 2),
                             order = 2),
        # alpha = guide_legend(override.aes = list(shape = c("NonMitoUp" = 21), size = c("NonMitoUp" = 3), 
        #                                          fill = "red", color = "black"),
        #                      order = 3),
        alpha = guide_none(),
        size = guide_none()
      )
    
    print(paste(TP, " -- ", KD))
    
    ggsave(filename = paste0("17_Volcanoplots from Korrelation mit KlinischeDaten/", "Volcanoplot",
                             TP, "__", KD, "__Mitochondrion", ".png"), 
           width = 20, height = 15, units = "cm", dpi = 200)
    
  }
}
