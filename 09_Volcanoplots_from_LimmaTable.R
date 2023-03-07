library(tidyverse)
library(scales)
library(ggrepel)
library(ggforce)
library(rlang)
library(openxlsx)

theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("01_ProteomicsNRE2_Limma_Analyse.RData")
LimmaAnalysisTable <- final_data_complete
rm(final_data_complete)

load("ProteomicsNRE2__Expression.RData")

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
MitoCarta_short$Protein <- LimmaAnalysisTable[grep(paste0(MitoCarta_short$Uniprot_ID, collapse = "|"), LimmaAnalysisTable$Uniprot_ID),
                                              "Gene_name"]


## Combine both Tables

Mito_Proteins <- unique(c(ProteomicsDataNRE2_Mitochondrion_GO$Protein,
                           MitoCarta_short$Protein))


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


for(TP in c("T5T0", "T2T0", "T4T5", "T4T2", "T5T2", "T4T0")){
  for (pval in c("__pValue", "__adj.pValue")) {
    
    Ratio_col <- paste0(TP, "_FC__Ratio")
    pval_col <- paste0(TP, "_FC__Ratio", pval)
    
    ## Explanation
    # first parse your expressions to a quosure with parse_expr from rlang (parse_quosure is depricted)
    # then unquote using !! (or UQ())
    # A quosure is an object which contains an expression and an environment
    
    RegulatedGenes_Count <- LimmaAnalysisTable %>% 
      select(Gene_name, Uniprot_ID, starts_with(Ratio_col)) %>% 
      mutate(Regulation = case_when(!!parse_expr(Ratio_col) >= 1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "up",
                                   !!parse_expr(Ratio_col) <= 1/1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "down",
                                   TRUE ~ "ns")) %>% 
      count(Regulation) %>% 
      mutate(Subtitles = paste0(Regulation, " (", n, ")"))
    
    # Subtitle <- RegulatedGenes_Count$Subtitles
    Subtitle <- RegulatedGenes_Count$n
    names(Subtitle) = RegulatedGenes_Count$Regulation
    
    MaxVal = LimmaAnalysisTable %>% select(all_of(Ratio_col)) %>% max() %>% ceiling()
    
    Ratios = LimmaAnalysisTable %>% select(Gene_name, all_of(Ratio_col), all_of(pval_col)) %>% 
      arrange(-!!parse_expr(Ratio_col)) %>% filter(!!parse_expr(pval_col) <= 0.05) %>%  head(5)

    LimmaAnalysisTable %>% 
      select(Gene_name, Uniprot_ID, starts_with(Ratio_col)) %>% 
      mutate(Regulation = case_when(!!parse_expr(Ratio_col) >= 1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "up",
                                   !!parse_expr(Ratio_col) <= 1/1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "down",
                                   TRUE ~ "ns")) %>% 
      ggplot(aes(x = .data[[Ratio_col]],
                 y = .data[[pval_col]],
                 )) +
      geom_point(aes(fill = Regulation), shape = 21, color = "black", size = 2) +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      # geom_vline(xintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = 1.2, linetype = "dashed") +
      geom_vline(xintercept = 1/1.2, linetype = "dashed") +
      geom_label_repel(data = Ratios, 
                       aes(alpha = "Red", label = Gene_name),
                       box.padding   = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'black',
                       min.segment.length = 0,
                       max.overlaps = 10,
                       force = 15,
                       seed = 1,
                       color = "black",
                       direction = "y",
                       nudge_x = MaxVal - Ratios[, Ratio_col],
                       alpha = 1) +
      labs(x ="Ratio",
           y = ifelse(pval == "__adj.pValue", "adj. pValue (BH)", "pValue"),
           title = TP,
           # subtitle = Subtitle,
           caption = ifelse(TP %in% c("T4T5", "T4T2", "T4T0"), "n = 24", "n = 25")) +
      scale_y_continuous(trans = reverselog_trans(10)) +
      scale_fill_manual("", values = colors, breaks = RegulatedGenes_Count$Regulation, labels = Subtitle) +
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
    
    print(paste0(TP, pval))
    
    ggsave(filename = paste0("09_Volcanoplots from LimmaTable/", "Volcanoplot",
                             TP, "__", pval, ".png"), 
           width = 20, height = 15, units = "cm", dpi = 200)
  }
}

#### PlotData with Mitos-----
# Add color
colors <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
shape <- c("Mito" = 22, "Non Mito" = 21)
size <- c("Mito" = 2.5, "Non Mito" = 2)


# Old Colors
# Non Mito & up: #BF5204


for(TP in c("T5T0", "T2T0", "T4T5", "T4T2", "T5T2", "T4T0")){
  for (pval in c("__pValue", "__adj.pValue")) {
    
    Ratio_col <- paste0(TP, "_FC__Ratio")
    pval_col <- paste0(TP, "_FC__Ratio", pval)
    
    ## Explanation
    # first parse your expressions to a quosure with parse_expr from rlang (parse_quosure is depricted)
    # then unquote using !! (or UQ())
    # A quosure is an object which contains an expression and an environment
    
    LimmaAnalysisTable_mutated <- LimmaAnalysisTable %>% 
      mutate(Regulation = case_when(!!parse_expr(Ratio_col) >= 1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "up",
                                    !!parse_expr(Ratio_col) <= 1/1.2 &  !!parse_expr(pval_col) <= 0.05 ~ "down",
                                    TRUE ~ "ns"),
             Mitochondrion = case_when(Gene_name %in% Mito_Proteins ~ "Mito",
                                       TRUE ~ "Non Mito"),
             Non_Mito_up = case_when(Regulation == "up" & Mitochondrion == "Non Mito" ~ "NonMitoUp",
                                     TRUE ~ "Not")) 
    
    RegulatedGenes_Count <- LimmaAnalysisTable_mutated %>% 
      select(Gene_name, Uniprot_ID, starts_with(Ratio_col), Regulation) %>% 
      count(Regulation) %>% 
      mutate(Subtitles = paste0(Regulation, " (", n, ")"))
    
    # Subtitle <- RegulatedGenes_Count$Subtitles
    Subtitle <- RegulatedGenes_Count$n
    names(Subtitle) = RegulatedGenes_Count$Regulation
    
    Non_Mito_up <- LimmaAnalysisTable_mutated %>%
      filter(Regulation == "up", Mitochondrion == "Non Mito") %>%
      nrow()
    
    if (Non_Mito_up > 0) {
      if (Non_Mito_up > 10) {
        
        Uniprot_IDs = LimmaAnalysisTable_mutated %>% 
          filter(., Non_Mito_up == "NonMitoUp") %>% 
          mutate(MostImportant = !!parse_expr(Ratio_col) - !!parse_expr(pval_col)) %>% 
          arrange(-MostImportant) %>% 
          mutate(MostImportant_Select = c(rep(TRUE, times = 10), rep(NA, times = nrow(.)-10))) %>% 
          filter(MostImportant_Select == TRUE) %>% pull(Uniprot_ID)
        
        Ratios_df = LimmaAnalysisTable_mutated %>% filter(Uniprot_ID %in% Uniprot_IDs) %>% 
          select(Uniprot_ID, all_of(Ratio_col))
        
        Ratios = Ratios_df[, 2]
        names(Ratios) = Ratios_df[, 1]
        
      } else {
        Ratios = LimmaAnalysisTable_mutated %>% 
          filter(., Non_Mito_up == "NonMitoUp") %>% pull(all_of(Ratio_col))
        
      }
      MaxVal = LimmaAnalysisTable %>% select(all_of(Ratio_col)) %>% max() %>% ceiling()
    }
    
    LimmaAnalysisTable_mutated %>% 
      ggplot(aes(x = .data[[Ratio_col]],
                 y = .data[[pval_col]])) +
      geom_point(aes(fill = Regulation, shape = Mitochondrion, size = Mitochondrion), color = "black") + {
        if (Non_Mito_up > 0) 
          geom_point(data = . %>% filter(Non_Mito_up == "NonMitoUp"), 
                     aes(alpha = "Red"), size = 2, shape = 21, fill = "red", color = "black")
      } + {
        if (Non_Mito_up > 0) 
          scale_alpha_manual("", breaks = c("Red"), values = 1, labels = "Non Mito & up")
      } + 
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      # geom_vline(xintercept = 1, linetype = "dashed") +
      geom_vline(xintercept = 1.2, linetype = "dashed") +
      geom_vline(xintercept = 1/1.2, linetype = "dashed") + {
        if (Non_Mito_up > 10)
          geom_label_repel(data = . %>% filter(Uniprot_ID %in% Uniprot_IDs), 
                           aes(alpha = "Red", label = Gene_name),
                           box.padding   = 0.35, 
                           point.padding = 0.5,
                           segment.color = 'black',
                           min.segment.length = 0,
                           max.overlaps = 10,
                           force = 15,
                           seed = 1,
                           color = "black",
                           direction = "y",
                           nudge_x = MaxVal - Ratios,
                           alpha = 1)
        
      } + {
        if (Non_Mito_up > 0 & Non_Mito_up <= 10)
          geom_label_repel(data = . %>% filter(Non_Mito_up == "NonMitoUp"), 
                           aes(alpha = "Red", label = Gene_name),
                           box.padding   = 0.35, 
                           point.padding = 0.5,
                           segment.color = 'black',
                           min.segment.length = 0,
                           max.overlaps = 10,
                           force = 15,
                           seed = 1,
                           color = "black",
                           direction = "y",
                           nudge_x = MaxVal - Ratios,
                           alpha = 1)
        
      } +
      labs(x ="Ratio",
           y = ifelse(pval == "__adj.pValue", "adj. pValue (BH)", "pValue"),
           title = TP,
           # subtitle = Subtitle,
           caption = ifelse(TP %in% c("T4T5", "T4T2", "T4T0"), "n = 24", "n = 25")) +
      scale_y_continuous(trans = reverselog_trans(10)) +
      scale_fill_manual("", values = colors, breaks = RegulatedGenes_Count$Regulation, labels = Subtitle) +
      scale_shape_manual("", values = shape) +
      scale_size_manual("", values = size) +
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
    
    
    
    
    print(paste0(TP, pval))
    
    ggsave(filename = paste0("09_Volcanoplots from LimmaTable/", "Volcanoplot",
                             TP, "__", pval, "__Mitochondrion", ".png"), 
           width = 20, height = 15, units = "cm", dpi = 200)
  }
}
