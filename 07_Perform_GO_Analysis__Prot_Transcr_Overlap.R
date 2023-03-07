library(openxlsx)
library(tidyverse)
library(gprofiler2)
theme_set(theme_bw())

## Path where Input File is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

load("06_Korrelation_ProteomicsNRE2_mit_TranscriptomicsNRE2.RData")

PVal_cols <- colnames(KorrelationTable_Results)[grep("pValue", colnames(KorrelationTable_Results))]
# PVal_cols <- PVal_cols[!str_detect(PVal_cols, "adj.pValue")]

## functions
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- paste0("10^-",l)
  # return this as an expression
  parse(text=l)
}

## gProfiler for both directions combined ---- 
for (i in PVal_cols) {
  Uniprot_ID <- KorrelationTable_Results[KorrelationTable_Results[, i] < 0.05, "Uniprot_ID"]
  
  
  gostres <- gost(query = Uniprot_ID,
                  organism = "hsapiens",
                  user_threshold = 0.05,
                  correction_method = "gSCS",
                  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                  evcodes = TRUE)
  if (!is.null(gostres)) {
    
    Results <- gostres$result
    
    Results$intersection <- sapply(Results$intersection, function(x){
      x <- unlist(str_split(x, pattern = ","))
      x <- gconvert(x, organism = "hsapiens")$name
      x <- paste0(x, collapse = ",")
    })
    
    Results <- Results %>% 
      select(source, term_name, term_id, p_value, term_size, query_size, intersection_size, effective_domain_size, intersection) %>% 
      mutate(negative_log10_of_adjusted_p_value = -log10(p_value),
             TermAndID = paste(term_name, " (", term_id, ")", sep = ""),
             AnnotatedGenes = cut(intersection_size,
                                  breaks = c(0,51,151, Inf),
                                  labels = c("<50", "50 - 150" , ">150")),
             AnnotatedGenesPercent = round(intersection_size / term_size *100, digits = 10))
    
    write.csv(Results, file = paste0("07_GO Analysis Pictures - Prot-Transcr-Overlap/", "gProfiler__", i, ".csv"))
    
    
    
    Top5_Source <- NULL
    for (Source in unique(Results$source)) {
      Top5 <- Results %>% filter(source %in% Source)
      Top5 <- Top5[1:5, "TermAndID"]
      
      Top5_Source <- c(Top5_Source, Top5)
    }
    
    Results_new <- Results %>% 
      filter(TermAndID %in% Top5_Source) %>% 
      arrange(match(source, c("REAC", "KEGG", "GO:MF", "GO:BP", "GO:CC")), AnnotatedGenesPercent) %>% 
      select(source, term_name, term_id, TermAndID, AnnotatedGenesPercent, p_value, negative_log10_of_adjusted_p_value)
    
    for (RowNR in 1:nrow(Results_new)) {
      if (nchar(Results_new[RowNR, "TermAndID"]) > 85) {
        TermID_Lenght <- nchar(Results_new[RowNR, "term_id"])
        New_Text_lenght <- 85 - TermID_Lenght - 6
        TermAndID_new <- paste(substr(Results_new[RowNR, "term_name"], 1, New_Text_lenght), "... (",
                               Results_new[RowNR, "term_id"], ")", sep = "")
        
        Results_new[RowNR, "TermAndID"] <- TermAndID_new
      }
    }
    
    Results_new$TermAndID <- factor(Results_new$TermAndID,
                                    levels = unique(Results_new$TermAndID))
    
    
    
    
    vlines <- table(Results_new$source)
    vlines <- vlines[c("REAC", "KEGG", "GO:MF", "GO:BP", "GO:CC")]
    vlines <- vlines[!is.na(names(vlines))]
    vlines <- cumsum(vlines)+0.5
    vlines <- vlines[-length(vlines)]
    
    y_line_color <- c("white", "white", "white", "white", "black")
    
    Results_new %>%   
      ggplot(aes(x = TermAndID,
                 y = AnnotatedGenesPercent)) +
      geom_point(aes(color = source, size = negative_log10_of_adjusted_p_value)) +
      ylab("Annotated Genes \n[% of all genes/term]")+
      xlab("") +
      geom_vline(xintercept = vlines, color = "lightgray") +
      coord_flip() +
      scale_color_manual("",
                         breaks = c("GO:CC", "GO:BP", "GO:MF", "KEGG", "REAC"),
                         labels = c("GO:CC", "GO:BP", "GO:MF", "KEGG", "REAC"),
                         values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) +
      scale_size("p-Value",
                 labels = fancy_scientific) +
      ggtitle(paste0(i, "\n(n = ", length(Uniprot_ID), ")")) +
      theme(
        axis.ticks.length = unit(0.25, units = "cm"),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 8, color = "black"),
        legend.text.align = 0,
        legend.title = element_text(size = 9, face = "bold"),
        title = element_text(size = 9, face = "bold"), 
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(),
        panel.grid.major.y = element_blank(),
      ) +
      guides(
        color = guide_legend(override.aes = list(size=4), order = 1),
        size = guide_legend(order = 2)
      )
    
    ggsave(paste0("07_GO Analysis Pictures - Prot-Transcr-Overlap/", "gProfiler__", i, ".png"), dpi = 300, units = "cm", 
           width = 24, height = 15)
  }
  
  print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", i))
}


## gProfiler for both directions separated ---- 
for (i in PVal_cols) {
  
  TP = unlist(str_split(i, "__"))[1]
  
  for (Direction in c("positiveCorrelation", "negativeCorrelation")) {
    
    
    if(Direction == "positiveCorrelation"){
      Uniprot_ID <- KorrelationTable_Results[KorrelationTable_Results[, i] < 0.05 &
                                               KorrelationTable_Results[, paste0(TP, "__pearsonR")] > 0,
                                             "Uniprot_ID"]
    } else if(Direction == "negativeCorrelation"){
      Uniprot_ID <- KorrelationTable_Results[KorrelationTable_Results[, i] < 0.05 &
                                               KorrelationTable_Results[, paste0(TP, "__pearsonR")] < 0,
                                             "Uniprot_ID"]
    }
    
    
    gostres <- gost(query = Uniprot_ID,
                    organism = "hsapiens",
                    user_threshold = 0.05,
                    correction_method = "gSCS",
                    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                    evcodes = TRUE)
    if (!is.null(gostres)) {
      
      Results <- gostres$result
      
      Results$intersection <- sapply(Results$intersection, function(x){
        x <- unlist(str_split(x, pattern = ","))
        x <- gconvert(x, organism = "hsapiens")$name
        x <- paste0(x, collapse = ",")
      })
      
      Results <- Results %>% 
        select(source, term_name, term_id, p_value, term_size, query_size, intersection_size, effective_domain_size, intersection) %>% 
        mutate(negative_log10_of_adjusted_p_value = -log10(p_value),
               TermAndID = paste(term_name, " (", term_id, ")", sep = ""),
               AnnotatedGenes = cut(intersection_size,
                                    breaks = c(0,51,151, Inf),
                                    labels = c("<50", "50 - 150" , ">150")),
               AnnotatedGenesPercent = round(intersection_size / term_size *100, digits = 10))
      
      write.csv(Results, file = paste0("07_GO Analysis Pictures - Prot-Transcr-Overlap/", 
                                       "gProfiler__", i, ifelse(Direction == "positiveCorrelation",
                                                                "__positiveCorrelation",
                                                                "__negativeCorrelation"),
                                       ".csv"))
      
      
      
      Top5_Source <- NULL
      for (Source in unique(Results$source)) {
        Top5 <- Results %>% filter(source %in% Source)
        Top5 <- Top5[1:5, "TermAndID"]
        
        Top5_Source <- c(Top5_Source, Top5)
      }
      
      Results_new <- Results %>% 
        filter(TermAndID %in% Top5_Source) %>% 
        arrange(match(source, c("REAC", "KEGG", "GO:MF", "GO:BP", "GO:CC")), AnnotatedGenesPercent) %>% 
        select(source, term_name, term_id, TermAndID, AnnotatedGenesPercent, p_value, negative_log10_of_adjusted_p_value)
      
      for (RowNR in 1:nrow(Results_new)) {
        if (nchar(Results_new[RowNR, "TermAndID"]) > 85) {
          TermID_Lenght <- nchar(Results_new[RowNR, "term_id"])
          New_Text_lenght <- 85 - TermID_Lenght - 6
          TermAndID_new <- paste(substr(Results_new[RowNR, "term_name"], 1, New_Text_lenght), "... (",
                                 Results_new[RowNR, "term_id"], ")", sep = "")
          
          Results_new[RowNR, "TermAndID"] <- TermAndID_new
        }
      }
      
      Results_new$TermAndID <- factor(Results_new$TermAndID,
                                      levels = unique(Results_new$TermAndID))
      
      
      
      
      vlines <- table(Results_new$source)
      vlines <- vlines[c("REAC", "KEGG", "GO:MF", "GO:BP", "GO:CC")]
      vlines <- vlines[!is.na(names(vlines))]
      vlines <- cumsum(vlines)+0.5
      vlines <- vlines[-length(vlines)]
      
      y_line_color <- c("white", "white", "white", "white", "black")
      
      Results_new %>%   
        ggplot(aes(x = TermAndID,
                   y = AnnotatedGenesPercent)) +
        geom_point(aes(color = source, size = negative_log10_of_adjusted_p_value)) +
        ylab("Annotated Genes \n[% of all genes/term]")+
        xlab("") +
        geom_vline(xintercept = vlines, color = "lightgray") +
        coord_flip() +
        scale_color_manual("",
                           breaks = c("GO:CC", "GO:BP", "GO:MF", "KEGG", "REAC"),
                           labels = c("GO:CC", "GO:BP", "GO:MF", "KEGG", "REAC"),
                           values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) +
        scale_size("p-Value",
                   labels = fancy_scientific) +
        ggtitle(paste0(i, "\n", Direction, " (n = ", length(Uniprot_ID), ")")) +
        theme(
          axis.ticks.length = unit(0.25, units = "cm"),
          axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 8, color = "black"),
          legend.text.align = 0,
          legend.title = element_text(size = 9, face = "bold"),
          title = element_text(size = 9, face = "bold"), 
          panel.grid.major.x = element_line(),
          panel.grid.minor.x = element_line(),
          panel.grid.major.y = element_blank(),
        ) +
        guides(
          color = guide_legend(override.aes = list(size=4), order = 1),
          size = guide_legend(order = 2)
        )
      
      ggsave(paste0("07_GO Analysis Pictures - Prot-Transcr-Overlap/", 
                    "gProfiler__", i, ifelse(Direction == "positiveCorrelation",
                                             "__positiveCorrelation",
                                             "__negativeCorrelation"), ".png"), 
             dpi = 300, units = "cm", width = 24, height = 15)
    }
  }
  print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", i))
}


