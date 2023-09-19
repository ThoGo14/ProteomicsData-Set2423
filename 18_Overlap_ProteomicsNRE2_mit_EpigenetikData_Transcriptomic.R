library(openxlsx)
library(tidyverse)
library(rlang)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")
load("E:/Auswertung mit R/00_Helpful Files/NRE2_All_TimePoints_Transcripts_FC_no_DABG_remove.RData")

NRE2_GenesTable <- readRDS("E:/Auswertung mit R/00_Helpful Files/NRE2_GenesTable_no_DABG_remove.RDS")

# Proteomics_Proteins <- row.names(ProteomicsDataNRE2[["ProteomicsData__Expression_log2"]])
Proteomics_Proteins <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]$Gene_name

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/00_Klinische_Daten.xlsx")

## Limma Protein
load("01_ProteomicsNRE2_Limma_Analyse.RData")
Limma_Protein <- final_data_complete
rm(final_data_complete)


## Limma Transcript
Limma_Transcripts <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/NRE1&NRE2 Combined MI + Limma/Limma Test - All Time Points/All_Timepoint_NRE2_withResponseGrouping.xlsx")
Limma_Transcripts <- left_join(Limma_Transcripts,
                               NRE2_GenesTable %>% select(Probeset, Unique),
                               by = "Probeset")

# Proteins_notFound <- read.xlsx("06_Proteins_NotFound_Table.xlsx")


NRE2_TimePoints_Baseintensity_noDABG <- readRDS("E:/Auswertung mit R/00_Helpful Files/NRE2_TimePoints_Baseintensity_noDABG.RDS")

## Selecting Proteins with different names ----
table(Proteomics_Proteins %in% NRE2_GenesTable$Gene.Symbol)
Proteins_found <- Proteomics_Proteins[Proteomics_Proteins %in% NRE2_GenesTable$Gene.Symbol]

HGNC_DataTable <- read.delim("E:/Auswertung mit R/00_Gene_Annotation/HGNC/HGNC_BioMart_Full.txt")

NotFound <- Proteomics_Proteins[!Proteomics_Proteins %in% NRE2_GenesTable$Gene.Symbol]
table(NotFound %in% HGNC_DataTable$Approved.symbol)

## Selecting Alias symbols and previous symbols for genes not found
NotFound_New <- list()
for(i in NotFound){
  Any_found = NULL
  
  RowNR = grep(paste0("^", i, "$"), HGNC_DataTable$Approved.symbol)
  
  AliasSymbols <- HGNC_DataTable[RowNR, "Alias.symbol"]
  for (x in AliasSymbols) {
    if (x %in% NRE2_GenesTable$Gene.Symbol) {
      Any_found = c(Any_found, x)
    }
  }
  
  Previous.symbol <- HGNC_DataTable[RowNR, "Previous.symbol"]
  for (x in Previous.symbol) {
    if (x %in% NRE2_GenesTable$Gene.Symbol) {
      Any_found = c(Any_found, x)
    }
  }
  
  
  NotFound_New[[i]] <- unique(Any_found)
}

## Creating Table with Protein Name and Gen Symbol
Proteins_notFound <- data.frame(Protein.Name = NA,
                                Gene.Symbol = NA)
for (i in names(NotFound_New)) {
  x = data.frame(Protein.Name = rep(i, length(NotFound_New[[i]])), 
                 Gene.Symbol = NotFound_New[[i]])
  
  Proteins_notFound <- rbind(Proteins_notFound, x)
}
Proteins_notFound <- Proteins_notFound[-1, ]



Proteins_with_Transcripts <- rbind(Proteins_notFound,
                                   cbind(Protein.Name = Proteins_found, 
                                         Gene.Symbol = Proteins_found))


Timepoints = data.frame(Epigenetik = c("T2T1", "T4T1", "T3T4"),
                        Proteomics = c("T2T0", "T5T0", "T4T5"))

Unsplitted = c("unsplitted_METH_EXP_p001.xlsx", "unsplitted_METH_EXP_p005.xlsx")
Splitted = c("LRE_RES_SPLIT_METH_EXP_p001.xlsx", "LRE_RES_SPLIT_METH_EXP_p005.xlsx")


Proteins_with_Transcripts_TP_list = list() 

for(EpiTable in c(Unsplitted, Splitted)){
  
  if (EpiTable %in% Unsplitted) {
    load("01_ProteomicsNRE2_Limma_Analyse.RData")
    Limma_Protein <- final_data_complete
    rm(final_data_complete)
  } else if (EpiTable %in% Splitted) {
    load("03_Limma_Analyse__Split_RESvsLRE__FC.RData")
    Limma_Protein <- final_data
    rm(final_data)
  }
  
  for(TP in Timepoints$Proteomics){
    
    if (EpiTable %in% Unsplitted) {
      LimmaCols = grep(paste0(TP, "_FC__FoldChange"), colnames(Limma_Protein))
    } else if (EpiTable %in% Splitted) {
      LimmaCols = grep(paste0(TP, "_FC___FoldChange"), colnames(Limma_Protein))
    }
    
    
    Proteins_with_Transcripts_TP = left_join(Proteins_with_Transcripts,
                                             Limma_Protein[, c(2, LimmaCols)] %>% rename(Protein.Name = Gene_name),
                                             by = "Protein.Name")
    
    EpiTable_Names <- names(loadWorkbook(paste0("E:/Auswertung mit R/Epigentik Analysen NRE2/Exel Files Auswertung/", EpiTable)))
    EpiTable_Names_TPs = grep(Timepoints$Epigenetik[Timepoints$Proteomics == TP], EpiTable_Names)
    
    for(EpiTable_Names_TP in EpiTable_Names_TPs){
      
      Epigenetik = read.xlsx(paste0("E:/Auswertung mit R/Epigentik Analysen NRE2/Exel Files Auswertung/", EpiTable),
                             sheet = EpiTable_Names_TP)
      
      
      if (grepl("_cor_dir", EpiTable_Names[EpiTable_Names_TP])) {
        # Proteins_with_Transcripts_TP <- Proteins_with_Transcripts_TP %>% 
        #   mutate(Transkr_Epig__Correct_direction = case_when((Protein.Name %in% Epigenetik$Gene.Symbol ~ "TRUE")))
        
        Correct_Sites = Epigenetik$IlmnID
        
      } else {
        
        Epigenetik_short <- Epigenetik %>% 
          select(Gene.Symbol, 
                 all_of(colnames(Epigenetik)[grep("mean_", colnames(Epigenetik))[1]:grep("IlmnID", colnames(Epigenetik))]),
                 UCSC_RefGene_Group_short) %>% 
          rename(Protein.Name = Gene.Symbol) 
        
        Proteins_with_Transcripts_TP = Proteins_with_Transcripts_TP %>% 
          filter(Protein.Name %in% Epigenetik_short$Protein.Name)
        
        Proteins_with_Transcripts_TP = left_join(Proteins_with_Transcripts_TP, 
                                                 Epigenetik_short,
                                                 by = "Protein.Name")
      }
    }
    Delta_col = grep("delta_", colnames(Proteins_with_Transcripts_TP), value = TRUE)
    
    if (EpiTable %in% Unsplitted) {
      Proteins_with_Transcripts_TP <- Proteins_with_Transcripts_TP %>%
        mutate(
          Transkr_Epig__Correct_direction = case_when(IlmnID %in% Correct_Sites ~ TRUE),
          Protein_Epig__Correct_direction = case_when(
            !!parse_expr(paste0(
              TP,
              "_FC__FoldChange"
            )) > 0 &
              !!parse_expr(Delta_col) > 0 &
              UCSC_RefGene_Group_short == "Body" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC__FoldChange"
            )) < 0 &
              !!parse_expr(Delta_col) < 0 &
              UCSC_RefGene_Group_short == "Body" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC__FoldChange"
            )) > 0 &
              !!parse_expr(Delta_col) < 0 &
              UCSC_RefGene_Group_short == "Promotor" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC__FoldChange"
            )) < 0 &
              !!parse_expr(Delta_col) > 0 &
              UCSC_RefGene_Group_short == "Promotor" ~ TRUE
          )
        ) 
    } else if(EpiTable %in% Splitted) {
      Proteins_with_Transcripts_TP <- Proteins_with_Transcripts_TP %>%
        mutate(
          Transkr_Epig__Correct_direction = case_when(IlmnID %in% Correct_Sites ~ TRUE),
          Protein_Epig__Correct_direction = case_when(
            !!parse_expr(paste0(
              TP,
              "_FC___FoldChange"
            )) < 0 &
              !!parse_expr(Delta_col) > 0 &
              UCSC_RefGene_Group_short == "Body" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC___FoldChange"
            )) > 0 &
              !!parse_expr(Delta_col) < 0 &
              UCSC_RefGene_Group_short == "Body" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC___FoldChange"
            )) < 0 &
              !!parse_expr(Delta_col) < 0 &
              UCSC_RefGene_Group_short == "Promotor" ~ TRUE,
            !!parse_expr(paste0(
              TP,
              "_FC___FoldChange"
            )) > 0 &
              !!parse_expr(Delta_col) > 0 &
              UCSC_RefGene_Group_short == "Promotor" ~ TRUE
          )
        )
    }
    
    
    Proteins_with_Transcripts_TP <- Proteins_with_Transcripts_TP %>% 
      rename(Correct__direction__of__Epigenomics__and__Transcriptomics = Transkr_Epig__Correct_direction,
             Correct__direction__of__Epigenomics__and__Proteomics = Protein_Epig__Correct_direction)
    
    if (grepl("001", EpiTable)) {
      Proteins_with_Transcripts_TP_list[[paste0(TP, "_001", ifelse(EpiTable %in% Unsplitted,
                                                                   "", "_Split_RES_LRE"))]] <- Proteins_with_Transcripts_TP
    } else {
      Proteins_with_Transcripts_TP_list[[paste0(TP, "_005", ifelse(EpiTable %in% Unsplitted,
                                                                   "", "_Split_RES_LRE"))]] <- Proteins_with_Transcripts_TP
    }
    
    print(paste0(format(Sys.time(), "%x %T"), " Finished: ", EpiTable, " - ", TP))
  }
  
}

## Notes and comments for xlsx ----
# Create 

Comment_1 <- c(X1 = paste("Correct direction of Epigenomics and Transcriptomics"),
               X2 = paste("Die Analyse zwischen Epigenomics und Transcriptomics wurde von Markus J\U00E4hnert gemacht"),
               X3 = NA_character_)

Comment_2 <- c(X1 = paste("Correct direction of Epigenomics and Proteomics"),
               X2 = paste("Die Analyse zwischen Epigenomics und Proteomics wurde von Thomas Goj gemacht"),
               X3 = NA_character_)

Comment_3 <- c(X1 = paste("Correct direction Regeln for unsplitted"),
               X2 = paste("Promotor: Prot./Transkr FC positiv + Epigenomics delta negativ \n",
                          "Promotor: Prot./Transkr FC negativ + Epigenomics delta positiv \n",
                          "Body: Prot./Transkr FC positiv + Epigenomics delta positiv \n",
                          "Body: Prot./Transkr FC negativ + Epigenomics delta negativ"),
               X3 = NA_character_)

Comment_4 <- c(
  X1 = paste("Correct Direction Regeln for RES-LRE splitted"),
  X2 = paste(
    "Promotor: Prot./Transkr FC negativ + Epigenomics delta negativ\n",
    "Promotor: Prot./Transkr FC positiv + Epigenomics delta positiv\n",
    "Body: Prot./Transkr FC negativ + Epigenomics delta positiv\n",
    "Body: Prot./Transkr FC positiv + Epigenomics delta negativ"
  ),
  X3 = NA_character_
)

Comments_List <- lapply(ls(pattern = "Comment_"), get)

Comments_File <- do.call(rbind, Comments_List)
names(Comments_File) <- NULL


## write xlsx ----

pValue_cols <- grep("pValue|rawP_|BH_", colnames(Proteins_with_Transcripts_TP))
Number_cols <- which(unlist(lapply(Proteins_with_Transcripts_TP, is.numeric)))
For_Border <- c(5)


Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")

wb <- createWorkbook()
# Sheets for comparisson
{
  for (sheet in names(Proteins_with_Transcripts_TP_list)) {
    
    colnames(Proteins_with_Transcripts_TP_list[[sheet]]) <- sapply(colnames(Proteins_with_Transcripts_TP_list[[sheet]]), 
                                                                   function (x) str_replace_all(x, "__", " "))
    addWorksheet(wb, sheetName = sheet)
    writeDataTable(wb, sheet = sheet, Proteins_with_Transcripts_TP_list[[sheet]])
    
    
    for (i in pValue_cols) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(Proteins_with_Transcripts_TP_list[[sheet]])+1,
                            rule = "<0.05", style = Cell_green,
                            type = "expression")
    }
    for (i in For_Border) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(Proteins_with_Transcripts_TP_list[[sheet]])+1,
                            rule = '!=""', style = Cell_border_right)
    }
    for (i in For_Border) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(Proteins_with_Transcripts_TP_list[[sheet]])+1,
                            rule = '=""', style = Cell_border_right)
    }
    for(i in Number_cols){
      addStyle(wb, sheet = sheet, style = Rounding_4, rows = 1:nrow(Proteins_with_Transcripts_TP_list[[sheet]])+1, cols = i)
    }
    
    freezePane(wb, sheet = sheet, firstCol = TRUE, firstRow = TRUE)
    setColWidths(wb, sheet = sheet, cols = 1:ncol(Proteins_with_Transcripts_TP_list[[sheet]]), widths = 16.5)
    addStyle(wb, sheet = sheet, style = style_linebreak, rows = 1, cols = 1:ncol(Proteins_with_Transcripts_TP_list[[sheet]]), 
             gridExpand = TRUE)
  }
}
# Sheet Comments
{
  addWorksheet(wb, "Comments")
  writeData(wb, sheet = "Comments", Comments_File, colNames = FALSE)
  setColWidths(wb, sheet = "Comments", cols = 1, widths = 25)
  setColWidths(wb, sheet = "Comments", cols = 2:3, widths = 100)
  addStyle(wb, sheet = "Comments", style = style_linebreak, rows = 1:nrow(Comments_File),
           cols = 0:ncol(Comments_File)+1, gridExpand = TRUE)
}
saveWorkbook(wb, file = "18_Overlap_ProteomicsNRE2_mit_EpigenetikData_Transcriptomic.xlsx", overwrite = TRUE)
