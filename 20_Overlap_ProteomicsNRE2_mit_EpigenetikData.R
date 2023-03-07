library(openxlsx)
library(tidyverse)
library(rlang)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")

# Limma Protein
load("01_ProteomicsNRE2_Limma_Analyse.RData")
Limma_Protein <- final_data_complete
rm(final_data_complete)

# Limma Epigeneitc
load("E:/Auswertung mit R/Epigentik Analysen NRE2/DataAnalysis/Epigenomics_Data_LimmaResults.RData")

## Pretesting how many Proteins are found in Epigenetics Data ----

# All timepoints have the same amount of CpG sites
# We filter all row which are no annotated a a specific gene
EpigenomicsData_Test = Epigenomics_Data_results[["T2_T0"]] %>% 
  filter(!UCSC_RefGene_Name_short == "")

Proteom_in_Epigenom = table(Limma_Protein$Gene_name %in% unique(EpigenomicsData_Test$UCSC_RefGene_Name_short))

Proteins_found = Limma_Protein$Gene_name[Limma_Protein$Gene_name %in% unique(EpigenomicsData_Test$UCSC_RefGene_Name_short)]

## Notes and comments for xlsx ----
# Create 
Comment_1 <- c(X1 = paste("Epigenetik Daten"),
               X2 = paste("In den Epigenetik Daten gibt es", nrow(EpigenomicsData_Test),
                          "zugeordnete CpG Sites. CpG Stellen ohne Gensymbol wurden gefiltert.",
                          "Zu diesen CpG Stellen finden wir", length(unique(EpigenomicsData_Test$UCSC_RefGene_Name_short)),
                          "Gene."),
               X3 = NA_character_)

Comment_2 <- c(X1 = paste("Proteomics Daten"),
               X2 = paste("In den Proteomics Daten gibt es", nrow(Limma_Protein),
                          "gemessene Proteine"),
               X3 = NA_character_)

Comment_3 <- c(X1 = paste("Overlap Epeigenetik und Proteomics"),
               X2 = paste("Wir finden von den Proteomics Daten insgesamt",
                          Proteom_in_Epigenom["TRUE"], "Proteine in den Epigenetik Daten",
                          "wieder und", Proteom_in_Epigenom["FALSE"], "nicht. Die Daten wurden",
                          "Anhand des Gensymbols gematcht."),
               X3 = NA_character_)

Comments_List <- lapply(ls(pattern = "Comment_"), get)

Comments_File <- do.call(rbind, Comments_List)
names(Comments_File) <- NULL

## Remove Test_Datatable so save space ----
rm(EpigenomicsData_Test, Proteom_in_Epigenom)


## Perform overlap of Proteom and epigenom ----
OverlapTable_raw = data.frame("Gene_Symbol" = Proteins_found)

OverlapTable_list = list()
OverlapTable_list_RData = list()


## Comparission of 
for(Timepoint in names(Epigenomics_Data_results)){
  
  TP_Proteom = str_remove(Timepoint, "_")
  # ProteomicsData = ProteomicsDataNRE2[[paste0("ProteomicsData__", TP_Proteom)]]
  
  EpigenomicsData = Epigenomics_Data_results[[Timepoint]] %>% filter(!UCSC_RefGene_Name_short == "")
  
  ## Select all cols needed for combination 
  
  # Epigenomics Data
  EpigenomicsData_short = EpigenomicsData %>%
    select(
      starts_with(c("mean", "sd", "delta")),
      contains(c("FC_", "rawP", "BH")),
      IlmnID,
      UCSC_RefGene_Group_short,
      UCSC_RefGene_Name_short
    ) %>%
    rename(Gene_Symbol = UCSC_RefGene_Name_short)
  
  # Proteomics Data
  Limma_Protein_short = Limma_Protein %>%
    select(
      Gene_name,
      Uniprot_ID,
      Description,
      Unique_Peptides,
      Valid_values,
      starts_with(paste0(TP_Proteom, "_FC__FoldChange"))
    ) %>%
    rename(Gene_Symbol = Gene_name)
  
  # Combine all Tables
  OverlapTable = 
    left_join(OverlapTable_raw,
              Limma_Protein_short,
              by = "Gene_Symbol") %>% 
    left_join(.,
              EpigenomicsData_short,
              by = "Gene_Symbol")
    
  
  ## Add Comments to Overlap Table (Significance, same Direction)
  
  # cols from Proteomics
  pValue_col = grep("__pValue", names(OverlapTable), value = TRUE)
  FC_col = grep("FoldChange$", names(OverlapTable), value = TRUE)
  
  # cols from Epigenomics
  delta_col = grep("delta", names(OverlapTable), value = TRUE)
  rawP_col = grep("rawP", names(OverlapTable), value = TRUE)
  
  # pValue < 0.05
  OverlapTable <- OverlapTable %>%
    mutate(Both__Significant__p005 = case_when(!!parse_expr(pValue_col) < 0.05 &
                                                !!parse_expr(rawP_col) < 0.05 ~ TRUE),
           Same__Direction__p005 = case_when(Both__Significant__p005 == TRUE &
                                              !!parse_expr(FC_col) > 0 &
                                              !!parse_expr(delta_col) > 0 &
                                              UCSC_RefGene_Group_short == "Body" ~ TRUE,
                                            Both__Significant__p005 == TRUE &
                                              !!parse_expr(FC_col) < 0 &
                                              !!parse_expr(delta_col) < 0 &
                                              UCSC_RefGene_Group_short == "Body" ~ TRUE,
                                            Both__Significant__p005 == TRUE &
                                              !!parse_expr(FC_col) > 0 &
                                              !!parse_expr(delta_col) < 0 &
                                              UCSC_RefGene_Group_short == "Promotor" ~ TRUE,
                                            Both__Significant__p005 == TRUE &
                                              !!parse_expr(FC_col) < 0 &
                                              !!parse_expr(delta_col) > 0 &
                                              UCSC_RefGene_Group_short == "Promotor" ~ TRUE)
    )
  
  # pValue < 0.05
  OverlapTable <- OverlapTable %>%
    mutate(Both__Significant__p001 = case_when(!!parse_expr(pValue_col) < 0.01 &
                                                !!parse_expr(rawP_col) < 0.01 ~ TRUE),
           Same__Direction__p001 = case_when(Both__Significant__p001 == TRUE &
                                              !!parse_expr(FC_col) > 0 &
                                              !!parse_expr(delta_col) > 0 &
                                              UCSC_RefGene_Group_short == "Body" ~ TRUE,
                                            Both__Significant__p001 == TRUE &
                                              !!parse_expr(FC_col) < 0 &
                                              !!parse_expr(delta_col) < 0 &
                                              UCSC_RefGene_Group_short == "Body" ~ TRUE,
                                            Both__Significant__p001 == TRUE &
                                              !!parse_expr(FC_col) > 0 &
                                              !!parse_expr(delta_col) < 0 &
                                              UCSC_RefGene_Group_short == "Promotor" ~ TRUE,
                                            Both__Significant__p001 == TRUE &
                                              !!parse_expr(FC_col) < 0 &
                                              !!parse_expr(delta_col) > 0 &
                                              UCSC_RefGene_Group_short == "Promotor" ~ TRUE)
    )
  
  
  OverlapTable_list_RData[[Timepoint]] <- OverlapTable
  
  SigGenes <- colSums(OverlapTable[grep("pValue|rawP|BH", names(OverlapTable))] < 0.05)
  names(OverlapTable)[grep("pValue|rawP|BH", names(OverlapTable))] <- 
    paste0(names(OverlapTable)[grep("pValue|rawP|BH", names(OverlapTable))], "__(", SigGenes, ")")
  
  TRUE_Cols <- colSums(OverlapTable[grep("__p005|__p001", names(OverlapTable))], na.rm = TRUE)
  names(OverlapTable)[grep("__p005|__p001", names(OverlapTable))] <- 
    paste0(names(OverlapTable)[grep("__p005|__p001", names(OverlapTable))], "__(", TRUE_Cols, ")")
  
    
  OverlapTable_list[[TP_Proteom]] <- OverlapTable
}

## Save RData List ----
save(OverlapTable_list_RData, file = "20_Overlap_ProteomicsNRE2_mit_EpigenetikData.RData")

## write xlsx ----

pValue_cols <- grep("pValue|rawP|BH", colnames(OverlapTable))
Number_cols <- which(unlist(lapply(OverlapTable, is.numeric)))
For_Border <- c(5, 8, 19)

Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")
Rounding_scientific = createStyle(numFmt = "0.00E+00")




wb <- createWorkbook()
# Sheets for comparisson
{
  for (sheet in names(OverlapTable_list)) {
    
    colnames(OverlapTable_list[[sheet]]) <- sapply(colnames(OverlapTable_list[[sheet]]), 
                                                                   function (x) str_replace_all(x, "__", " "))
    addWorksheet(wb, sheetName = sheet)
    writeDataTable(wb, sheet = sheet, OverlapTable_list[[sheet]])
    
    
    for (i in pValue_cols) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(OverlapTable_list[[sheet]])+1,
                            rule = "<0.05", style = Cell_green,
                            type = "expression")
    }
    for (i in For_Border) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(OverlapTable_list[[sheet]])+1,
                            rule = '!=""', style = Cell_border_right)
    }
    for (i in For_Border) {
      conditionalFormatting(wb, sheet = sheet, cols = i, 
                            rows = 1:nrow(OverlapTable_list[[sheet]])+1,
                            rule = '=""', style = Cell_border_right)
    }
    for(i in Number_cols){
      addStyle(wb, sheet = sheet, style = Rounding_4, rows = 1:nrow(OverlapTable_list[[sheet]])+1, cols = i)
    }
    
    
    ## Selecting all pValue-cells which have values <0.0001
    Below_0001 <- lapply(OverlapTable_list[[sheet]][, pValue_cols], function (x) x < 0.0001)
    names(Below_0001) <- which(colnames(OverlapTable_list[[sheet]]) %in% names(Below_0001))
    for(i in names(Below_0001)){
      
      if (all(Below_0001[[i]] == FALSE)) {
        Below_0001[[i]] <- NULL
      } else {
        Below_0001[[i]] <- which(row.names(OverlapTable) %in% row.names(OverlapTable)[Below_0001[[i]]])
      }
    }
    for (col in names(Below_0001)) {
      for (i in Below_0001[[col]]) {
        addStyle(wb, sheet = sheet, style = Rounding_scientific, rows = i+1, cols = col)
      }
    }
    
    freezePane(wb, sheet = sheet, firstCol = TRUE, firstRow = TRUE)
    setColWidths(wb, sheet = sheet, cols = 1:ncol(OverlapTable_list[[sheet]]), widths = 16.5)
    addStyle(wb, sheet = sheet, style = style_linebreak, rows = 1, cols = 1:ncol(OverlapTable_list[[sheet]]), 
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
saveWorkbook(wb, file = "20_Overlap_ProteomicsNRE2_mit_EpigenetikData.xlsx", overwrite = TRUE)
