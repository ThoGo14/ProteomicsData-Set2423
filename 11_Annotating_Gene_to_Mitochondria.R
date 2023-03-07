library(tidyverse)
library(openxlsx)
library(rlang)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Read Data ----

load("ProteomicsNRE2__Expression.RData")

load("05_Korrelation_ProteomicsNRE2_mit_KlinischenDaten.RData")

#### Annotated by GO-Terms
Mitochondrion_ChildTerms <- read.delim("E:/Auswertung mit R/00_Gene_Annotation/GeneOntology/ChildTerms__Mitochondrion.txt")

AerobicRespiration_ChildTerms <- read.delim("E:/Auswertung mit R/00_Gene_Annotation/GeneOntology/ChildTerms__AerobicRespiration.txt")

Proteomics_Protein_Info_GO <- ProteomicsDataNRE2[["Proteomics_Protein_Info_GO"]]


unique(Proteomics_Protein_Info_GO$GO.Domain)

## Select Mitochondrion (GO:0005739) -----

## It is not needed to select first the Genes which belong to mitochondrion and then
## mitochondrion + child terms and filter for those genes. Some genes might belong
## to a child term, but are somehow not annotated to mitochondrion itself

# ProteomicsDataNRE2_Mitochondrion_GO <- Proteomics_Protein_Info_GO %>% 
#   filter(GO.Domain == "cellular_component", GO.ID == "GO:0005739")
# 
# Mitochondrion_Genes <- unique(ProteomicsDataNRE2_Mitochondrion_GO$Gene_name)

ProteomicsDataNRE2_MitoSelected <- Proteomics_Protein_Info_GO %>% 
  filter(GO.ID %in% c(Mitochondrion_ChildTerms$GO.ID, "GO:0005739")) %>% 
  mutate(Unique_Name = paste(Gene_name, GO.ID, sep = "__"))

for (Unique_Name in unique(ProteomicsDataNRE2_MitoSelected$Unique_Name)) {
  Evidence_codes <- grep(Unique_Name, ProteomicsDataNRE2_MitoSelected$Unique_Name)
  ProteomicsDataNRE2_MitoSelected[Evidence_codes, "GO.Evidence_code"] <- paste0(ProteomicsDataNRE2_MitoSelected$GO.Evidence_code[Evidence_codes],
                                                                                collapse = " | ")
}


## aerobic respiration (GO:0009060) -----

## It is not needed to select first the Genes which belong to aerobicResp and then
## aerobicResp + child terms and filter for those genes. Some genes might belong
## to a child term, but are somehow not annotated to aerobicResp itself

# ProteomicsDataNRE2_aerobicResp_GO <- Proteomics_Protein_Info_GO %>% 
#   filter(GO.Domain == "biological_process", GO.ID == "GO:0009060")
# 
# aerobicRespiration_Genes <- unique(ProteomicsDataNRE2_aerobicResp_GO$Gene_name)

ProteomicsDataNRE2_aeroRespSelected <- Proteomics_Protein_Info_GO %>% 
  filter(GO.ID %in% c(AerobicRespiration_ChildTerms$GO.ID, "GO:0009060")) %>% 
  mutate(Unique_Name = paste(Gene_name, GO.ID, sep = "__"))

for (Unique_Name in unique(ProteomicsDataNRE2_aeroRespSelected$Unique_Name)) {
  Evidence_codes <- grep(Unique_Name, ProteomicsDataNRE2_aeroRespSelected$Unique_Name)
  ProteomicsDataNRE2_aeroRespSelected[Evidence_codes, "GO.Evidence_code"] <- paste0(ProteomicsDataNRE2_aeroRespSelected$GO.Evidence_code[Evidence_codes],
                                                                                    collapse = " | ")
}

# ProteomicsDataNRE2_aeroRespSelected <- distinct(ProteomicsDataNRE2_aeroRespSelected, Unique_Name, .keep_all = TRUE)

# ProteomicsDataNRE2_aeroRespSelected$Unique_Name <- NULL

## Conmbine Mitochondrion and 

ProteomicsDataNRE2_MitoSelected <- rbind(ProteomicsDataNRE2_MitoSelected,
                                         ProteomicsDataNRE2_aeroRespSelected)

ProteomicsDataNRE2_MitoSelected <- distinct(ProteomicsDataNRE2_MitoSelected, Unique_Name, .keep_all = TRUE)
 
ProteomicsDataNRE2_MitoSelected$Unique_Name <- NULL
ProteomicsDataNRE2_MitoSelected[, 1] <- NULL

## Add Limma Analysis Significant Proteins ----
ProteomicsDataNRE2_MitoSelected[, c("LimmaAnalysis__Significant")] <- NA

load("01_ProteomicsNRE2_Limma_Analyse.RData")
LimmaAnalysis <- final_data_complete
rm(final_data_complete)

LimmaAnalysis <- LimmaAnalysis[, c("Uniprot_ID", grep("__FoldChange", names(LimmaAnalysis), value = TRUE))]
LimmaTP <- sapply(names(LimmaAnalysis[-1]), function (x) unlist(str_split(x, "_FC__"))[1]) %>% unique()

for (TP in LimmaTP) {
  LimmaAnalysis_short <- LimmaAnalysis %>% select(Uniprot_ID, starts_with(TP)) %>% 
    filter(!!parse_expr(paste0(TP, "_FC__FoldChange__pValue")) < 0.05)
  
  ProteomicsDataNRE2_MitoSelected <- left_join(ProteomicsDataNRE2_MitoSelected,
                                               LimmaAnalysis_short,
                                               by = "Uniprot_ID")
}

## Add Respiration Significant Proteins ----
ProteomicsDataNRE2_MitoSelected[, c("Respiration__Significant")] <- NA
# for(TP in names(KorrelationTable_Results_List)){
#   ProteomicsDataNRE2_MitoSelected[, paste0(TP, "__pValue")] <- NA
#   ProteomicsDataNRE2_MitoSelected[, paste0(TP, "__pearsonR")] <- NA
# }

names(KorrelationTable_Results_List) %in% LimmaTP

for(TP in LimmaTP){
  KorrelationTable_Results <- KorrelationTable_Results_List[[TP]]
  
  ## We need to unquote the paste0() expression and then use the special assign operator :=
  SignifGenes <- KorrelationTable_Results %>% 
    filter(SM_MODP.S_FC__pValue < 0.05) %>% 
    select(Gene_name, SM_MODP.S_FC__pValue, SM_MODP.S_FC__pearsonR) %>% 
    rename(!!paste0(TP, "__pValue") := SM_MODP.S_FC__pValue,
           !!paste0(TP, "__pearsonR") := SM_MODP.S_FC__pearsonR) 
  
  ProteomicsDataNRE2_MitoSelected <- left_join(ProteomicsDataNRE2_MitoSelected,
                                               SignifGenes,
                                               by = "Gene_name")
}



## Significant Limma and Reespiration -----
ProteomicsDataNRE2_MitoSelected[, c("Comparission__LimmaAnalysis__and__Respiration")] <- NA


for (TP in LimmaTP) {
  ProteomicsDataNRE2_MitoSelected <- ProteomicsDataNRE2_MitoSelected %>% 
    mutate(!!parse_expr(paste0(TP, "__LimmaAnalysis__and__Respiration__significant")) := 
             case_when(
               !is.na(!!parse_expr(paste0(TP, "_FC__FoldChange__pValue"))) &
                 !is.na(!!parse_expr(paste0(TP, "__pValue"))) ~ TRUE),
           !!parse_expr(paste0(TP, "__Direction__LimmaAnalysis____Respiration")) :=
             case_when(!is.na(!!parse_expr(paste0(TP, "__LimmaAnalysis__and__Respiration__significant"))) &
                         !!parse_expr(paste0(TP, "_FC__FoldChange")) > 0 &
                         !!parse_expr(paste0(TP, "__pearsonR")) > 0 ~ "up | up",
                       !is.na(!!parse_expr(paste0(TP, "__LimmaAnalysis__and__Respiration__significant"))) &
                         !!parse_expr(paste0(TP, "_FC__FoldChange")) < 0 &
                         !!parse_expr(paste0(TP, "__pearsonR")) < 0 ~ "down | down",
                       !is.na(!!parse_expr(paste0(TP, "__LimmaAnalysis__and__Respiration__significant"))) &
                         !!parse_expr(paste0(TP, "_FC__FoldChange")) < 0 &
                         !!parse_expr(paste0(TP, "__pearsonR")) > 0 ~ "down | up",
                       !is.na(!!parse_expr(paste0(TP, "__LimmaAnalysis__and__Respiration__significant"))) &
                         !!parse_expr(paste0(TP, "_FC__FoldChange")) > 0 &
                         !!parse_expr(paste0(TP, "__pearsonR")) < 0 ~ "up | down")
           )
  
}
names(ProteomicsDataNRE2_MitoSelected) <- sapply(names(ProteomicsDataNRE2_MitoSelected), function (x) str_replace_all(x,
                                                                                                                      "____",
                                                                                                                      "__|__"))


## Write xlsx ----
## Table ProteomicsDataNRE2_MitoSelected
numeric_cols <- which(unlist(lapply(ProteomicsDataNRE2_MitoSelected, is.numeric)))
pValue_cols <- grep("__pValue", colnames(ProteomicsDataNRE2_MitoSelected))

For_Border <- c(grep("GO.Evidence_code", names(ProteomicsDataNRE2_MitoSelected)),
                grep("__adj.pValue", names(ProteomicsDataNRE2_MitoSelected)),
                grep("__pearsonR", names(ProteomicsDataNRE2_MitoSelected)),
                grep("__|__", names(ProteomicsDataNRE2_MitoSelected), fixed = TRUE))


Cells_for_Orange <- c(grep("LimmaAnalysis__Significant", names(ProteomicsDataNRE2_MitoSelected)),
                      grep("Respiration__Significant", names(ProteomicsDataNRE2_MitoSelected)),
                      grep("Comparission__LimmaAnalysis__and__Respiration", names(ProteomicsDataNRE2_MitoSelected)))

SigGenes <- colSums(ProteomicsDataNRE2_MitoSelected[, pValue_cols] < 0.05, na.rm = T)

colnames(ProteomicsDataNRE2_MitoSelected)[pValue_cols] <- paste0(colnames(ProteomicsDataNRE2_MitoSelected)[pValue_cols], 
                                                          "__(", SigGenes, ")")

colnames(ProteomicsDataNRE2_MitoSelected) <- sapply(colnames(ProteomicsDataNRE2_MitoSelected), function (x){
  x <- str_replace_all(x, "__", " ")
})


# Keine unter 0.0001 gefunden
# Below_0001 <- lapply(ProteomicsDataNRE2_MitoSelected[, pValue_cols], function (x) x < 0.0001)
# names(Below_0001) <- which(colnames(ProteomicsDataNRE2_MitoSelected) %in% names(Below_0001))
# for(i in names(Below_0001)){
#   
#   if (all(Below_0001[[i]] == FALSE | is.na(Below_0001[[i]]))) {
#     Below_0001[[i]] <- NULL
#   } else {
#     Below_0001[[i]] <- which(row.names(ProteomicsDataNRE2_MitoSelected) %in% row.names(ProteomicsDataNRE2_MitoSelected)[Below_0001[[i]]])
#   }
# }


wb <- createWorkbook()
addWorksheet(wb, "MitoSelected_GO-Terms")
writeDataTable(wb, sheet = 1, ProteomicsDataNRE2_MitoSelected)


Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#FF962D") ## medium orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")
Rounding_scientific = createStyle(numFmt = "0.00E+00")

## Sheet 1
{
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = 1, cols = i, 
                          rows = 1:nrow(ProteomicsDataNRE2_MitoSelected)+1,
                          rule = '!=""', style = Cell_border_right)
  }
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = 1, cols = i, 
                          rows = 1:nrow(ProteomicsDataNRE2_MitoSelected)+1,
                          rule = '=""', style = Cell_border_right)
  }
  freezePane(wb, sheet = 1, firstActiveCol = 3, firstActiveRow = 2)
  setColWidths(wb, sheet = 1, cols = 1:ncol(ProteomicsDataNRE2_MitoSelected), widths = 16.5)
  addStyle(wb, sheet = 1, style = style_linebreak, rows = 1, cols = 1:ncol(ProteomicsDataNRE2_MitoSelected), gridExpand = TRUE)
  for(i in Cells_for_Orange) {
    conditionalFormatting(wb, sheet = 1, cols = i, 
                          rows = 1,
                          rule = '!=""', style = Cell_orange)
  }
  
  for (i in numeric_cols){
    addStyle(wb, sheet = 1, style = Rounding_4, rows = 1:nrow(ProteomicsDataNRE2_MitoSelected)+1, 
             cols = i, gridExpand = TRUE)
  }
  # for (col in names(Below_0001)) {
  #   for (i in Below_0001[[col]]) {
  #     addStyle(wb, sheet = 1, style = Rounding_scientific, rows = i+1, cols = col)
  #   }
  # }
}

## Add Sheets for evidence codes and explanation
{
  Evidence_WB <- loadWorkbook("E:/Auswertung mit R/00_Gene_Annotation/GeneOntology/00_Evidence_Code.xlsx")
  lapply(names(Evidence_WB), function(s) {
    dt <- read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/GeneOntology/00_Evidence_Code.xlsx", sheet = s)
    addWorksheet(wb , sheetName = s)
    writeDataTable(wb, s, dt)
  })
  setColWidths(wb, sheet = "Evidence_codes", cols = 1:3, widths = "auto")
  setColWidths(wb, sheet = "Evidence_Explanation", cols = 1, widths = "auto")
  setColWidths(wb, sheet = "Evidence_Explanation", cols = 2, widths = 120)
  addStyle(wb, sheet = "Evidence_Explanation", style = style_linebreak, rows = 1:28, 
           cols = 1:2, gridExpand = TRUE)
  
}
saveWorkbook(wb, file = "11_Annotating_Gene_to_Mitochondria.xlsx", overwrite = TRUE)

