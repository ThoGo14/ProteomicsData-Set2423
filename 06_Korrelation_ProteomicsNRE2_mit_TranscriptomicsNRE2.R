library(openxlsx)
library(tidyverse)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")
load("E:/Auswertung mit R/00_Helpful Files/NRE2_All_TimePoints_Transcripts_FC_no_DABG_remove.RData")

NRE2_GenesTable <- readRDS("E:/Auswertung mit R/00_Helpful Files/NRE2_GenesTable_no_DABG_remove.RDS")

Proteomics_Proteins <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]$Gene_name

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/20210315_Klinische_Daten.xlsx")

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

## functions ----
Return_ProteinTable <- function(ProtTable){
  ProtTable_names = row.names(ProtTable)
  ProtTable = t(ProtTable)
  ProtTable = data.frame(MyoID = row.names(ProtTable), ProtTable)
  ProtTable$MyoID = sapply(ProtTable$MyoID, function(x) {
    x = unlist(str_split(x, "__"))[1]
    return(x)
  })
  colnames(ProtTable) = c("MyoID", ProtTable_names)
  return(ProtTable)
}

Return_ProteinTable_BL <- function(ProtTable, TimePoint){
  ProtTable_names = row.names(ProtTable)
  ProtTable = ProtTable[, grep(TimePoint, colnames(ProtTable))]
  ProtTable = t(ProtTable)
  ProtTable = data.frame(MyoID = row.names(ProtTable), ProtTable)
  ProtTable$MyoID = sapply(ProtTable$MyoID, function(x) {
    x = unlist(str_split(x, "_"))[1]
    return(x)
  })
  colnames(ProtTable) = c("MyoID", ProtTable_names)
  return(ProtTable)
}
## Create Table - FC ----
Welche_TP <- c("T2T0", "T4T5", "T5T0", "T4T2")

Proteins_with_Transcripts <- rbind(Proteins_notFound,
                                   cbind(Protein.Name = Proteins_found, 
                                         Gene.Symbol = Proteins_found))
row.names(Proteins_with_Transcripts) <- NULL

KorrelationTable_Results <- NRE2_GenesTable[grep(paste0("^", Proteins_with_Transcripts$Gene.Symbol, "$",
                                                        collapse = "|"),
                                                 NRE2_GenesTable$Gene.Symbol), 
                                            ] %>% 
  select(-Cytoband, -DABG, -Probeset) %>% 
  left_join(., Proteins_with_Transcripts,
            by = "Gene.Symbol") %>% 
  left_join(., 
            ProteomicsDataNRE2[["Proteomics_Protein_Info"]] %>% select(Gene_name, Uniprot_ID) %>% 
              rename(Protein.Name = Gene_name),
            by = "Protein.Name") %>% 
  mutate(Unique = str_replace_all(Unique, "-", "\\."),
         Gen_Protein = paste(Unique, Uniprot_ID, sep = "___")) %>% 
  filter(!duplicated(Gen_Protein))
  



for (TP in Welche_TP){
  KorrelationTable_Results[, paste0(TP, "__pValue")] <- NA
  KorrelationTable_Results[, paste0(TP, "__adj.pValue__Sex+Age")] <- NA
  KorrelationTable_Results[, paste0(TP, "__adj.pValue__Sex+Age+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(TP, "__adj.pValue__Sex+Age+BMI")] <- NA
  KorrelationTable_Results[, paste0(TP, "__adj.pValue__Sex+Age+BMI+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(TP, "__pearsonR")] <- NA
}

KorrelationTable_Results[, "Limma Analysis"] <- NA

## Add Limma Protein
Limma_Protein <- Limma_Protein[, c(2, (which(colnames(Limma_Protein) == "linearGroupFC")+1L):
                                      (which(colnames(Limma_Protein) == "linearGroupRatios")-1L))] %>% 
  select(Gene_name, starts_with(Welche_TP)) %>% 
  rename(Protein.Name = Gene_name) %>% 
  rename_if(is.numeric, ~paste0("Limma__Protein__", .))


KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      Limma_Protein,
                                      by = "Protein.Name")

## Add Limma Transcripts
Limma_Transcripts <- Limma_Transcripts %>% select(Gene.Symbol, ends_with("__All")) %>% 
  rename_if(is.numeric, ~paste0("Limma__Transcripts__", .))

KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      Limma_Transcripts,
                                      by = "Gene.Symbol")


KorrelationTable_Results <- KorrelationTable_Results %>% distinct(Gen_Protein, .keep_all = TRUE)

for (TP in Welche_TP){
  
  if (TP == "T2T0") {
    ProteinTable <- Return_ProteinTable(ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]])
    
    Table_For_Korrelation <- left_join(NRE2_BL.ACvsBL.RE, 
                                       ProteinTable, 
                                       by = "MyoID")
    
    Table_For_Korrelation <- left_join(Table_For_Korrelation, KlinischeDaten, by = "MyoID")
    
    
  } else if (TP == "T4T5") {
    ProteinTable <- Return_ProteinTable(ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]])
    
    Table_For_Korrelation <- left_join(NRE2_TR.ACvsTR.RE, 
                                       ProteinTable, 
                                       by = "MyoID")
    
    Table_For_Korrelation <- left_join(Table_For_Korrelation, KlinischeDaten, by = "MyoID")
    
    
  } else if (TP == "T5T0") {
    ProteinTable <- Return_ProteinTable(ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]])
    
    Table_For_Korrelation <- left_join(NRE2_TR.REvsBL.RE, 
                                       ProteinTable, 
                                       by = "MyoID")
    
    Table_For_Korrelation <- left_join(Table_For_Korrelation, KlinischeDaten, by = "MyoID")
    
    
  } else if (TP == "T4T2") {
    ProteinTable <- Return_ProteinTable(ProteomicsDataNRE2[[paste0("ProteomicsData__", TP)]])
    
    Table_For_Korrelation <- left_join(NRE2_TR.ACvsBL.AC, 
                                       ProteinTable, 
                                       by = "MyoID")
    
    Table_For_Korrelation <- left_join(Table_For_Korrelation, KlinischeDaten, by = "MyoID")
    
  }
  
  
  
  ColNR_pValue <- grep(paste0(TP, "__pValue"), colnames(KorrelationTable_Results))
  ColNR_pearsonR <- grep(paste0(TP, "__pearsonR"), colnames(KorrelationTable_Results))
  ColNR_adj.pValue1 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results))[1]
  ColNR_adj.pValue2 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results))[2]
  ColNR_adj.pValue3 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results))[3]
  ColNR_adj.pValue4 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results))[4]
  
  for(Gen_unique in KorrelationTable_Results$Gen_Protein){
    RowNR <- grep(paste0("^", Gen_unique, "$"), KorrelationTable_Results$Gen_Protein)
    Gen <- unlist(str_split(Gen_unique, "___"))[1]
    Protein <- KorrelationTable_Results[RowNR, "Uniprot_ID"]
    
    Correlation <- cor.test(Table_For_Korrelation[, Gen], 
                            Table_For_Korrelation[, Protein])
    
    LM1 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0))
    
    LM2 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$ISIMATS_BL))
    
    LM3 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0))
    
    LM4 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0 +
                        Table_For_Korrelation$ISIMATS_BL))
    
    
    
    KorrelationTable_Results[RowNR, ColNR_pValue] <- Correlation[["p.value"]]
    KorrelationTable_Results[RowNR, ColNR_pearsonR] <- Correlation[["estimate"]][["cor"]]
    KorrelationTable_Results[RowNR, ColNR_adj.pValue1] <- LM1[["coefficients"]][2, 4]
    KorrelationTable_Results[RowNR, ColNR_adj.pValue2] <- LM2[["coefficients"]][2, 4]
    KorrelationTable_Results[RowNR, ColNR_adj.pValue3] <- LM3[["coefficients"]][2, 4]
    KorrelationTable_Results[RowNR, ColNR_adj.pValue4] <- LM4[["coefficients"]][2, 4]
    
    
    ## Adding Limma
    
  }
  
  print(paste0(format(Sys.time(), "%x %T"), " Finished: ", TP))
  
  
}

KorrelationTable_Results$Gen_Protein <- NULL

## Create Table - Baseline ----
TP_Table <- data.frame(TP_old = c("1", "2", "3", "4"),
                       TP_new = c("pre_resting", "pre_acuteexercise",
                                  "post_acuteexercise", "post_resting"))

Proteins_with_Transcripts <- rbind(Proteins_notFound,
                                   cbind(Protein.Name = Proteins_found, 
                                         Gene.Symbol = Proteins_found))

KorrelationTable_Results_BL <- NRE2_GenesTable[grep(paste0("^", Proteins_with_Transcripts$Gene.Symbol, "$",
                                                           collapse = "|"),
                                                    NRE2_GenesTable$Gene.Symbol), ] %>% 
  left_join(., Proteins_with_Transcripts,
            by = "Gene.Symbol") %>% 
  left_join(., 
            ProteomicsDataNRE2[["Proteomics_Protein_Info"]] %>% select(Gene_name, Uniprot_ID) %>% 
              rename(Protein.Name = Gene_name),
            by = "Protein.Name") %>%
  mutate(Unique = str_replace_all(Unique, "-" , "\\."),
         Gen_Protein = paste(Unique, Uniprot_ID, sep = "___")) %>% 
  filter(!duplicated(Gen_Protein))



for (TP in TP_Table$TP_new){
  KorrelationTable_Results_BL[, paste0(TP, "__pValue")] <- NA
  KorrelationTable_Results_BL[, paste0(TP, "__adj.pValue__Sex+Age")] <- NA
  KorrelationTable_Results_BL[, paste0(TP, "__adj.pValue__Sex+Age+ISI_BL")] <- NA
  KorrelationTable_Results_BL[, paste0(TP, "__adj.pValue__Sex+Age+BMI")] <- NA
  KorrelationTable_Results_BL[, paste0(TP, "__adj.pValue__Sex+Age+BMI+ISI_BL")] <- NA
  KorrelationTable_Results_BL[, paste0(TP, "__pearsonR")] <- NA
  
}


for (TP in TP_Table$TP_new){
  
  
  ProteinTable <- Return_ProteinTable_BL(ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]], TP)
  
  Transcript_TP <- TP_Table[TP_Table$TP_new == TP, 1]
  TranscrTable <- NRE2_TimePoints_Baseintensity_noDABG %>% filter(TimePoint == Transcript_TP)
  
  Table_For_Korrelation <- left_join(TranscrTable, 
                                     ProteinTable, 
                                     by = "MyoID")
  
  Table_For_Korrelation <- left_join(Table_For_Korrelation, KlinischeDaten, by = "MyoID")
  
  
  
  
  ColNR_pValue <- grep(paste0(TP, "__pValue"), colnames(KorrelationTable_Results_BL))
  ColNR_pearsonR <- grep(paste0(TP, "__pearsonR"), colnames(KorrelationTable_Results_BL))
  ColNR_adj.pValue1 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results_BL))[1]
  ColNR_adj.pValue2 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results_BL))[2]
  ColNR_adj.pValue3 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results_BL))[3]
  ColNR_adj.pValue4 <- grep(paste0(TP, "__adj.pValue"), colnames(KorrelationTable_Results_BL))[4]
  
  for(Gen_unique in KorrelationTable_Results_BL$Gen_Protein){
    RowNR <- grep(paste0("^", Gen_unique, "$"), KorrelationTable_Results_BL$Gen_Protein)
    Gen <- unlist(str_split(Gen_unique, "___"))[1]
    Protein <- KorrelationTable_Results_BL[RowNR, "Uniprot_ID"]
    
    Correlation <- cor.test(Table_For_Korrelation[, Gen], 
                            Table_For_Korrelation[, Protein])
    
    LM1 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0))
    
    LM2 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$ISIMATS_BL))
    
    LM3 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0))
    
    LM4 <- summary(lm(Table_For_Korrelation[, Gen] ~ Table_For_Korrelation[, Protein] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0 +
                        Table_For_Korrelation$ISIMATS_BL))
    
    
    
    KorrelationTable_Results_BL[RowNR, ColNR_pValue] <- Correlation[["p.value"]]
    KorrelationTable_Results_BL[RowNR, ColNR_pearsonR] <- Correlation[["estimate"]][["cor"]]
    KorrelationTable_Results_BL[RowNR, ColNR_adj.pValue1] <- LM1[["coefficients"]][2, 4]
    KorrelationTable_Results_BL[RowNR, ColNR_adj.pValue2] <- LM2[["coefficients"]][2, 4]
    KorrelationTable_Results_BL[RowNR, ColNR_adj.pValue3] <- LM3[["coefficients"]][2, 4]
    KorrelationTable_Results_BL[RowNR, ColNR_adj.pValue4] <- LM4[["coefficients"]][2, 4]
    
  }
  
  print(paste0(format(Sys.time(), "%x %T"), " Finished: ", TP))
  
  
}

KorrelationTable_Results_BL$Gen_Protein <- NULL



## Save Data ----

save(KorrelationTable_Results, KorrelationTable_Results_BL, file = "06_Korrelation_ProteomicsNRE2_mit_TranscriptomicsNRE2.RData")
# load("06_Korrelation_ProteomicsNRE2_mit_TranscriptomicsNRE2.RData")

## Write xlsx ----

## FC Table
pValue_cols <- grep("pValue|P.Value", colnames(KorrelationTable_Results))
# For_Border <- seq(6, ncol(KorrelationTable_Results), by = 6)
For_Border <- c(6, 
                grep("__pearsonR", colnames(KorrelationTable_Results)), 
                grep("__FoldChange__adj.pValue", colnames(KorrelationTable_Results)),
                grep("__Adj.P.Value$", colnames(KorrelationTable_Results)),
                grep("__FC__All$", colnames(KorrelationTable_Results)))


## Selecting all pValue-cells which have values <0.0001
Below_0001 <- lapply(KorrelationTable_Results[, pValue_cols], function (x) x < 0.0001)
names(Below_0001) <- which(colnames(KorrelationTable_Results) %in% names(Below_0001))
for(i in names(Below_0001)){
  
  if (all(Below_0001[[i]] == FALSE)) {
    Below_0001[[i]] <- NULL
  } else {
    Below_0001[[i]] <- which(row.names(KorrelationTable_Results) %in% row.names(KorrelationTable_Results)[Below_0001[[i]]])
  }
}

SigGenes <- colSums(KorrelationTable_Results[, pValue_cols] < 0.05)

colnames(KorrelationTable_Results)[pValue_cols] <- paste0(colnames(KorrelationTable_Results)[pValue_cols], 
                                                          "__(", SigGenes, ")")

colnames(KorrelationTable_Results) <- sapply(colnames(KorrelationTable_Results), function (x){
  x <- str_remove(x, "__All")
  x <- str_replace_all(x, "__", " ")
})


## Baseline Table
pValue_cols2 <- grep("pValue", colnames(KorrelationTable_Results_BL))
For_Border2 <- seq(9, ncol(KorrelationTable_Results_BL), by = 6)


## Selecting all pValue-cells which have values <0.0001
Below_0001_BL <- lapply(KorrelationTable_Results_BL[, pValue_cols2], function (x) x < 0.0001)
names(Below_0001_BL) <- which(colnames(KorrelationTable_Results_BL) %in% names(Below_0001_BL))
for(i in names(Below_0001_BL)){
  
  if (all(Below_0001_BL[[i]] == FALSE)) {
    Below_0001_BL[[i]] <- NULL
  } else {
    Below_0001_BL[[i]] <- which(row.names(KorrelationTable_Results_BL) %in% row.names(KorrelationTable_Results_BL)[Below_0001_BL[[i]]])
  }
}

SigGenes <- colSums(KorrelationTable_Results_BL[, pValue_cols2] < 0.05)

colnames(KorrelationTable_Results_BL)[pValue_cols2] <- paste0(colnames(KorrelationTable_Results_BL)[pValue_cols2], 
                                                              "__(", SigGenes, ")")

colnames(KorrelationTable_Results_BL) <- sapply(colnames(KorrelationTable_Results_BL), function (x){
  x <- str_replace_all(x, "__", " ")
})


wb <- createWorkbook()
addWorksheet(wb, "Prot_Transcr_FC")
writeDataTable(wb, sheet = 1, KorrelationTable_Results)
addWorksheet(wb, "Prot_Transcr_BL")
writeDataTable(wb, sheet = 2, KorrelationTable_Results_BL)


Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#FF962D") ## medium orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")
Rounding_scientific = createStyle(numFmt = "0.00E+00")

## Sheet 1
for (i in pValue_cols) {
  conditionalFormatting(wb, sheet = 1, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results)+1,
                        rule = "<0.05", style = Cell_green,
                        type = "expression")
}
for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results)+1,
                        rule = '!=""', style = Cell_border_right)
}
for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results)+1,
                        rule = '=""', style = Cell_border_right)
}


freezePane(wb, sheet = 1, firstCol = TRUE, firstRow = TRUE)
setColWidths(wb, sheet = 1, cols = 1:ncol(KorrelationTable_Results), widths = 16.5)
addStyle(wb, sheet = 1, style = style_linebreak, rows = 1, cols = 1:ncol(KorrelationTable_Results), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = Rounding_4, rows = 1:nrow(KorrelationTable_Results)+1, 
         cols = 5:ncol(KorrelationTable_Results), gridExpand = TRUE)
for (col in names(Below_0001)) {
  for (i in Below_0001[[col]]) {
    addStyle(wb, sheet = 1, style = Rounding_scientific, rows = i+1, cols = col)
  }
}

## Sheet 2
for (i in pValue_cols2) {
  conditionalFormatting(wb, sheet = 2, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_BL)+1,
                        rule = "<0.05", style = Cell_green,
                        type = "expression")
}
for (i in For_Border2) {
  conditionalFormatting(wb, sheet = 2, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_BL)+1,
                        rule = '!=""', style = Cell_border_right)
}
for (i in For_Border2) {
  conditionalFormatting(wb, sheet = 2, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_BL)+1,
                        rule = '=""', style = Cell_border_right)
}


freezePane(wb, sheet = 2, firstCol = TRUE, firstRow = TRUE)
setColWidths(wb, sheet = 2, cols = 1:ncol(KorrelationTable_Results_BL), widths = 16.5)
addStyle(wb, sheet = 2, style = style_linebreak, rows = 1, cols = 1:ncol(KorrelationTable_Results_BL), gridExpand = TRUE)
addStyle(wb, sheet = 2, style = Rounding_4, rows = 1:nrow(KorrelationTable_Results_BL)+1, 
         cols = 5:ncol(KorrelationTable_Results_BL), gridExpand = TRUE)
for (col in names(Below_0001_BL)) {
  for (i in Below_0001_BL[[col]]) {
    addStyle(wb, sheet = 2, style = Rounding_scientific, rows = i+1, cols = col)
  }
}
saveWorkbook(wb, file = "06_Korrelation_ProteomicsNRE2_mit_TranscriptomicsNRE2.xlsx", overwrite = TRUE)






