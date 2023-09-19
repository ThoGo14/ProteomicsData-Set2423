library(openxlsx)
library(tidyverse)


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

# row.names(Proteomics_Protein_Info) <- Proteomics_Protein_Info$Gene.name

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/00_Klinische_Daten.xlsx")

#### Select Mitochondrion (GO:0005739) -----
Proteomics_Protein_Info_GO <- ProteomicsDataNRE2[["Proteomics_Protein_Info_GO"]]

ProteomicsDataNRE2_Mitochondrion_GO <- Proteomics_Protein_Info_GO %>% 
  filter(GO.Domain == "cellular_component", GO.ID == "GO:0005739") %>% 
  select(Gene.name, GO.Term) %>% 
  distinct(Gene.name, .keep_all = TRUE)
  
ProteomicsDataNRE2_Mitochondrion_GO <- rename(ProteomicsDataNRE2_Mitochondrion_GO,
                                              Protein = Gene.name,
                                              !!paste0("Mitochondrion (GO:0005739) (",
                                                     nrow(ProteomicsDataNRE2_Mitochondrion_GO), ")") := GO.Term)


#### Select Mitochondria from MitoCarta -----

MitoCarta <- read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/MitoCarta/Human.MitoCarta3.0.xlsx",
                       sheet = "A Human MitoCarta3.0")

MitoCarta_short <- MitoCarta %>% 
  select(UniProt, MitoCarta3.0_MitoPathways, MitoCarta3.0_SubMitoLocalization) %>% 
  filter(UniProt %in% Proteomics_Protein_Info$Uniprot.ID) %>% 
  mutate(MitoCarta3.0_MitoPathways = ifelse(MitoCarta3.0_MitoPathways == "0",
                                            NA, MitoCarta3.0_MitoPathways))
  
MitoCarta_short <- rename(MitoCarta_short, 
                          Uniprot.ID = UniProt,
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
  filter(HGNC.Symbol %in% Proteomics_Protein_Info$Gene.name) %>% 
  distinct(HGNC.Symbol, .keep_all = T)

MitoEvidenceIMPI_short <- rename(MitoEvidenceIMPI_short, 
                                 Protein = HGNC.Symbol,
                                 !!paste0("MitoEvidence_IMPI.Class (", 
                                          nrow(MitoEvidenceIMPI_short), ")")  := IMPI.Class)




## Korrelationsanalyse ----
Timepoints <- c("T2T0", "T4T5", "T5T0", "T4T2", "T5T2", "T4T0")
Welche_Proteine <- c("MB")
All_Proteins <- Proteomics_Protein_Info$Gene.name

KorrelationTable_Results <- Proteomics_Protein_Info
colnames(KorrelationTable_Results)[1] <- "Protein"


KorrelationTable_Results[, "Gene Annotations"] <- NA

## GO Term
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      ProteomicsDataNRE2_Mitochondrion_GO,
                                      by = "Protein")

## MitoCarta
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      MitoCarta_short,
                                      by = "Uniprot.ID")

## MitoEvidenceIMPI
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      MitoEvidenceIMPI_short,
                                      by = "Protein")

KorrelationTable_Results[, "Gene Korrelations"] <- NA

for (Prot in Welche_Proteine){
  KorrelationTable_Results[, paste0(Prot, "__pValue")] <- NA
  KorrelationTable_Results[, paste0(Prot, "__adj.pValue__Sex+Age")] <- NA
  KorrelationTable_Results[, paste0(Prot, "__adj.pValue__Sex+Age+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(Prot, "__adj.pValue__Sex+Age+BMI")] <- NA
  KorrelationTable_Results[, paste0(Prot, "__adj.pValue__Sex+Age+BMI+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(Prot, "__pearsonR")] <- NA
  
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
  
  for(Prot in Welche_Proteine){
    for(Protein in All_Proteins){
      if (Protein == Prot) {
        next
      }
      
      RowNR <- grep(paste0("^", Protein, "$"), KorrelationTable_Results$Protein)
      ColNR_pValue <- grep(paste0(Prot, "__pValue"), colnames(KorrelationTable_Results))
      ColNR_pearsonR <- grep(paste0(Prot, "__pearsonR"), colnames(KorrelationTable_Results))
      ColNR_adj.pValue1 <- grep(paste0(Prot, "__adj.pValue"), colnames(KorrelationTable_Results))[1]
      ColNR_adj.pValue2 <- grep(paste0(Prot, "__adj.pValue"), colnames(KorrelationTable_Results))[2]
      ColNR_adj.pValue3 <- grep(paste0(Prot, "__adj.pValue"), colnames(KorrelationTable_Results))[3]
      ColNR_adj.pValue4 <- grep(paste0(Prot, "__adj.pValue"), colnames(KorrelationTable_Results))[4]
      
      
      Correlation <- cor.test(Table_For_Korrelation[, Protein], 
                              Table_For_Korrelation[, Prot])
      
      LM1 <- summary(lm(Table_For_Korrelation[, Protein] ~ Table_For_Korrelation[, Prot] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0))
      
      LM2 <- summary(lm(Table_For_Korrelation[, Protein] ~ Table_For_Korrelation[, Prot] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$ISIMATS_BL))
      
      LM3 <- summary(lm(Table_For_Korrelation[, Protein] ~ Table_For_Korrelation[, Prot] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0))
      
      LM4 <- summary(lm(Table_For_Korrelation[, Protein] ~ Table_For_Korrelation[, Prot] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0 +
                          Table_For_Korrelation$ISIMATS_BL))
      
      
      
      KorrelationTable_Results[RowNR, ColNR_pValue] <- Correlation[["p.value"]]
      KorrelationTable_Results[RowNR, ColNR_pearsonR] <- Correlation[["estimate"]][["cor"]]
      KorrelationTable_Results[RowNR, ColNR_adj.pValue1] <- LM1[["coefficients"]][2, 4]
      KorrelationTable_Results[RowNR, ColNR_adj.pValue2] <- LM2[["coefficients"]][2, 4]
      KorrelationTable_Results[RowNR, ColNR_adj.pValue3] <- LM3[["coefficients"]][2, 4]
      KorrelationTable_Results[RowNR, ColNR_adj.pValue4] <- LM4[["coefficients"]][2, 4]
    }
  }
  
  KorrelationTable_Results_List[[TP]] <- KorrelationTable_Results
  print(paste0(format(Sys.time(), "%x %T"), "  Finished:  ", TP))
}




## Save Data ----

save(KorrelationTable_Results_List, 
     file = "13_Korrelation_ProteomicsNRE2_mit_Einzelnen_Proteinen.RData")

# load("13_Korrelation_ProteomicsNRE2_mit_Einzelnen_Proteinen.RData")

## Write xlsx ----
## Table with FC
pValue_cols <- grep("pValue", colnames(KorrelationTable_Results))
For_Border <- c(4, 9, grep("__pearsonR", colnames(KorrelationTable_Results)))

# selecting all coll which are numeric
Number_cols <- which(colnames(KorrelationTable_Results) %in% 
                       colnames(KorrelationTable_Results)[unlist(lapply(KorrelationTable_Results, is.numeric))])

Below_0001_list <- list()
for(i in names(KorrelationTable_Results_List)){
  SigGenes <- colSums(KorrelationTable_Results_List[[i]][, pValue_cols] < 0.05)
  
  #### Selecting all pValue-cells which have values <0.0001 -- No Values below 0.0001 -----
  Below_0001 <- lapply(KorrelationTable_Results_List[[i]][, pValue_cols], function (x) x < 0.0001)
  names(Below_0001) <- which(colnames(KorrelationTable_Results_List[[i]]) %in% names(Below_0001))
  for(x in names(Below_0001)){

    if (all(Below_0001[[x]] == FALSE)) {
      Below_0001[[x]] <- NULL
    } else {
      Below_0001[[x]] <- which(Below_0001[[x]] != FALSE)
    }
  }
  Below_0001_list[[i]] <- if (length(Below_0001) != 0) Below_0001
  #### -------------------------------------------------------------------
  
  colnames(KorrelationTable_Results_List[[i]])[pValue_cols] <- paste0(colnames(KorrelationTable_Results_List[[i]])[pValue_cols], 
                                                                      "__(", SigGenes, ")")
  
  colnames(KorrelationTable_Results_List[[i]]) <- sapply(colnames(KorrelationTable_Results_List[[i]]), function (x){
    x <- str_replace_all(x, "__", " ")
  })
}

Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#FF962D") ## medium orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")
Rounding_scientific = createStyle(numFmt = "0.00E+00")


  
wb <- createWorkbook()
addWorksheet(wb, "T5T0")
writeDataTable(wb, sheet = 1, KorrelationTable_Results_List[["T5T0"]])
addWorksheet(wb, "T2T0")
writeDataTable(wb, sheet = 2, KorrelationTable_Results_List[["T2T0"]])
addWorksheet(wb, "T4T5")
writeDataTable(wb, sheet = 3, KorrelationTable_Results_List[["T4T5"]])
addWorksheet(wb, "T4T2")
writeDataTable(wb, sheet = 4, KorrelationTable_Results_List[["T4T2"]])
addWorksheet(wb, "T5T2")
writeDataTable(wb, sheet = 5, KorrelationTable_Results_List[["T5T2"]])
addWorksheet(wb, "T4T0")
writeDataTable(wb, sheet = 6, KorrelationTable_Results_List[["T4T0"]])


## Sheet 1-6
for(sheetNR in names(KorrelationTable_Results_List)){
  for (i in pValue_cols) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = "<0.05", style = Cell_green,
                          type = "expression")
    # conditionalFormatting(wb, sheet = sheetNR, cols = i, 
    #                       rows = 1:nrow(KorrelationTable_Results)+1,
    #                       rule = ">=0.0001", style = Rounding_4,
    #                       type = "between"
    #                       )
    # conditionalFormatting(wb, sheet = sheetNR, cols = i, 
    #                       rows = 1:nrow(KorrelationTable_Results)+1,
    #                       rule = "<0.0001", style = Rounding_scientific,
    #                       type = "expression"
    #                       )
  }
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = '!=""', style = Cell_border_right)
  }
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = '=""', style = Cell_border_right)
  }
  freezePane(wb, sheet = sheetNR, firstCol = TRUE, firstRow = TRUE)
  setColWidths(wb, sheet = sheetNR, cols = 1:ncol(KorrelationTable_Results), widths = 16.5)
  addStyle(wb, sheet = sheetNR, style = style_linebreak, rows = 1, cols = 1:ncol(KorrelationTable_Results), gridExpand = TRUE)
  for(i in Number_cols){
    addStyle(wb, sheet = sheetNR, style = Rounding_4, rows = 1:nrow(KorrelationTable_Results)+1, cols = i)
  }
  for (sheetNR in names(Below_0001_list)) {
    for (Col in names(Below_0001_list[[sheetNR]])){
      for (i in Below_0001_list[[sheetNR]][[Col]]) {
        addStyle(wb, sheet = sheetNR, style = Rounding_scientific, rows = i+1, cols = Col)
      }
    }
  }
}


saveWorkbook(wb, file = "13_Korrelation_ProteomicsNRE2_mit_Einzelnen_Proteinen.xlsx", overwrite = TRUE)

# rm(i, wb)
# rm(pValue_cols)
# rm(For_Border)




