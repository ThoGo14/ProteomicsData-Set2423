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

# row.names(Proteomics_Protein_Info) <- Proteomics_Protein_Info$Gene_name

KlinischeDaten <- read.xlsx("E:/Auswertung mit R/00_Helpful Files/00_Klinische_Daten.xlsx")

## Load Limma-analysis table for FC
load("01_ProteomicsNRE2_Limma_Analyse.RData")
Limma_FC <- final_data_complete[, c(2, grep("__FoldChange$", colnames(final_data_complete)))] %>% 
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
  select(Uniprot_ID, GO.Term) %>% 
  distinct(Uniprot_ID, .keep_all = TRUE)

ProteomicsDataNRE2_Mitochondrion_GO <- rename(ProteomicsDataNRE2_Mitochondrion_GO,
                                              !!paste0("Mitochondrion (GO:0005739) (",
                                                       nrow(ProteomicsDataNRE2_Mitochondrion_GO), ")") := GO.Term)


#### Select Mitochondria from MitoCarta -----

MitoCarta <- openxlsx::read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/MitoCarta/Human.MitoCarta3.0.xlsx",
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

MitoEvidenceIMPI <- openxlsx::read.xlsx("E:/Auswertung mit R/00_Gene_Annotation/MitoEvidenceIMPI/impi-2021-q4pre-20211001-dist_0.xlsx",
                                        sheet = "IMPI-2021Q4pre")

MitoEvidenceIMPI$HGNC.Symbol <- ifelse(is.na(MitoEvidenceIMPI$HGNC.Symbol),
                                       MitoEvidenceIMPI$Symbol, MitoEvidenceIMPI$HGNC.Symbol)

MitoEvidenceIMPI_short <- MitoEvidenceIMPI %>% 
  select(HGNC.Symbol, IMPI.Class) %>% 
  filter(HGNC.Symbol %in% Proteomics_Protein_Info$Gene_name) %>% 
  distinct(HGNC.Symbol, .keep_all = T)

MitoEvidenceIMPI_short <- rename(MitoEvidenceIMPI_short, 
                                 Uniprot_ID = HGNC.Symbol,
                                 !!paste0("MitoEvidence_IMPI.Class (", 
                                          nrow(MitoEvidenceIMPI_short), ")")  := IMPI.Class)




## Korrelationsanalyse ----
Timepoints <- c("T2T0", "T4T5", "T5T0", "T4T2", "T5T2", "T4T0")
Welche_KD <- c("ISI_FC", "SM_MODP.S_FC", "IAS_ergo_BW_FC")
Welche_KD_TTest <- c("SEX")
All_Proteins <- Proteomics_Protein_Info$Uniprot_ID

KorrelationTable_Results <- Proteomics_Protein_Info[, -1]

KorrelationTable_Results[, "Limma__Analysis"] <- NA
KorrelationTable_Results <- left_join(KorrelationTable_Results, 
                                      Limma_FC,
                                      by = "Gene_name")

KorrelationTable_Results[, "Correlation__Analysis"] <- NA

for (KD in Welche_KD){
  
  KorrelationTable_Results[, paste0(KD, "__pValue")] <- NA
  KorrelationTable_Results[, paste0(KD, "__adj.pValue__Sex+Age")] <- NA
  KorrelationTable_Results[, paste0(KD, "__adj.pValue__Sex+Age+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(KD, "__adj.pValue__Sex+Age+BMI")] <- NA
  KorrelationTable_Results[, paste0(KD, "__adj.pValue__Sex+Age+BMI+ISI_BL")] <- NA
  KorrelationTable_Results[, paste0(KD, "__pearsonR")] <- NA
  
}

KorrelationTable_Results[, "TTest_Analysis"] <- NA
for (KD in Welche_KD_TTest){
  
  KorrelationTable_Results[, paste0(KD, "__pValue")] <- NA
  
}
KorrelationTable_Results[, "Gene Annotations"] <- NA

## GO Term
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      ProteomicsDataNRE2_Mitochondrion_GO,
                                      by = "Uniprot_ID")

## MitoCarta
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      MitoCarta_short,
                                      by = "Uniprot_ID")

## MitoEvidenceIMPI
KorrelationTable_Results <- left_join(KorrelationTable_Results,
                                      MitoEvidenceIMPI_short,
                                      by = "Uniprot_ID")

KorrelationTable_Results <- KorrelationTable_Results %>% arrange(Gene_name) %>% 
  filter(!duplicated(Uniprot_ID))


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
    for(Uniprot_ID in All_Proteins){
      RowNR <- grep(paste0("^", Uniprot_ID, "$"), KorrelationTable_Results$Uniprot_ID)
      ColNR_pValue <- grep(paste0(KD, "__pValue"), colnames(KorrelationTable_Results))
      
      ColNR_pearsonR <- grep(paste0(KD, "__pearsonR"), colnames(KorrelationTable_Results))
      ColNR_adj.pValue1 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results))[1]
      ColNR_adj.pValue2 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results))[2]
      ColNR_adj.pValue3 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results))[3]
      ColNR_adj.pValue4 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results))[4]
      
      
      Correlation <- cor.test(Table_For_Korrelation[, Uniprot_ID], 
                              Table_For_Korrelation[, KD])
      
      LM1 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0))
      
      LM2 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$ISIMATS_BL))
      
      LM3 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD] +
                          Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0))
      
      LM4 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD] +
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
  for(KD in Welche_KD_TTest){
    for(Uniprot_ID in All_Proteins){
      RowNR <- grep(paste0("^", Uniprot_ID, "$"), KorrelationTable_Results$Uniprot_ID)
      ColNR_pValue <- grep(paste0(KD, "__pValue"), colnames(KorrelationTable_Results))
      
      if (KD == "SEX") {
        Group_0 <- Table_For_Korrelation %>% select(all_of(Uniprot_ID), SEX) %>% filter(SEX == 0) %>% pull(all_of(Uniprot_ID))
        Group_1 <- Table_For_Korrelation %>% select(all_of(Uniprot_ID), SEX) %>% filter(SEX == 1) %>% pull(all_of(Uniprot_ID))
        
        TTest <- t.test(Group_0, Group_1, paired = FALSE)
      }
      
      KorrelationTable_Results[RowNR, ColNR_pValue] <- TTest[["p.value"]]
    }
  }
  
  KorrelationTable_Results_List[[TP]] <- KorrelationTable_Results
  print(paste0(format(Sys.time(), "%x %T"), "  Finished:  ", TP))
}



## Korrelationsanalyse - Baseline ----
TP_Table <- data.frame(TP_old = c("T0", "T2", "T4", "T5"),
                       TP_new = c("pre_resting", "pre_acuteexercise",
                                  "post_acuteexercise", "post_resting"))


Welche_KD <- c("ISIMATS_BL__T0", "ISIMATS_TR__T5", 
               "SM_MODP.S_V0__T0", "SM_MODP.S_V0__T2", "SM_MODP.S_V1__T5", "SM_MODP.S_V1__T4",
               "IAS_ergo_BW_V0__T0", "IAS_ergo_BW_V0__T2", "IAS_ergo_BW_V1__T5", "IAS_ergo_BW_V1__T4")

Welche_KD_TTest <- c("SEX__T0", "SEX__T2", "SEX__T4", "SEX__T5")
All_Proteins <- Proteomics_Protein_Info$Uniprot_ID

KorrelationTable_Results_Baseline <- Proteomics_Protein_Info[, -1]

KorrelationTable_Results_Baseline[ "Limma__Analysis"] <- NA
KorrelationTable_Results_Baseline <- left_join(KorrelationTable_Results_Baseline, 
                                               Limma_FC,
                                               by = "Gene_name")

KorrelationTable_Results_Baseline[, "Correlation__Analysis"] <- NA

for (KD in Welche_KD){
  KorrelationTable_Results_Baseline[, paste0(KD, "__pValue")] <- NA
  KorrelationTable_Results_Baseline[, paste0(KD, "__adj.pValue__Sex+Age")] <- NA
  KorrelationTable_Results_Baseline[, paste0(KD, "__adj.pValue__Sex+Age+BMI")] <- NA
  # KorrelationTable_Results_Baseline[, paste0(KD, "__adj.pValue__Sex+Age+ISI_BL")] <- NA
  # KorrelationTable_Results_Baseline[, paste0(KD, "__adj.pValue__Sex+Age+BMI+ISI_BL")] <- NA
  KorrelationTable_Results_Baseline[, paste0(KD, "__pearsonR")] <- NA
  
}

KorrelationTable_Results_Baseline[, "TTest_Analysis"] <- NA
for (KD in Welche_KD_TTest){
  
  KorrelationTable_Results_Baseline[, paste0(KD, "__pValue")] <- NA
  
}

KorrelationTable_Results_Baseline[, "Gene Annotations"] <- NA

## GO Term
KorrelationTable_Results_Baseline <- left_join(KorrelationTable_Results_Baseline,
                                               ProteomicsDataNRE2_Mitochondrion_GO,
                                               by = "Uniprot_ID")

## MitoCarta
KorrelationTable_Results_Baseline <- left_join(KorrelationTable_Results_Baseline,
                                               MitoCarta_short,
                                               by = "Uniprot_ID")

## MitoEvidenceIMPI
KorrelationTable_Results_Baseline <- left_join(KorrelationTable_Results_Baseline,
                                               MitoEvidenceIMPI_short,
                                               by = "Uniprot_ID")

KorrelationTable_Results_Baseline <- KorrelationTable_Results_Baseline %>% arrange(Gene_name) %>% 
  filter(!duplicated(Uniprot_ID))

for(KD in Welche_KD){
  
  TP = unlist(str_split(KD, "__"))[2]
  TP = TP_Table[grep(TP, TP_Table$TP_old), "TP_new"]
  
  TP_Cols <- grep(TP, colnames(ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]))
  
  ProteinTable <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]][, TP_Cols]
  colnames(ProteinTable) <- sapply(colnames(ProteinTable), function(x){
    x = unlist(str_split(x, "_"))[1]
    return(x)
  })
  ProteinTable <- data.frame("MyoID" = colnames(ProteinTable), t(ProteinTable), check.names = F)
  
  Table_For_Korrelation <- right_join(KlinischeDaten, ProteinTable, by = "MyoID")
  
  KD_new <- unlist(str_split(KD, "__"))[1]
  
  for(Uniprot_ID in All_Proteins){
    RowNR <- grep(paste0("^", Uniprot_ID, "$"), KorrelationTable_Results_Baseline$Uniprot_ID)
    ColNR_pValue <- grep(paste0(KD, "__pValue"), colnames(KorrelationTable_Results_Baseline))
    ColNR_pearsonR <- grep(paste0(KD, "__pearsonR"), colnames(KorrelationTable_Results_Baseline))
    ColNR_adj.pValue1 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results_Baseline))[1]
    ColNR_adj.pValue2 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results_Baseline))[2]
    # ColNR_adj.pValue3 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results_Baseline))[3]
    # ColNR_adj.pValue4 <- grep(paste0(KD, "__adj.pValue"), colnames(KorrelationTable_Results_Baseline))[4]
    
    
    Correlation <- cor.test(Table_For_Korrelation[, Uniprot_ID], 
                            Table_For_Korrelation[, KD_new])
    
    LM1 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD_new] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0))
    
    
    
    LM2 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD_new] +
                        Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0))
    
    
    # LM3 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD_new] +
    #                     Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$ISIMATS_BL))
    # 
    # LM4 <- summary(lm(Table_For_Korrelation[, Uniprot_ID] ~ Table_For_Korrelation[, KD_new] +
    #                     Table_For_Korrelation$SEX + Table_For_Korrelation$AGE_V0 + Table_For_Korrelation$BMI_V0 +
    #                     Table_For_Korrelation$ISIMATS_BL))
    
    KorrelationTable_Results_Baseline[RowNR, ColNR_pValue] <- Correlation[["p.value"]]
    KorrelationTable_Results_Baseline[RowNR, ColNR_pearsonR] <- Correlation[["estimate"]][["cor"]]
    KorrelationTable_Results_Baseline[RowNR, ColNR_adj.pValue1] <- LM1[["coefficients"]][2, 4]
    KorrelationTable_Results_Baseline[RowNR, ColNR_adj.pValue2] <- LM2[["coefficients"]][2, 4]
    # KorrelationTable_Results_Baseline[RowNR, ColNR_adj.pValue3] <- round(LM3[["coefficients"]][2, 4], digits = 4)
    # KorrelationTable_Results_Baseline[RowNR, ColNR_adj.pValue4] <- round(LM4[["coefficients"]][2, 4], digits = 4)
  }
  
  print(paste0(format(Sys.time(), "%x %T"), "  Finished:  ", KD))
}

for(KD in Welche_KD_TTest){
  
  
  TP = unlist(str_split(KD, "__"))[2]
  TP = TP_Table[grep(TP, TP_Table$TP_old), "TP_new"]
  
  TP_Cols <- grep(TP, colnames(ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]))
  
  ProteinTable <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]][, TP_Cols]
  colnames(ProteinTable) <- sapply(colnames(ProteinTable), function(x){
    x = unlist(str_split(x, "_"))[1]
    return(x)
  })
  ProteinTable <- data.frame("MyoID" = colnames(ProteinTable), t(ProteinTable), check.names = F)
  
  Table_For_Korrelation <- right_join(KlinischeDaten, ProteinTable, by = "MyoID")
  
  KD_new <- unlist(str_split(KD, "__"))[1]
  
  
  for(Uniprot_ID in All_Proteins){
    RowNR <- grep(paste0("^", Uniprot_ID, "$"), KorrelationTable_Results_Baseline$Uniprot_ID)
    ColNR_pValue <- grep(paste0(KD, "__pValue"), colnames(KorrelationTable_Results_Baseline))
    
    if (KD_new == "SEX") {
      Group_0 <- Table_For_Korrelation %>% select(all_of(Uniprot_ID), SEX) %>% filter(SEX == 0) %>% pull(all_of(Uniprot_ID))
      Group_1 <- Table_For_Korrelation %>% select(all_of(Uniprot_ID), SEX) %>% filter(SEX == 1) %>% pull(all_of(Uniprot_ID))
      
      TTest <- t.test(Group_0, Group_1, paired = FALSE)
    }
    
    KorrelationTable_Results_Baseline[RowNR, ColNR_pValue] <- TTest[["p.value"]]
  }
}


## Save Data ----

save(KorrelationTable_Results_List, KorrelationTable_Results_Baseline, file = "05_Korrelation_ProteomicsNRE2_mit_KlinischenDaten.RData")


## Write xlsx ----
## Table with FC
pValue_cols <- grep("pValue", colnames(KorrelationTable_Results))
For_Border <- c(5, 12, grep("__pearsonR", colnames(KorrelationTable_Results)),
                grep("Gene Annotations", names(KorrelationTable_Results))-1, ncol(KorrelationTable_Results))
Limma_cols <- grep(" FoldChange ", colnames(KorrelationTable_Results))


# selecting all coll which are numeric
Number_cols <- which(colnames(KorrelationTable_Results) %in% 
                       colnames(KorrelationTable_Results)[unlist(lapply(KorrelationTable_Results, is.numeric))])

Number_cols <- setdiff(Number_cols, Limma_cols)

Below_0001_list <- list()
for(i in names(KorrelationTable_Results_List)){
  SigGenes <- colSums(KorrelationTable_Results_List[[i]][, pValue_cols] < 0.05)
  
  
  
  #### Selecting all pValue-cells which have values <0.0001 -- No Values below 0.0001 -----
  # Below_0001 <- lapply(KorrelationTable_Results_List[[i]][, pValue_cols], function (x) x < 0.0001)
  # names(Below_0001) <- which(colnames(KorrelationTable_Results_List[[i]]) %in% names(Below_0001))
  # for(x in names(Below_0001)){
  #   
  #   if (all(Below_0001[[x]] == FALSE)) {
  #     Below_0001[[x]] <- NULL
  #   } else {
  # Below_0001[[x]] <- which(row.names(KorrelationTable_Results_List[[x]]) %in%
  #                            row.names(KorrelationTable_Results_List[[x]])[Below_0001[[x]]])
  #   }
  # }
  # Below_0001_list[[i]] <- if (length(Below_0001 != 0)) Below_0001
  #### -------------------------------------------------------------------
  
  colnames(KorrelationTable_Results_List[[i]])[pValue_cols] <- paste0(colnames(KorrelationTable_Results_List[[i]])[pValue_cols], 
                                                                      "__(", SigGenes, ")")
  
  
  colnames(KorrelationTable_Results_List[[i]]) <- sapply(colnames(KorrelationTable_Results_List[[i]]), function (x){
    x <- str_replace_all(x, "__", " ")
  })
}



  
for(i in 1:length(KorrelationTable_Results_List)){  
  
}

## Table with Basline Correlations
pValue_cols2 <- grep("pValue", colnames(KorrelationTable_Results_Baseline))
For_Border2 <- c(5, 12, grep("__pearsonR", colnames(KorrelationTable_Results_Baseline)), 
                 grep("Gene Annotations", names(KorrelationTable_Results_Baseline))-1, ncol(KorrelationTable_Results_Baseline))
Limma_cols2 <- grep(" FoldChange ", colnames(KorrelationTable_Results_Baseline))

SigGenes <- colSums(KorrelationTable_Results_Baseline[, pValue_cols2] < 0.05)

Number_cols2 <- which(colnames(KorrelationTable_Results_Baseline) %in% 
                       colnames(KorrelationTable_Results_Baseline)[unlist(lapply(KorrelationTable_Results_Baseline, is.numeric))])

Number_cols2 <- setdiff(Number_cols2, Limma_cols2)
# Selecting all pValue-cells which have values <0.0001
Below_0001_2 <- lapply(KorrelationTable_Results_Baseline[, pValue_cols2], function (x) x < 0.0001)
names(Below_0001_2) <- which(colnames(KorrelationTable_Results_Baseline) %in% names(Below_0001_2))
for(x in names(Below_0001_2)){

  if (all(Below_0001_2[[x]] == FALSE)) {
    Below_0001_2[[x]] <- NULL
  } else {
    Below_0001_2[[x]] <- which(row.names(KorrelationTable_Results_Baseline) %in% 
                                 row.names(KorrelationTable_Results_Baseline)[Below_0001_2[[x]]])
  }
}

colnames(KorrelationTable_Results_Baseline)[pValue_cols2] <- paste0(colnames(KorrelationTable_Results_Baseline)[pValue_cols2], 
                                                                   "__(", SigGenes, ")")

colnames(KorrelationTable_Results_Baseline) <- sapply(colnames(KorrelationTable_Results_Baseline), function (x){
    x <- str_replace_all(x, "__", " ")
  })

  
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
addWorksheet(wb, "Baselines")
writeDataTable(wb, sheet = 7, KorrelationTable_Results_Baseline)


## Sheet 1-6
for(sheetNR in 1:length(KorrelationTable_Results_List)){
  for (i in pValue_cols) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = "<0.05", style = Cell_green,
                          type = "expression")
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = "<=1", style = Rounding_4,
                          type = "expression"
                          )
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = "<0.0001", style = Rounding_scientific,
                          type = "expression"
                          )
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
  for (col in names(Limma_cols)){
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(KorrelationTable_Results)+1,
                          rule = '<WhichFC', style = Cell_border_right)
  }
  freezePane(wb, sheet = sheetNR, firstCol = TRUE, firstRow = TRUE)
  setColWidths(wb, sheet = sheetNR, cols = 1:ncol(KorrelationTable_Results), widths = 16.5)
  addStyle(wb, sheet = sheetNR, style = style_linebreak, rows = 1, cols = 1:ncol(KorrelationTable_Results), gridExpand = TRUE)
  for(i in Number_cols){
    addStyle(wb, sheet = sheetNR, style = Rounding_4, rows = 1:nrow(KorrelationTable_Results)+1, cols = i)
  }
}

## Sheet 7
for (i in pValue_cols2) {
  conditionalFormatting(wb, sheet = 7, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_Baseline)+1,
                        rule = "<0.05", style = Cell_green,
                        type = "expression")
  conditionalFormatting(wb, sheet = 7, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_Baseline)+1,
                        rule = "<=1", style = Rounding_4,
                        type = "expression"
  )
  conditionalFormatting(wb, sheet = 7, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_Baseline)+1,
                        rule = "<0.0001", style = Rounding_scientific,
                        type = "expression"
  )
}
for (i in For_Border2) {
  conditionalFormatting(wb, sheet = 7, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_Baseline)+1,
                        rule = '!=""', style = Cell_border_right)
}
for (i in For_Border2) {
  conditionalFormatting(wb, sheet = 7, cols = i, 
                        rows = 1:nrow(KorrelationTable_Results_Baseline)+1,
                        rule = '=""', style = Cell_border_right)
}

freezePane(wb, sheet = 7, firstCol = TRUE, firstRow = TRUE)
setColWidths(wb, sheet = 7, cols = 1:ncol(KorrelationTable_Results_Baseline), widths = 16.5)
addStyle(wb, sheet = 7, style = style_linebreak, rows = 1, cols = 1:ncol(KorrelationTable_Results_Baseline), gridExpand = TRUE)
for(i in Number_cols2){
  addStyle(wb, sheet = 7, style = Rounding_4, rows = 1:nrow(KorrelationTable_Results_Baseline)+1, cols = i)
}
for (col in names(Below_0001_2)) {
  for (i in Below_0001_2[[col]]) {
    addStyle(wb, sheet = 7, style = Rounding_scientific, rows = i+1, cols = col)
  }
}
saveWorkbook(wb, file = "05_Korrelation_ProteomicsNRE2_mit_KlinischenDaten.xlsx", overwrite = TRUE)

# rm(i, wb)
# rm(pValue_cols)
# rm(For_Border)




