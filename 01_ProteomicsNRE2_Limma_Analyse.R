options(java.parameters = "-Xmx8000m")
library(xlsx)
library(openxlsx)
library(tidyverse)
library(limma)
library(openxlsx)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")

ProteomicsData__Raw.Log2.Impute <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]
ProteomicsData__Metadata <- as.matrix(ProteomicsDataNRE2[["ProteomicsData__Metadata"]]) %>% t()
colnames(ProteomicsData__Metadata) <- ProteomicsData__Metadata["New_Colname", ]

ProteomicsData__Expression <- 2^ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]

Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]] %>% arrange(Gene_name)
row.names(Proteomics_Protein_Info) <- Proteomics_Protein_Info$Uniprot_ID



## Create design matrices for limma ----
condition <- ProteomicsData__Metadata[c("Timepoint", "MyoID"), ]

design_list <- list()
PairedSubjects_list <- list()
TP_Comparisson <- c("T5T0_FC", "T2T0_FC", "T4T5_FC", "T4T2_FC", "T5T2_FC", "T4T0_FC")
for (i in TP_Comparisson) {
  if (i == "T5T0_FC") {
    
    a = grep("post_resting", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "post_resting"
    
    b = grep("pre_resting", condition["Timepoint", ])
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "pre_resting"
    
  } else if (i == "T2T0_FC") {
    
    a = grep("pre_acuteexercise", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "pre_acuteexercise"
    
    b = grep("pre_resting", condition["Timepoint", ])
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "pre_resting"
    
  } else if (i == "T4T5_FC") {
    
    ## M113 has no post_acuteexercise values
    
    a = grep("post_acuteexercise", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "post_acuteexercise"
    
    M113 = grep("113", condition["MyoID", ])            ## grep M113 cols
    b = grep("post_resting", condition["Timepoint", ])
    b = b[!b %in% M113]                                 ## remove M113 from selected cols
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "post_resting"
    
  } else if (i == "T4T2_FC") {
    
    ## M113 has no post_acuteexercise values
    
    a = grep("post_acuteexercise", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "post_acuteexercise"
    
    M113 = grep("113", condition["MyoID", ])            ## grep M113 cols
    b = grep("pre_acuteexercise", condition["Timepoint", ])
    b = b[!b %in% M113]                                 ## remove M113 from selected cols
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "pre_acuteexercise"
    
  } else if (i == "T5T2_FC") {
    
    a = grep("post_resting", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "post_resting"
    
    b = grep("pre_acuteexercise", condition["Timepoint", ])
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "pre_acuteexercise"
    
  } else if (i == "T4T0_FC") {
    
    ## M113 has no post_acuteexercise values
    
    a = grep("post_acuteexercise", condition["Timepoint", ])
    a_subject = paste(condition["MyoID", a], condition["Timepoint", a], sep = "__")
    a_Timepoint = "post_acuteexercise"
    
    M113 = grep("113", condition["MyoID", ])            ## grep M113 cols
    b = grep("pre_resting", condition["Timepoint", ])
    b = b[!b %in% M113]                                 ## remove M113 from selected cols
    b_subject = paste(condition["MyoID", b], condition["Timepoint", b], sep = "__")
    b_Timepoint = "pre_resting"
    
  }
  
  batch_block <- as.vector(condition[2, c(a,b)])
  conditions <- condition[-2, c(a,b)]
  conditions <- factor(conditions, levels = c(b_Timepoint, a_Timepoint))
  
  
  design <- model.matrix(~batch_block + conditions)
  design_list[[i]] <- design
  
  PairedSubjects_list[["cols"]][[i]] <- c(a, b)
  PairedSubjects_list[["unique"]][[i]] <- c(a_subject, b_subject)
  
}


## Calculate limma ----

tt.list <- list()
for (i in TP_Comparisson) {
  
  design <- design_list[[i]]
  PairedSubjects <- PairedSubjects_list[["cols"]][[i]]
  fit <- lmFit(ProteomicsData__Raw.Log2.Impute[, PairedSubjects], design = design)
  fit <- eBayes(fit)
  
  tt <- topTable(fit, colnames(design)[ncol(design)], number = Inf, p.value = 1)
  
  tt.list[[i]] <- tt
  
  xlsx::write.xlsx2(tt, paste("01_Limma analysis Tables/tt.", i, ".xlsx", sep = ""))
  
}

## create dataset with FC ,ratio and linear data ---- 

#order data
data <- ProteomicsData__Raw.Log2.Impute[order(rownames(ProteomicsData__Raw.Log2.Impute)), ]
Log2Data <- rep(NA, dim(data)[1])

LinearData <- rep(NA, dim(data)[1])
data_lin <- data.matrix(ProteomicsData__Expression[order(rownames(ProteomicsData__Raw.Log2.Impute)), ])

#linear group average
linearGroupAverages <- rep(NA,dim(data_lin)[1] )

average <- matrix(nrow = dim(data_lin)[1], ncol = length(unique(condition[1, ])))
count<-1
for( i in unique(condition[1, ])){
  
  identify_columns <- which(colnames(data_lin) %in% names(which(condition[1, ] == i)))
  average[, count] <- rowMeans(data_lin[, identify_columns])
  
  count<-count + 1
}
colnames(average) <- unique(condition[1, ])

#linear group ratios & FC
linearGroupRatios <- rep(NA, dim(data)[1])
ratio <- matrix(nrow = dim(data)[1], ncol = length(tt.list) * 3)
Colnames <- NULL
for (i in names(tt.list)) {
  Colnames <- c(Colnames, paste(i, c("Ratio", "Ratio__pValue", "Ratio__adj.pValue"), sep = "__"))
}
colnames(ratio) <- Colnames
row.names(ratio) <- row.names(data)
rm(Colnames)

linearGroupFC <- rep(NA, dim(data)[1])
fc <- matrix(nrow = dim(data)[1], ncol = length(tt.list) * 3)
Colnames <- NULL
for (i in names(tt.list)) {
  Colnames <- c(Colnames, paste(i, c("FoldChange", "FoldChange__pValue", "FoldChange__adj.pValue"), sep = "__"))
}
colnames(fc) <- Colnames
row.names(fc) <- row.names(data)
rm(Colnames)


#only for pre-definde comparisons in the same order

# expression <- expression[rownames(data_lin),] ## order rows like in data_lin (same as in average)
for(i in TP_Comparisson){
  
  # LogFC aus tt Tabellen
  tt <- tt.list[[i]]
  tt <- tt[rownames(data), ]
  
  ColNR <- grep(i, colnames(ratio))
  
  ratio[, ColNR[1]] <- tt$logFC
  ratio[, ColNR[1]] <- round(2^ratio[, ColNR[1]], digits = 6)
  
  ratio[, ColNR[2]] <- tt$P.Value
  ratio[, ColNR[3]] <- tt$adj.P.Val
  
  
  for( u in 1:dim(ratio)[1]){
    if(ratio[u, ColNR[1]]>1){
      fc[u, ColNR[1]]<-round(ratio[u, ColNR[1]], 6)
    }else{
      fc[u, ColNR[1]]<-round(-1/ratio[u, ColNR[1]], 6)
    }
  }
  
  fc[, ColNR[2]] <- tt$P.Value
  fc[, ColNR[3]] <- tt$adj.P.Val
  
}

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




## Write xlsx ----
Proteomics_Protein_Info <- Proteomics_Protein_Info[row.names(data), ]

final_data <- cbind.data.frame("Uniprot_ID" = row.names(data),
                               linearGroupFC , fc,
                               linearGroupRatios, ratio,
                               linearGroupAverages, average)

final_data <- right_join(Proteomics_Protein_Info,
                                  final_data, by = "Uniprot_ID")

final_data[, "Gene Annotations"] <- NA

final_data <- left_join(final_data,
                        ProteomicsDataNRE2_Mitochondrion_GO,
                        by = "Uniprot_ID")

final_data <- left_join(final_data,
                        MitoCarta_short,
                        by = "Uniprot_ID")

final_data <- left_join(final_data,
                        MitoEvidenceIMPI_short,
                        by = "Uniprot_ID")

final_data_complete <- final_data

final_data_raw <- cbind.data.frame(Log2Data, data,
                                   LinearData, data_lin)



save(final_data_complete, file = "01_ProteomicsNRE2_Limma_Analyse.RData")

if(table(final_data_complete$Uniprot_ID == row.names(final_data_raw))["TRUE"] == nrow(final_data_complete)){
  final_data_complete <- cbind(final_data_complete, final_data_raw)
}

# Arrnage alphabetical
Proteomics_Protein_Info <- arrange(Proteomics_Protein_Info, Gene_name)
final_data_complete <- final_data_complete[Proteomics_Protein_Info$Uniprot_ID, ]

Cells_for_Orange <- c("linearGroupFC", "linearGroupRatios", "linearGroupAverages", "Log2Data", "LinearData", "Gene Annotations")

## selecting all coll which are numeric
Number_cols <- which(colnames(final_data_complete) %in% colnames(final_data_complete)[unlist(lapply(final_data_complete, is.numeric))])

pValue_cols <- grep("pValue", colnames(final_data_complete))
For_Border <- grep("adj.pValue", colnames(final_data_complete))
SigGenes <- colSums(final_data_complete[, pValue_cols] < 0.05)

## Selecting all pValue-cells which have values <0.0001
Below_0001 <- lapply(final_data_complete[, pValue_cols], function (x) x < 0.0001)
names(Below_0001) <- which(colnames(final_data_complete) %in% names(Below_0001))
for(i in names(Below_0001)){
  
  if (all(Below_0001[[i]] == FALSE)) {
    Below_0001[[i]] <- NULL
  } else {
    Below_0001[[i]] <- which(row.names(final_data_complete) %in% row.names(final_data_complete)[Below_0001[[i]]])
  }
}



colnames(final_data_complete)[pValue_cols] <- paste0(colnames(final_data_complete)[pValue_cols], "__(", SigGenes, ")")

colnames(final_data_complete) <- sapply(colnames(final_data_complete), function (x){
  x <- str_replace_all(x, "__", " ")
})


wb <- createWorkbook()
addWorksheet(wb, "Limma Analysis")
writeData(wb, sheet = 1, final_data_complete)

Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
# Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700") ## lightning orange
Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#FF962D") ## medium orange
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt = "0.0000")
Rounding_scientific = createStyle(numFmt = "0.00E+00")


for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data_complete)+1,
                        rule = '!=""', style = Cell_border_right)
}
for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data_complete)+1,
                        rule = '=""', style = Cell_border_right)
}
for (i in pValue_cols) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data_complete)+1,
                        rule = "<0.05", style = Cell_green,
                        type = "expression")
}
for (i in Cells_for_Orange){
  conditionalFormatting(wb, sheet = 1, cols = 1:ncol(final_data_complete),
                        rows = 1,
                        rule = i, style = Cell_orange,
                        type = "contains")
}
freezePane(wb, sheet = 1, firstCol = TRUE, firstRow = TRUE)
setColWidths(wb, sheet = 1, cols = 1:ncol(final_data_complete), widths = 16.5)
addStyle(wb, sheet = 1, style = style_linebreak, rows = 1, cols = 1:ncol(final_data_complete), gridExpand = TRUE)
for(i in Number_cols){
  addStyle(wb, sheet = 1, style = Rounding_4, rows = 1:nrow(final_data_complete)+1, cols = i)
}
for (col in names(Below_0001)) {
  for (i in Below_0001[[col]]) {
    addStyle(wb, sheet = 1, style = Rounding_scientific, rows = i+1, cols = col)
  }
}
saveWorkbook(wb, file = "01_ProteomicsNRE2_Limma_Analyse.xlsx", overwrite = TRUE)

# rm(i, wb)
# rm(pValue_cols)
# rm(For_Border)




