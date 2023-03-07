options(java.parameters = "-Xmx8000m")
library(xlsx)
library(openxlsx)
library(tidyverse)
library(limma)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")


ProteomicsData__Metadata <- as.matrix(ProteomicsDataNRE2[["ProteomicsData__Metadata"]]) %>% t()
colnames(ProteomicsData__Metadata) <- ProteomicsData__Metadata["New_Colname", ]

ProteomicsData__Raw.Log2.Impute <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]

Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]]
row.names(Proteomics_Protein_Info) <- Proteomics_Protein_Info$Uniprot_ID



## Create design matrices for limma ----
condition <- ProteomicsData__Metadata[c("Timepoint", "MyoID", "Rating_Quantile"), ]

design_list <- list()
Subjects_list <- list()
TP_Comparisson <- c("T5T0_FC", "T2T0_FC", "T4T5_FC", "T4T2_FC", "T5T2_FC", "T4T0_FC")
for (i in TP_Comparisson) {
  
  All_Subjectes <- colnames(ProteomicsDataNRE2[[paste0("ProteomicsData__", str_remove(i, "_FC"))]])
  
  a_group = "RES"
  a = unique(condition["MyoID", grep(a_group, condition["Rating_Quantile", ])])
  
  a_subject = NULL
  for(x in a){
    a_subject = c(a_subject ,grep(x, All_Subjectes))
  }
  a_subject = All_Subjectes[a_subject]
  

  b_group = "LRE"
  b = unique(condition["MyoID", grep(b_group, condition["Rating_Quantile", ])])
  b_subject = NULL
  for(x in b){
    b_subject = c(b_subject ,grep(x, All_Subjectes))
  }
  b_subject = All_Subjectes[b_subject]
  

  
  
  # batch_block<-as.vector(condition[2, c(a,b)])
  conditions <- vector(length = length(All_Subjectes))
  names(conditions) <- All_Subjectes
  conditions[names(conditions) %in% a_subject] <- a_group
  conditions[names(conditions) %in% b_subject] <- b_group
  
  conditions <- factor(conditions, levels = c(b_group, a_group))
 
  
  design <- model.matrix(~0 + conditions)
  colnames(design)=levels(as.factor(conditions))
  
  design_list[[i]] <- design
  
  Subjects_list[["LRE-RES"]][[i]] <- c(conditions)

}


## Calculate limma ----

tt.list <- list()
for (i in TP_Comparisson) {
  
  design <- design_list[[i]]
  Subjects <- Subjects_list[["LRE-RES"]][[i]]
  Expression <- ProteomicsDataNRE2[[paste0("ProteomicsData__", str_remove(i, "_FC"))]][, names(Subjects)[!is.na(Subjects)]]
  fit <- lmFit(Expression, design = design)
  
  cont.matrix <- makeContrasts(A=RES-LRE, levels=design)
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  tt=topTable(fit2,"A",number=Inf,p.value=1)
  
  tt.list[[i]] <- tt
  
  xlsx::write.xlsx2(tt, paste("01_Limma analysis Tables/tt.", i, "_LRE-RES__Quantile.xlsx", sep = ""))
  
}

##create dataset with FC ,ratio and linear data ---- 

#order data
data_raw_complete  = NULL
for (TP in TP_Comparisson){
  
  Subjects <- Subjects_list[["LRE-RES"]][[TP]]
  data = ProteomicsDataNRE2[[paste0("ProteomicsData__", str_remove(TP, "_FC"))]][, names(Subjects)[!is.na(Subjects)]]
  data = data[order(row.names(ProteomicsData__Raw.Log2.Impute)), ]
  
  conditions = Subjects_list[["LRE-RES"]][[TP]]
  for (x in names(conditions)) {
    colnames(data)[colnames(data) == x] = paste(colnames(data)[colnames(data) == x], 
                                                ifelse(conditions[x] == 1, "LRE", "RES"), sep = "_")
  }
  
  # data_lin <- 2^data
  # data_raw_complete = cbind(data_raw_complete, "x"=paste0(TP, "_log2"), data, "x"=paste0(TP, "_linear"), data_lin)
  if (is.null(data_raw_complete))
  {
    data_raw_complete = cbind("x"=paste0(TP, "_log2"), data)
  } else {
    data_raw_complete = cbind(data_raw_complete, "x"=paste0(TP, "_log2"), data)
  }
}
Missing_col_name = grep("x", colnames(data_raw_complete))
colnames(data_raw_complete)[Missing_col_name] <- paste(data_raw_complete[1, Missing_col_name])
data_raw_complete[, Missing_col_name] <- NA
rm(Missing_col_name)

data_raw_complete <- as.matrix(apply(data_raw_complete, MARGIN = 2, FUN = as.numeric))

# data <- ProteomicsData__Raw.Log2.Impute[order(rownames(ProteomicsData__Raw.Log2.Impute)), ]
# Log2Data <- rep(NA, dim(data)[1])


# #linear group average
# linearGroupAverages <- rep(NA,dim(data_lin)[1] )
# 
# average <- matrix(nrow = dim(data_lin)[1], ncol = length(unique(condition[1, ])))
# count<-1
# for( i in unique(condition[1, ])){
#   
#   identify_columns <- which(colnames(data_lin) %in% names(which(condition[1, ] == i)))
#   average[, count] <- rowMeans(data_lin[, identify_columns])
#   
#   count<-count + 1
# }
# colnames(average) <- unique(condition[1, ])

#linear group ratios & FC
linearGroupRatios_RESvsLRE <- rep(NA, dim(data)[1])
ratio <- matrix(nrow = dim(data)[1], ncol = length(tt.list) * 3)
Colnames <- NULL
for (i in names(tt.list)) {
  Colnames <- c(Colnames, paste(i, c("Ratio", "Ratio___pValue", "Ratio___adj.pValue"), sep = "___"))
}
colnames(ratio) <- Colnames
row.names(ratio) <- row.names(data)
rm(Colnames)

linearGroupFC_RESvsLRE <- rep(NA, dim(data)[1])
fc <- matrix(nrow = dim(data)[1], ncol = length(tt.list) * 3)
Colnames <- NULL
for (i in names(tt.list)) {
  Colnames <- c(Colnames, paste(i, c("FoldChange", "FoldChange___pValue", "FoldChange___adj.pValue"), sep = "___"))
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
  select(Gene_name, GO.Term) %>% 
  distinct(Gene_name, .keep_all = TRUE)

ProteomicsDataNRE2_Mitochondrion_GO <- rename(ProteomicsDataNRE2_Mitochondrion_GO,
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
                                 Gene_name = HGNC.Symbol,
                                 !!paste0("MitoEvidence_IMPI.Class (", 
                                          nrow(MitoEvidenceIMPI_short), ")")  := IMPI.Class)




## Write xlsx ----
# Proteomics_Protein_Info <- Proteomics_Protein_Info[row.names(data), ]

final_data <- cbind.data.frame("Uniprot_ID" = row.names(data),
                               linearGroupFC_RESvsLRE , fc,
                               linearGroupRatios_RESvsLRE, ratio,
                               data_raw_complete)


final_data <- right_join(Proteomics_Protein_Info,
                         final_data, by = "Uniprot_ID")

final_data[, "Gene Annotations"] <- NA

final_data <- left_join(final_data,
                        ProteomicsDataNRE2_Mitochondrion_GO,
                        by = "Gene_name")

final_data <- left_join(final_data,
                        MitoCarta_short,
                        by = "Uniprot_ID")

final_data <- left_join(final_data,
                        MitoEvidenceIMPI_short,
                        by = "Gene_name")

final_data <- arrange(final_data, Gene_name)


save(final_data, file = "04_Limma_Analyse__Split_RESvsLRE__Quantile.RData")


## selecting all coll which are numeric
Number_cols <- which(colnames(final_data) %in% colnames(final_data)[unlist(lapply(final_data, is.numeric))])

pValue_cols <- grep("pValue", colnames(final_data))
For_Border <- grep("adj.pValue", colnames(final_data))
SigGenes <- colSums(final_data[, pValue_cols] < 0.05)

colnames(final_data)[pValue_cols] <- paste0(colnames(final_data)[pValue_cols], "___(", SigGenes, ")")

colnames(final_data) <- sapply(colnames(final_data), function (x){
  x <- str_replace_all(x, "___", " ")
})


wb <- createWorkbook()
addWorksheet(wb, "Limma Analysis")
writeDataTable(wb, sheet = 1, final_data)

Cell_green <- createStyle(fontColour = "#006100", bgFill = "#c6efce")
Cell_orange <- createStyle(fontColour = "#000000", bgFill = "#ff6700")
Cell_border_right <- createStyle(border = "right", borderStyle = "medium")
style_linebreak <- createStyle(wrapText = TRUE)
Rounding_4 = createStyle(numFmt="0.0000")


for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data)+1,
                        rule = '!=""', style = Cell_border_right)
}
for (i in For_Border) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data)+1,
                        rule = '=""', style = Cell_border_right)
}
for (i in pValue_cols) {
  conditionalFormatting(wb, sheet = 1, cols = i,
                        rows = 1:nrow(final_data)+1,
                        rule = "<0.05", style = Cell_green,
                        type = "expression")
}
freezePane(wb, sheet = 1, firstCol = TRUE, firstRow = TRUE)
setColWidths(wb, sheet = 1, cols = 1:ncol(final_data), widths = 16.5)
addStyle(wb, sheet = 1, style = style_linebreak, rows = 1, cols = 1:ncol(final_data), gridExpand = TRUE)
for(i in Number_cols){
  addStyle(wb, sheet = 1, style = Rounding_4, rows = 1:nrow(final_data)+1, cols = i, gridExpand = TRUE)
}
saveWorkbook(wb, file = "04_Limma_Analyse__Split_RESvsLRE__Quantile.xlsx", overwrite = TRUE)

rm(i, wb)
rm(pValue_cols)
rm(For_Border)




