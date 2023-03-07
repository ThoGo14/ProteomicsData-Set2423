library(openxlsx)
library(tidyverse)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


load("01_ProteomicsNRE2_Limma_Analyse.RData")

#### Timepoint comparissions ----
Which_Timepoints <- c("T4T2__T5T2", "T4T2__T5T0")


Overlap_Results_List <- list()


for(TP in Which_Timepoints){
  TP1 = unlist(str_split(TP, "__"))[1]
  TP2 = unlist(str_split(TP, "__"))[2]
  
  TP1_Cols <- colnames(final_data_complete)[grep(paste0(TP1, "_FC__FoldChange"), colnames(final_data_complete))]
  TP2_Cols <- colnames(final_data_complete)[grep(paste0(TP2, "_FC__FoldChange"), colnames(final_data_complete))]
  
  
  SigProteins_overlap <- intersect(final_data_complete[final_data_complete[, TP1_Cols[2]] < 0.05, "Gene_name"],
                                   final_data_complete[final_data_complete[, TP2_Cols[2]] < 0.05, "Gene_name"])
  

  Overlap_Results <- data.frame(Gene_name = SigProteins_overlap)
  Overlap_Results <- left_join(Overlap_Results,
                               final_data_complete[, c("Gene_name", TP1_Cols)],
                               by = "Gene_name")
  Overlap_Results <- left_join(Overlap_Results,
                               final_data_complete[, c("Gene_name", TP2_Cols)],
                               by = "Gene_name")
  
  FC1 = 1.1
  FC2 = 1.2
  
  Overlap_Results[, paste0("Both FC>", FC1, " & same direction")] <- 
    (Overlap_Results[, TP1_Cols[1]] > FC1 & Overlap_Results[, TP2_Cols[1]] > FC1) | 
    (Overlap_Results[, TP1_Cols[1]] < -FC1 & Overlap_Results[, TP2_Cols[1]] < -FC1)
  
  Overlap_Results[, paste0("Both FC>", FC2, " & same direction")] <- 
    (Overlap_Results[, TP1_Cols[1]] > FC2 & Overlap_Results[, TP2_Cols[1]] > FC2) | 
    (Overlap_Results[, TP1_Cols[1]] < -FC2 & Overlap_Results[, TP2_Cols[1]] < -FC2)
  
  
  ## Change Colnames
  
  colnames(Overlap_Results)[1] <- paste0(colnames(Overlap_Results)[1], " (", nrow(Overlap_Results), ")")
  for (Col in c(8:9)) {
    TRUE_cols = table(Overlap_Results[, Col])["TRUE"]
    colnames(Overlap_Results)[Col] <- paste0(colnames(Overlap_Results)[Col], " (", TRUE_cols, ")")
  }
  
  
  ## Save in List
  Overlap_Results_List[[TP]] <- Overlap_Results
}


## Save Data ----

save(Overlap_Results_List, 
     file = "14_Overlap_of_Limma_Timepoints.RData")

# load("14_Overlap_of_Limma_Timepoints.RData")

## Write xlsx ----
## Table with FC
pValue_cols <- grep("pValue", colnames(Overlap_Results))
For_Border <- c(4, 9, grep("__pearsonR", colnames(Overlap_Results)))

# selecting all coll which are numeric
Number_cols <- which(colnames(Overlap_Results) %in% 
                       colnames(Overlap_Results)[unlist(lapply(Overlap_Results, is.numeric))])

Below_0001_list <- list()
for(i in names(Overlap_Results_List)){
  SigGenes <- colSums(Overlap_Results_List[[i]][, pValue_cols] < 0.05)
  
  #### Selecting all pValue-cells which have values <0.0001 -- No Values below 0.0001 -----
  Below_0001 <- lapply(Overlap_Results_List[[i]][, pValue_cols], function (x) x < 0.0001)
  names(Below_0001) <- which(colnames(Overlap_Results_List[[i]]) %in% names(Below_0001))
  for(x in names(Below_0001)){
    
    if (all(Below_0001[[x]] == FALSE)) {
      Below_0001[[x]] <- NULL
    } else {
      Below_0001[[x]] <- which(Below_0001[[x]] != FALSE)
    }
  }
  Below_0001_list[[i]] <- if (length(Below_0001) != 0) Below_0001
  #### -------------------------------------------------------------------
  
  colnames(Overlap_Results_List[[i]])[pValue_cols] <- paste0(colnames(Overlap_Results_List[[i]])[pValue_cols], 
                                                                      "__(", SigGenes, ")")
  
  colnames(Overlap_Results_List[[i]]) <- sapply(colnames(Overlap_Results_List[[i]]), function (x){
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

for(sheetNR in names(Overlap_Results_List)){
  addWorksheet(wb, sheetName = sheetNR)
  writeDataTable(wb, sheet = sheetNR, Overlap_Results_List[[sheetNR]])
  
  for (i in pValue_cols) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(Overlap_Results_List[[sheetNR]])+1,
                          rule = "<0.05", style = Cell_green,
                          type = "expression")
  }
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(Overlap_Results_List[[sheetNR]])+1,
                          rule = '!=""', style = Cell_border_right)
  }
  for (i in For_Border) {
    conditionalFormatting(wb, sheet = sheetNR, cols = i, 
                          rows = 1:nrow(Overlap_Results_List[[sheetNR]])+1,
                          rule = '=""', style = Cell_border_right)
  }
  freezePane(wb, sheet = sheetNR, firstCol = TRUE, firstRow = TRUE)
  setColWidths(wb, sheet = sheetNR, cols = 1:ncol(Overlap_Results_List[[sheetNR]]), widths = 16.5)
  addStyle(wb, sheet = sheetNR, style = style_linebreak, rows = 1, cols = 1:ncol(Overlap_Results_List[[sheetNR]]), gridExpand = TRUE)
  for(i in Number_cols){
    addStyle(wb, sheet = sheetNR, style = Rounding_4, rows = 1:nrow(Overlap_Results_List[[sheetNR]])+1, cols = i)
  }
  if (sheetNR %in% names(Below_0001_list)) {
    for (Col in names(Below_0001_list[[sheetNR]])){
      for (i in Below_0001_list[[sheetNR]][[Col]]) {
        addStyle(wb, sheet = sheetNR, style = Rounding_scientific, rows = i+1, cols = Col)
      }
    }
  }
}


saveWorkbook(wb, file = "14_Overlap_of_Limma_Timepoints.xlsx", overwrite = TRUE)

# rm(i, wb)
# rm(pValue_cols)
# rm(For_Border)




