library(openxlsx)
library(tidyverse)
library(ggsignif)
library(gridExtra)
theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## Load Data and Convert Dataframes in matrix ----
load("ProteomicsNRE2__Expression.RData")
LimmaAnalysis <- read.xlsx("01_ProteomicsNRE2_Limma_Analyse.xlsx")

ProteomicsData__Metadata <- as.matrix(ProteomicsDataNRE2[["ProteomicsData__Metadata"]]) %>% t()
colnames(ProteomicsData__Metadata) <- ProteomicsData__Metadata["New_Colname", ]

Proteomics_Protein_Info <- ProteomicsDataNRE2[["Proteomics_Protein_Info"]] %>% arrange(Gene_name)


ProteomicsData__Raw.Log2.Impute <- ProteomicsDataNRE2[["ProteomicsData__Raw.Log2.Impute"]]
ProteomicsData__Raw.Log2.Impute <- ProteomicsData__Raw.Log2.Impute[Proteomics_Protein_Info$Uniprot_ID, ]


Resp_Group = ProteomicsData__Metadata["MyoID", ]
names(Resp_Group) = ProteomicsData__Metadata["Rating_Quantile", ]
Resp_Group <- Resp_Group[!duplicated(Resp_Group)]
Resp_Group <- data.frame(Group = names(Resp_Group), MyoID = Resp_Group)


## Prepare Data as long table for Plots -----
ProteomicsData__Expression_long <- data.frame("Gene" = row.names(ProteomicsData__Raw.Log2.Impute), 
                                                   ProteomicsData__Raw.Log2.Impute)

ProteomicsData__Expression_long <- gather(data = ProteomicsData__Expression_long, 
                                               key = "Colname", value = "value",
                                               2:ncol(ProteomicsData__Expression_long))

ProteomicsData__Expression_long$MyoID <- sapply(ProteomicsData__Expression_long$Colname, function (x){
  x = unlist(str_split(x, "_"))[1]
  return(x)
})
ProteomicsData__Expression_long$Timepoint <- sapply(ProteomicsData__Expression_long$Colname, function (x){
  if (grepl("pre_resting", x)) {
    return("BL.RE")
  } else if (grepl("pre_acuteexercise", x)) {
    return("BL.AC")
  } else if (grepl("post_resting", x)) {
    return("TR.RE")
  } else if (grepl("post_acuteexercise", x)) {
    return("TR.AC")
  }
})

ProteomicsData__Expression_long <- left_join(ProteomicsData__Expression_long, Resp_Group,
                                             by = "MyoID")

##----------------------------------------------------------------------------------------------------####
## Plotting Data -----
Plots = list()

All_Proteins = unique(ProteomicsData__Expression_long$Gene)

T5T0_FC_Col <- grep("^T5T0_FC.FoldChange$", colnames(LimmaAnalysis))
T5T0_pValue_Col <- grep("^T5T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T2T0_FC_Col <- grep("^T2T0_FC.FoldChange$", colnames(LimmaAnalysis))
T2T0_pValue_Col <- grep("^T2T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T5_FC_Col <- grep("^T4T5_FC.FoldChange$", colnames(LimmaAnalysis))
T4T5_pValue_Col <- grep("^T4T5_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T2_FC_Col <- grep("^T4T2_FC.FoldChange$", colnames(LimmaAnalysis))
T4T2_pValue_Col <- grep("^T4T2_FC.FoldChange.pValue", colnames(LimmaAnalysis))


for (Protein in All_Proteins) {
  
  Protein_NR <- grep(paste0("^", Protein, "$"), All_Proteins)
  
  Limma_Row <- grep(paste0("^", Protein, "$"), LimmaAnalysis$Uniprot_ID)
  Gene_Name <- LimmaAnalysis[Limma_Row, "Gene_name"]
  
  Data_for_Plot <- ProteomicsData__Expression_long %>% filter(Gene == Protein) %>% 
    mutate(Timepoint = factor(Timepoint, levels = c("BL.RE", "BL.AC", "TR.RE", "TR.AC")))
  
  
  AnnotationSignif <- data.frame("Start" = c("BL.RE", "TR.RE", "BL.RE", "BL.AC"),
                                 "End" = c("BL.AC", "TR.AC", "TR.RE", "TR.AC"),
                                 "y_height" = NA,
                                 "pValue" = NA,
                                 "FC" = NA,
                                 "Label_chr" = NA)
  
  AnnotationSignif[1, "pValue"] <- round(LimmaAnalysis[Limma_Row, T2T0_pValue_Col], digits = 4)
  AnnotationSignif[2, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T5_pValue_Col], digits = 4)
  AnnotationSignif[3, "pValue"] <- round(LimmaAnalysis[Limma_Row, T5T0_pValue_Col], digits = 4)
  AnnotationSignif[4, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T2_pValue_Col], digits = 4)
  
  AnnotationSignif[1, "FC"] <- round(LimmaAnalysis[Limma_Row, T2T0_FC_Col], digits = 2)
  AnnotationSignif[2, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T5_FC_Col], digits = 2)
  AnnotationSignif[3, "FC"] <- round(LimmaAnalysis[Limma_Row, T5T0_FC_Col], digits = 2)
  AnnotationSignif[4, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T2_FC_Col], digits = 2)
  
  
  AnnotationSignif[1:2, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 1/5
  AnnotationSignif[3, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 2/5
  AnnotationSignif[4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 3/5
  
  
  AnnotationSignif$Label_chr <- sapply(AnnotationSignif$pValue, function(x){
    if (x < 0.001) {
      x <- as.character("<0.001")
    } else {
      x <- as.character(x)
    }
    return(x)
  })
  AnnotationSignif$Label_chr <- paste(AnnotationSignif$Label_chr, "\nFC ", AnnotationSignif$FC, sep = "")
  

  p <- Data_for_Plot %>%  
    ggplot(aes(x = Timepoint,
               y = value)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 3, position = position_jitter(width = 0.2, seed = 1)) +
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") +
    geom_signif(data = AnnotationSignif,
                aes(xmin = Start, xmax = End, y_position = y_height, annotations = Label_chr), 
                manual = T, vjust = 0.6, textsize = 5.5) +
    scale_x_discrete(breaks = c("BL.RE", "BL.AC", "TR.RE", "TR.AC"),
                     labels = c("BL.RE\n(T0)", "BL.AC\n(T2)", "TR.RE\n(T5)", "TR.AC\n(T4)")) +
    labs(title = paste0(Gene_Name, " (", Protein, ")"), y = "Log2[counts]") +
    theme(
      axis.text = element_text(size = 16, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 22, vjust = 2),
      legend.position = "none",
      # legend.text = element_text(size = 16),
      # legend.title = element_text(size = 15, face = "bold"),
      plot.title = element_text(face = "bold", size = 25),
      axis.ticks.length=unit(.25, "cm"),
    )
  
  
  Plots[[Protein]] <- p
  
  if ((Protein_NR %% 50) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Protein_NR, " / ", length(All_Proteins)))
  }
}


Plot_Number = 1
GGplots_List <- list()
Gene_Slide_Nr <- data.frame("Gene" = All_Proteins,
                            "Slide" = NA)

Blank <- Data_for_Plot %>%  
  ggplot(aes(x = Timepoint,
             y = value)) + 
  geom_blank() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  )
Blank


for(i in seq(1, length(Plots), by = 4)){
  
  GGplots_List[[Plot_Number]] <- grid.arrange(if((i) %in% 1:length(Plots)) Plots[[i]] else Blank,
                                              if((i+1) %in% 1:length(Plots)) Plots[[i+1]] else Blank,
                                              if((i+2) %in% 1:length(Plots)) Plots[[i+2]] else Blank,
                                              if((i+3) %in% 1:length(Plots)) Plots[[i+3]] else Blank)
  
  if ((i) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+1) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+1]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+2) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+2]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+3) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+3]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  
  
  if ((Plot_Number %% 25) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Plot_Number, " / ", length(seq(1, length(Plots), by = 4))))
  }
  
  Plot_Number = Plot_Number + 1
}


## Save as PDF ----
Filepictures <- getwd()
PDF_Name <- "02_ProteomicsNRE2_Raw_Data_plots.PDF"

pdf(file = paste(Filepictures, PDF_Name, sep = "/"), width = 15, height = 15)
for(i in 1:length(GGplots_List)){
  grid.arrange(GGplots_List[[i]])
}
dev.off()

write.xlsx(Gene_Slide_Nr, "02_ProteomicsNRE2_Raw_Data_plots.xlsx", overwrite = TRUE)

##----------------------------------------------------------------------------------------------------####
## Plotting Data - RES-LRE Quantile -----

Plots = list()

All_Proteins = unique(ProteomicsData__Expression_long$Gene)

T5T0_FC_Col <- grep("^T5T0_FC.FoldChange$", colnames(LimmaAnalysis))
T5T0_pValue_Col <- grep("^T5T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T2T0_FC_Col <- grep("^T2T0_FC.FoldChange$", colnames(LimmaAnalysis))
T2T0_pValue_Col <- grep("^T2T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T5_FC_Col <- grep("^T4T5_FC.FoldChange$", colnames(LimmaAnalysis))
T4T5_pValue_Col <- grep("^T4T5_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T2_FC_Col <- grep("^T4T2_FC.FoldChange$", colnames(LimmaAnalysis))
T4T2_pValue_Col <- grep("^T4T2_FC.FoldChange.pValue", colnames(LimmaAnalysis))


for (Protein in All_Proteins) {
  
  Protein_NR <- grep(paste0("^", Protein, "$"), All_Proteins)
  
  Limma_Row <- grep(paste0("^", Protein, "$"), LimmaAnalysis$Uniprot_ID)
  Gene_Name <- LimmaAnalysis[Limma_Row, "Gene_name"]
  
  
  Data_for_Plot <- ProteomicsData__Expression_long %>% filter(Gene == Protein) %>% 
    mutate(Timepoint = factor(Timepoint, levels = c("BL.RE", "BL.AC", "TR.RE", "TR.AC")))
  
  
  AnnotationSignif <- data.frame("Start" = c("BL.RE", "TR.RE", "BL.RE", "BL.AC"),
                                 "End" = c("BL.AC", "TR.AC", "TR.RE", "TR.AC"),
                                 "y_height" = NA,
                                 "pValue" = NA,
                                 "FC" = NA,
                                 "Label_chr" = NA)
  
  AnnotationSignif[1, "pValue"] <- round(LimmaAnalysis[Limma_Row, T2T0_pValue_Col], digits = 4)
  AnnotationSignif[2, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T5_pValue_Col], digits = 4)
  AnnotationSignif[3, "pValue"] <- round(LimmaAnalysis[Limma_Row, T5T0_pValue_Col], digits = 4)
  AnnotationSignif[4, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T2_pValue_Col], digits = 4)
  
  AnnotationSignif[1, "FC"] <- round(LimmaAnalysis[Limma_Row, T2T0_FC_Col], digits = 2)
  AnnotationSignif[2, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T5_FC_Col], digits = 2)
  AnnotationSignif[3, "FC"] <- round(LimmaAnalysis[Limma_Row, T5T0_FC_Col], digits = 2)
  AnnotationSignif[4, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T2_FC_Col], digits = 2)
  
  
  AnnotationSignif[1:2, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 1/5
  AnnotationSignif[3, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 2/5
  AnnotationSignif[4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 3/5
  
  
  
  
  AnnotationSignif$Label_chr <- sapply(AnnotationSignif$pValue, function(x){
    if (x < 0.001) {
      x <- as.character("<0.001")
    } else {
      x <- as.character(x)
    }
    return(x)
  })
  AnnotationSignif$Label_chr <- paste(AnnotationSignif$Label_chr, "\nFC ", AnnotationSignif$FC, sep = "")
  
  p <- Data_for_Plot %>%  
    ggplot(aes(x = Timepoint,
               y = value)) + 
    geom_point(aes(color = Group, shape = Group, fill = Group, group = Group), size = 3,
               position = position_dodge(width = 0.25)) +
    geom_signif(data = AnnotationSignif,
                aes(xmin = Start, xmax = End, y_position = y_height, annotations = Label_chr), 
                manual = T, vjust = 0.6, textsize = 5.5) +
    scale_fill_manual(breaks = c("RES", "LRE", "INT"), values = c("#ab273c", "#000088", "#ffffff")) +
    scale_color_manual(breaks = c("RES", "LRE", "INT"), values = c("#ab273c", "#000088", "#000000")) +
    scale_shape_manual(breaks = c("RES", "LRE", "INT"), values = c(21,21,1)) +
    scale_x_discrete(breaks = c("BL.RE", "BL.AC", "TR.RE", "TR.AC"),
                     labels = c("BL.RE\n(T0)", "BL.AC\n(T2)", "TR.RE\n(T5)", "TR.AC\n(T4)")) +
    labs(title = paste0(Gene_Name, " (", Protein, ")"), y = "Log2[counts]") +
    theme(
      axis.text = element_text(size = 16, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 22, vjust = 2),
      legend.position = c(0.08, 0.93),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 15, face = "bold"),
      legend.background=element_rect(fill = alpha("white", 0)),
      plot.title = element_text(face = "bold", size = 25),
      axis.ticks.length=unit(.25, "cm"),
    ) 
  
  Plots[[Protein]] <- p
  
  if ((Protein_NR %% 50) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Protein_NR, " / ", length(All_Proteins)))
  }
}


Plot_Number = 1
GGplots_List <- list()
Gene_Slide_Nr <- data.frame("Gene" = All_Proteins,
                            "Slide" = NA)

Blank <- Data_for_Plot %>%  
  ggplot(aes(x = Timepoint,
             y = value)) + 
  geom_blank() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  )



for(i in seq(1, length(Plots), by = 4)){
  
  GGplots_List[[Plot_Number]] <- grid.arrange(if((i) %in% 1:length(Plots)) Plots[[i]] else Blank,
                                              if((i+1) %in% 1:length(Plots)) Plots[[i+1]] else Blank,
                                              if((i+2) %in% 1:length(Plots)) Plots[[i+2]] else Blank,
                                              if((i+3) %in% 1:length(Plots)) Plots[[i+3]] else Blank)
  
  if ((i) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+1) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+1]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+2) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+2]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  if ((i+3) %in% 1:length(Plots)) {
    Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+3]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  }
  
  
  if ((Plot_Number %% 25) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Plot_Number, " / ", length(seq(1, length(Plots), by = 4))))
  }
  
  Plot_Number = Plot_Number + 1
}


## Save as PDF - RES-LRE Quantile ----
Filepictures <- getwd()
PDF_Name <- "02_ProteomicsNRE2_Raw_Data_plots__ByResponder_Quantile.PDF"

pdf(file = paste(Filepictures, PDF_Name, sep = "/"), width = 15, height = 15)
for(i in 1:length(GGplots_List)){
  grid.arrange(GGplots_List[[i]])
}
dev.off()

# write.xlsx(Gene_Slide_Nr, "02_ProteomicsNRE2_Raw_Data_plots.xlsx", overwrite = TRUE)

##----------------------------------------------------------------------------------------------------####
## Plotting Data - T2T4 and T5T0 -----
Plots = list()

All_Proteins = unique(ProteomicsData__Expression_long$Gene)

T5T0_FC_Col <- grep("^T5T0_FC.FoldChange$", colnames(LimmaAnalysis))
T5T0_pValue_Col <- grep("^T5T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T2T0_FC_Col <- grep("^T2T0_FC.FoldChange$", colnames(LimmaAnalysis))
T2T0_pValue_Col <- grep("^T2T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T5_FC_Col <- grep("^T4T5_FC.FoldChange$", colnames(LimmaAnalysis))
T4T5_pValue_Col <- grep("^T4T5_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T2_FC_Col <- grep("^T4T2_FC.FoldChange$", colnames(LimmaAnalysis))
T4T2_pValue_Col <- grep("^T4T2_FC.FoldChange.pValue", colnames(LimmaAnalysis))


for (Protein in All_Proteins) {
  
  Protein_NR <- grep(paste0("^", Protein, "$"), All_Proteins)
  
  Limma_Row <- grep(paste0("^", Protein, "$"), LimmaAnalysis$Uniprot_ID)
  Gene_Name <- LimmaAnalysis[Limma_Row, "Gene_name"]
  
  
  Data_for_Plot <- ProteomicsData__Expression_long %>% filter(Gene == Protein) %>% 
    mutate(Timepoint = factor(Timepoint, levels = c("BL.RE", "BL.AC", "TR.RE", "TR.AC")))
  
  
  AnnotationSignif <- data.frame("Start" = c("BL.RE", "TR.RE", "BL.RE", "BL.AC"),
                                 "End" = c("BL.AC", "TR.AC", "TR.RE", "TR.AC"),
                                 "y_height" = NA,
                                 "pValue" = NA,
                                 "FC" = NA,
                                 "Label_chr" = NA)
  
  AnnotationSignif[1, "pValue"] <- round(LimmaAnalysis[Limma_Row, T2T0_pValue_Col], digits = 4)
  AnnotationSignif[2, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T5_pValue_Col], digits = 4)
  AnnotationSignif[3, "pValue"] <- round(LimmaAnalysis[Limma_Row, T5T0_pValue_Col], digits = 4)
  AnnotationSignif[4, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T2_pValue_Col], digits = 4)
  
  AnnotationSignif[1, "FC"] <- round(LimmaAnalysis[Limma_Row, T2T0_FC_Col], digits = 2)
  AnnotationSignif[2, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T5_FC_Col], digits = 2)
  AnnotationSignif[3, "FC"] <- round(LimmaAnalysis[Limma_Row, T5T0_FC_Col], digits = 2)
  AnnotationSignif[4, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T2_FC_Col], digits = 2)
  
  
  AnnotationSignif[1:4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 1/5
  # AnnotationSignif[3, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 2/5
  # AnnotationSignif[4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 3/5
  # 
  
  AnnotationSignif$Label_chr <- sapply(AnnotationSignif$pValue, function(x){
    if (x < 0.001) {
      x <- as.character("<0.001")
    } else {
      x <- as.character(x)
    }
    return(x)
  })
  AnnotationSignif$Label_chr <- paste(AnnotationSignif$Label_chr, "\nFC ", AnnotationSignif$FC, sep = "")
  
  for(TP in c("T5T0", "T2T4")){
    p <- Data_for_Plot %>%  
      filter(if(TP == "T5T0") Timepoint %in% c("BL.RE", "TR.RE") else Timepoint %in% c("BL.AC", "TR.AC")) %>% 
      ggplot(aes(x = Timepoint,
                 y = value)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(size = 3, position = position_jitter(width = 0.2, seed = 1)) +
      stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                   width = .75, linetype = "dashed") +
      geom_signif(data = AnnotationSignif %>% 
                    filter(if(TP == "T5T0") Start == "BL.RE" & End == "TR.RE" else Start == "BL.AC" & End == "TR.AC"),
                  aes(xmin = Start, xmax = End, y_position = y_height, annotations = Label_chr), 
                  manual = T, vjust = 0.6, textsize = 5.5) +
      scale_x_discrete(breaks = c("BL.RE", "BL.AC", "TR.RE", "TR.AC"),
                       labels = c("BL.RE\n(T0)", "BL.AC\n(T2)", "TR.RE\n(T5)", "TR.AC\n(T4)")) +
      labs(title = paste0(Gene_Name, " (", Protein, ")"), y = "Log2[counts]") +
      theme(
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, vjust = 2),
        legend.position = "none",
        # legend.text = element_text(size = 16),
        # legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(face = "bold", size = 25),
        axis.ticks.length=unit(.25, "cm"),
      )
    
    p
    
    Plots[[paste(Protein, TP, sep = "__")]] <- p
  }
  
  if ((Protein_NR %% 50) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Protein_NR, " / ", length(All_Proteins)))
  }
}


Plot_Number = 1
GGplots_List <- list()
# Gene_Slide_Nr <- data.frame("Gene" = All_Proteins,
#                             "Slide" = NA)

Blank <- Data_for_Plot %>%  
  ggplot(aes(x = Timepoint,
             y = value)) + 
  geom_blank() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  )
# Blank


for(i in seq(1, length(Plots), by = 8)){
  
  GGplots_List[[Plot_Number]] <- grid.arrange(if((i) %in% 1:length(Plots)) Plots[[i]] else Blank,
                                              if((i+1) %in% 1:length(Plots)) Plots[[i+1]] else Blank,
                                              if((i+2) %in% 1:length(Plots)) Plots[[i+2]] else Blank,
                                              if((i+3) %in% 1:length(Plots)) Plots[[i+3]] else Blank,
                                              if((i+4) %in% 1:length(Plots)) Plots[[i+4]] else Blank,
                                              if((i+5) %in% 1:length(Plots)) Plots[[i+5]] else Blank,
                                              if((i+6) %in% 1:length(Plots)) Plots[[i+6]] else Blank,
                                              if((i+7) %in% 1:length(Plots)) Plots[[i+7]] else Blank,
                                              ncol = 4)
  
  # if ((i) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+1) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+1]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+2) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+2]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+3) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+3]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  
  
  if ((Plot_Number %% 25) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Plot_Number, " / ", length(seq(1, length(Plots), by = 8))))
  }
  
  Plot_Number = Plot_Number + 1
}


## Save as PDF - T2T4 and T5T0  ----
Filepictures <- getwd()
PDF_Name <- "02_ProteomicsNRE2_Raw_Data_plots_TPsplitted.PDF"

pdf(file = paste(Filepictures, PDF_Name, sep = "/"), width = 15, height = 15)
for(i in 1:length(GGplots_List)){
  grid.arrange(GGplots_List[[i]])
}
dev.off()


##----------------------------------------------------------------------------------------------------####
## Plotting Data - T2T4 and T5T0 - with lines -----
Plots = list()

All_Proteins = unique(ProteomicsData__Expression_long$Gene)

T5T0_FC_Col <- grep("^T5T0_FC.FoldChange$", colnames(LimmaAnalysis))
T5T0_pValue_Col <- grep("^T5T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T2T0_FC_Col <- grep("^T2T0_FC.FoldChange$", colnames(LimmaAnalysis))
T2T0_pValue_Col <- grep("^T2T0_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T5_FC_Col <- grep("^T4T5_FC.FoldChange$", colnames(LimmaAnalysis))
T4T5_pValue_Col <- grep("^T4T5_FC.FoldChange.pValue", colnames(LimmaAnalysis))

T4T2_FC_Col <- grep("^T4T2_FC.FoldChange$", colnames(LimmaAnalysis))
T4T2_pValue_Col <- grep("^T4T2_FC.FoldChange.pValue", colnames(LimmaAnalysis))


for (Protein in All_Proteins) {
  
  Protein_NR <- grep(paste0("^", Protein, "$"), All_Proteins)
  
  Limma_Row <- grep(paste0("^", Protein, "$"), LimmaAnalysis$Uniprot_ID)
  Gene_Name <- LimmaAnalysis[Limma_Row, "Gene_name"]
  
  
  Data_for_Plot <- ProteomicsData__Expression_long %>% filter(Gene == Protein) %>% 
    mutate(Timepoint = factor(Timepoint, levels = c("BL.RE", "BL.AC", "TR.RE", "TR.AC")))
  
  
  AnnotationSignif <- data.frame("Start" = c("BL.RE", "TR.RE", "BL.RE", "BL.AC"),
                                 "End" = c("BL.AC", "TR.AC", "TR.RE", "TR.AC"),
                                 "y_height" = NA,
                                 "pValue" = NA,
                                 "FC" = NA,
                                 "Label_chr" = NA)
  
  AnnotationSignif[1, "pValue"] <- round(LimmaAnalysis[Limma_Row, T2T0_pValue_Col], digits = 4)
  AnnotationSignif[2, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T5_pValue_Col], digits = 4)
  AnnotationSignif[3, "pValue"] <- round(LimmaAnalysis[Limma_Row, T5T0_pValue_Col], digits = 4)
  AnnotationSignif[4, "pValue"] <- round(LimmaAnalysis[Limma_Row, T4T2_pValue_Col], digits = 4)
  
  AnnotationSignif[1, "FC"] <- round(LimmaAnalysis[Limma_Row, T2T0_FC_Col], digits = 2)
  AnnotationSignif[2, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T5_FC_Col], digits = 2)
  AnnotationSignif[3, "FC"] <- round(LimmaAnalysis[Limma_Row, T5T0_FC_Col], digits = 2)
  AnnotationSignif[4, "FC"] <- round(LimmaAnalysis[Limma_Row, T4T2_FC_Col], digits = 2)
  
  
  AnnotationSignif[1:4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 1/5
  # AnnotationSignif[3, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 2/5
  # AnnotationSignif[4, "y_height"] <- max(Data_for_Plot$value) + (range(Data_for_Plot$value)[2] - range(Data_for_Plot$value)[1])* 3/5
  # 
  
  AnnotationSignif$Label_chr <- sapply(AnnotationSignif$pValue, function(x){
    if (x < 0.001) {
      x <- as.character("<0.001")
    } else {
      x <- as.character(x)
    }
    return(x)
  })
  AnnotationSignif$Label_chr <- paste(AnnotationSignif$Label_chr, "\nFC ", AnnotationSignif$FC, sep = "")
  
  for(TP in c("T5T0", "T2T4")){
    p <- Data_for_Plot %>%  
      filter(if(TP == "T5T0") Timepoint %in% c("BL.RE", "TR.RE") else Timepoint %in% c("BL.AC", "TR.AC")) %>% 
      ggplot(aes(x = Timepoint,
                 y = value)) + 
      geom_point(size = 3) +
      geom_line(aes(group = MyoID)) +
      # stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
      #              width = .75, linetype = "dashed") +
      geom_signif(data = AnnotationSignif %>% 
                    filter(if(TP == "T5T0") Start == "BL.RE" & End == "TR.RE" else Start == "BL.AC" & End == "TR.AC"),
                  aes(xmin = Start, xmax = End, y_position = y_height, annotations = Label_chr), 
                  manual = T, vjust = 0.6, textsize = 5.5) +
      scale_x_discrete(breaks = c("BL.RE", "BL.AC", "TR.RE", "TR.AC"),
                       labels = c("BL.RE\n(T0)", "BL.AC\n(T2)", "TR.RE\n(T5)", "TR.AC\n(T4)")) +
      labs(title = paste0(Gene_Name, " (", Protein, ")"), y = "Log2[counts]") +
      theme(
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, vjust = 2),
        legend.position = "none",
        # legend.text = element_text(size = 16),
        # legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(face = "bold", size = 25),
        axis.ticks.length=unit(.25, "cm"),
      )
    
    p
    
    Plots[[paste(Protein, TP, sep = "__")]] <- p
  }
  
  if ((Protein_NR %% 50) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Protein_NR, " / ", length(All_Proteins)))
  }
}


Plot_Number = 1
GGplots_List <- list()
# Gene_Slide_Nr <- data.frame("Gene" = All_Proteins,
#                             "Slide" = NA)

Blank <- Data_for_Plot %>%  
  ggplot(aes(x = Timepoint,
             y = value)) + 
  geom_blank() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  )
# Blank


for(i in seq(1, length(Plots), by = 8)){
  
  GGplots_List[[Plot_Number]] <- grid.arrange(if((i) %in% 1:length(Plots)) Plots[[i]] else Blank,
                                              if((i+1) %in% 1:length(Plots)) Plots[[i+1]] else Blank,
                                              if((i+2) %in% 1:length(Plots)) Plots[[i+2]] else Blank,
                                              if((i+3) %in% 1:length(Plots)) Plots[[i+3]] else Blank,
                                              if((i+4) %in% 1:length(Plots)) Plots[[i+4]] else Blank,
                                              if((i+5) %in% 1:length(Plots)) Plots[[i+5]] else Blank,
                                              if((i+6) %in% 1:length(Plots)) Plots[[i+6]] else Blank,
                                              if((i+7) %in% 1:length(Plots)) Plots[[i+7]] else Blank,
                                              ncol = 4)
  
  # if ((i) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+1) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+1]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+2) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+2]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  # if ((i+3) %in% 1:length(Plots)) {
  #   Gene_Slide_Nr[grep(paste0("^", names(Plots)[[i+3]], "$"), Gene_Slide_Nr$Gene), "Slide"] <- Plot_Number
  # }
  
  
  if ((Plot_Number %% 25) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Plot_Number, " / ", length(seq(1, length(Plots), by = 8))))
  }
  
  Plot_Number = Plot_Number + 1
}


## Save as PDF - T2T4 and T5T0 -with lines ----
Filepictures <- getwd()
PDF_Name <- "02_ProteomicsNRE2_Raw_Data_plots_TPsplitted_WithLines.PDF"

pdf(file = paste(Filepictures, PDF_Name, sep = "/"), width = 15, height = 15)
for(i in 1:length(GGplots_List)){
  grid.arrange(GGplots_List[[i]])
}
dev.off()

