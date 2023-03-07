library(tidyverse)
library(ggsignif)
library(gridExtra)
# library(doParallel)
theme_set(theme_classic())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load Data
load("03_Limma_Analyse__Split_RESvsLRE__FC.RData")

TPs <- grep("FC_log2$", names(final_data), value = TRUE)

Plots = list()


  
for (Protein in final_data$Uniprot_ID) {
  
  for (TP_FC_log2 in TPs) {
    
    TP_FC = str_remove(TP_FC_log2, "_log2")
    # TP = str_remove(TP_FC, "_FC")
    
    final_data_short <- final_data %>% select(Gene_name, Uniprot_ID, contains(TP_FC))
    Gene_Name = final_data_short[which(Protein == final_data_short$Uniprot_ID), "Gene_name"]
    
    Table_for_plot <- final_data_short %>% 
      select(Gene_name, Uniprot_ID, contains(paste0("__", TP_FC))) %>% 
      filter(Uniprot_ID == Protein) %>% 
      pivot_longer(., cols = contains(paste0("__", TP_FC)),
                   names_to = c("Subject", "RespGroup"),
                   names_sep = "__",
                   values_to = "value") %>% 
      mutate(RespGroup = str_remove(RespGroup, paste0(TP_FC, "_")))
    
    
    
    AnnotationSignif <- data.frame("Start" = c("LRE"),
                                   "End" = c("RES"),
                                   "y_height" = NA,
                                   "pValue" = NA,
                                   "FC" = NA,
                                   "Label_chr" = NA)
    
    AnnotationSignif[1, "pValue"] <- round(final_data_short[which(Protein == final_data_short$Uniprot_ID), 
                                                            paste0(TP_FC, "___FoldChange___pValue")], 
                                           digits = 4)
    AnnotationSignif[1, "FC"] <- round(final_data_short[which(Protein == final_data_short$Uniprot_ID), 
                                                        paste0(TP_FC, "___FoldChange")], 
                                       digits = 2)
    
    AnnotationSignif[1, "y_height"] <- max(Table_for_plot$value) + (range(Table_for_plot$value)[2] - range(Table_for_plot$value)[1])* 1/5
    
    AnnotationSignif$Label_chr <- sapply(AnnotationSignif$pValue, function(x){
      if (x < 0.001) {
        x <- as.character("<0.001")
      } else {
        x <- as.character(x)
      }
      return(x)
    })
    AnnotationSignif$Label_chr <- paste(AnnotationSignif$Label_chr, "\nFC ", AnnotationSignif$FC, sep = "")
    
    Plots[[paste(Protein, TP_FC, sep = "__")]] <-
      Table_for_plot %>% 
      ggplot(aes(x = RespGroup,
                 y = value)) +
      geom_boxplot(aes(fill = RespGroup), alpha = 0.3, outlier.shape = NA) +
      geom_point(aes(color = RespGroup), size = 3, position = position_jitter(width = 0.2, seed = 1)) +
      stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                   width = .75, linetype = "dashed") +
      geom_signif(data = AnnotationSignif,
                  aes(xmin = Start, xmax = End, y_position = y_height, annotations = Label_chr), 
                  manual = T, vjust = 0.6, textsize = 5.5) +
      scale_fill_manual(breaks = c("RES", "LRE"), values = c("#0072B2", "#D55E00")) +
      scale_color_manual(breaks = c("RES", "LRE"), values = c("#0072B2", "#D55E00")) +
      labs(title = paste0(Gene_Name, " (", Protein, ")"), 
           subtitle = TP_FC,
           y = "Log2[counts]") +
      theme(
        axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22, vjust = 2),
        legend.position = "none",
        # legend.text = element_text(size = 16),
        # legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(face = "bold", size = 25),
        plot.subtitle = element_text(size = 18),
        axis.ticks.length=unit(.25, "cm"),
      ) +
      coord_cartesian(clip = "off")
    
    if (length(Plots) %% 500 == 0) {
      print(paste0(format(Sys.time(), "%x %T"), " Finished: ", 
                   length(Plots), "/", 
                   length(final_data_short$Uniprot_ID)*length(TPs)))
      # break
    }
  }
  
}


Blank <- Table_for_plot %>% 
  ggplot(aes(x = RespGroup,
             V = value)) +
  geom_blank() +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )

# Plot_Number = 1
GGplots_List <- list()
# Gene_Slide_Nr <- data.frame("Gene" = All_Proteins,
#                             "Slide" = NA)


# no_cores <- detectCores() - 1
# no_cores <- 2
# registerDoParallel(cores=no_cores)
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)

for(i in seq(1, length(Plots), by = 12)) {
  
  Plot_Number = which(i == seq(1, length(Plots), by = 12))
  
  GGplots_List[[Plot_Number]] <- grid.arrange(if((i) %in% 1:length(Plots)) Plots[[i]] else Blank,
                                              if((i+1) %in% 1:length(Plots)) Plots[[i+1]] else Blank,
                                              if((i+2) %in% 1:length(Plots)) Plots[[i+2]] else Blank,
                                              if((i+3) %in% 1:length(Plots)) Plots[[i+3]] else Blank,
                                              if((i+4) %in% 1:length(Plots)) Plots[[i+4]] else Blank,
                                              if((i+5) %in% 1:length(Plots)) Plots[[i+5]] else Blank,
                                              if((i+6) %in% 1:length(Plots)) Plots[[i+6]] else Blank,
                                              if((i+7) %in% 1:length(Plots)) Plots[[i+7]] else Blank,
                                              if((i+8) %in% 1:length(Plots)) Plots[[i+8]] else Blank,
                                              if((i+9) %in% 1:length(Plots)) Plots[[i+9]] else Blank,
                                              if((i+10) %in% 1:length(Plots)) Plots[[i+10]] else Blank,
                                              if((i+11) %in% 1:length(Plots)) Plots[[i+11]] else Blank,
                                              ncol = 6)
  
  
  
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
  
  
  if ((Plot_Number %% 50) == 0) {
    print(paste0(format(Sys.time(), "%x %T"), "  Finished: ", Plot_Number, " / ", length(seq(1, length(Plots), by = 12))))
  }
  # break
  
  # Plot_Number = Plot_Number + 1
}

## Save as PDF - T2T4 and T5T0  ----
Filepictures <- getwd()
PDF_Name <- "19_Limma_Split_by_FC_plots.PDF"

pdf(file = paste(Filepictures, PDF_Name, sep = "/"), width = 25, height = 13)
for(i in 1:length(GGplots_List)){
  grid.arrange(GGplots_List[[i]])
  
}
dev.off()

## Save RData
save(Plots, file = "19_Limma_Split_by_FC_plots__singlePlots.RData")

### Save Single Plot
final_data %>% filter(Gene_name == "AKR1C3") 

Needed_Plot = Plots[["P42330__T5T0_FC"]]
plot(Needed_Plot)
