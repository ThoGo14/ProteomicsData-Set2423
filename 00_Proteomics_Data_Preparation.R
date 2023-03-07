#---------------------------------------------------------------------------------------------------------#
#
#---------------------------------------------------------------------------------------------------------#
httr::set_config(httr::config(ssl_verifypeer = FALSE))
# library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
library(biomaRt)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ProteomicsData__Raw <- read.xlsx("220426-NRE-proteome-PD_incl_fractions-to-Thomas.xlsx", sheet = "data_ForR") %>% 
  mutate(KEEP = TRUE)
ProteomicsData__Metadata <- read.xlsx("220426-NRE-proteome-PD_incl_fractions-to-Thomas.xlsx", sheet = "Subjects")

## Log2 transformation
ProteomicsData__Raw.Log2 <- ProteomicsData__Raw %>% 
  mutate_at(vars(contains("Abundances")), .funs = log2)


## Data imputation ----
impute_data = function(df, width, downshift){
  
  # Create new columns indicating whether the values are imputed
  Log2.names = grep("^Abundances", names(df), value = TRUE)
  impute.names = paste0("Imputed__", Log2.names)
  
  df[impute.names] = lapply(Log2.names, function(x) !is.finite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[Log2.names] = lapply(Log2.names,
                          function(x) {
                            temp = df[[x]]
                            temp[!is.finite(temp)] = NA
                            temp.sd = width * sd(temp[df$KEEP], na.rm = TRUE)   # shrink sd width
                            temp.mean = mean(temp[df$KEEP], na.rm = TRUE) - 
                              downshift * sd(temp[df$KEEP], na.rm = TRUE)   # shift mean of imputed values
                            n.missing = sum(is.na(temp))
                            temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                            return(temp)
                          })
  width <<- width
  downshift <<- downshift
  return(df)
}

ProteomicsData__Raw.Log2.Impute = impute_data(ProteomicsData__Raw.Log2,
                                              width = 0.5,
                                              downshift = 1.5)

## Plot graph for imputed Data
Proteomics_imputed_plot <- function(df) { 
  Proteomics1 <- df %>% 
    dplyr::select(starts_with("Abundances")) %>% 
    pivot_longer(cols = starts_with("Abundances"),
                 names_to = "Names",
                 values_to = "Values")
  
  Proteomics2 <- df %>% 
    dplyr::select(starts_with("Imputed")) %>% 
    pivot_longer(cols = starts_with("Imputed"),
                 names_to = "Names_imputed",
                 values_to = "Imputed")
  
  Proteomics3 <- cbind(Proteomics1, Proteomics2)
  Proteomics3 <- Proteomics3 %>% 
    separate(Names, c("Abundances", "F_Nr", "SampleID", "Timepoint"), sep = "__") %>% 
    mutate(SampleID = str_replace(SampleID, "Sample_", "S_"))
  
  
  Plot <- Proteomics3 %>% 
    ggplot() +
    geom_histogram(data = . %>% filter(Imputed == FALSE), 
                   aes(x = Values, fill = Imputed), alpha = 0.5, bins = 40) +
    geom_histogram(data = . %>% filter(Imputed == TRUE), 
                   aes(x = Values, fill = Imputed), alpha = 0.5, bins = 40) +
    facet_grid(SampleID ~ Timepoint) +
    ggtitle(paste0("width = ", width, ", downshift = ", downshift)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(Plot, filename = paste0("00_Data Preparation/",
                                 "Proteomics_imputed__w", width, "__d", downshift, ".png"),
         width = 20, height = 30, units = "cm", dpi = 200)
  
  print(Plot)
}

Proteomics_imputed_plot(ProteomicsData__Raw.Log2.Impute)
## Adding Missing Gene_names ----
UniprotData <- read.delim("E:/Auswertung mit R/00_Gene_Annotation/Uniprot/HUMAN_9606_idmapping.dat", header = F)
UniprotData <- rename(UniprotData, Uniprot.ID = V1, Attribute = V2,  Value = V3)

Proteomics_UniprotIDs <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID

UniprotData_selected <- UniprotData %>% 
  filter(Uniprot.ID %in% Proteomics_UniprotIDs, Attribute == "Gene_Name") %>% 
  distinct(Uniprot.ID, .keep_all = T) %>% 
  dplyr::select(Uniprot.ID, Value) %>% 
  rename(Gene_name = Value,
         Uniprot_ID = Uniprot.ID) 

ProteomicsData__Raw.Log2.Impute <- left_join(ProteomicsData__Raw.Log2.Impute, 
                                             UniprotData_selected,
                                             by = "Uniprot_ID")
ProteomicsData__Raw.Log2.Impute <- relocate(ProteomicsData__Raw.Log2.Impute,
                                            Gene_name,
                                            .before = Uniprot_ID)
## Not found Uniprot_ID
ProteomicsData__Raw.Log2.Impute$Uniprot_ID[is.na(ProteomicsData__Raw.Log2.Impute$Gene_name)]
# [1] "P02769"  "P0DOX7"  "Q56UQ5"  "His6-FC" "P0DOX8" 
# P02769 == ALBU_BOVIN
# P0DOX7 == IGK_HUMAN (no Gene name) --> IGKC?
# Q56UQ5 == TPT1L_HUMAN --> LOC121627959
# His6-FC
# P0DOX8 == IGL1_HUMAN --> ??

## Remove contaminations ----
ProteomicsData__Raw.Log2.Impute <- ProteomicsData__Raw.Log2.Impute[-grep("P02769|His6-FC",
                                                                         ProteomicsData__Raw.Log2.Impute$Uniprot_ID), ]

## Change Colnames of Data tables ----
### Imputed Data
ProteomicsData__Metadata$Old_Colname = grep("^Abundance",names(ProteomicsData__Raw.Log2.Impute), value = T)
ProteomicsData__Metadata$New_Colname = paste(ProteomicsData__Metadata$MyoID,
                                             ProteomicsData__Metadata$Timepoint,
                                             sep = "__")

names(ProteomicsData__Raw.Log2.Impute)[grep("^Abundance",
                                            names(ProteomicsData__Raw.Log2.Impute))] <- paste(ProteomicsData__Metadata$MyoID,
                                                                                              ProteomicsData__Metadata$Timepoint,
                                                                                              sep = "__")
names(ProteomicsData__Raw.Log2.Impute)[grep("^Imputed",
                                            names(ProteomicsData__Raw.Log2.Impute))] <- paste("Imputed",
                                                                                              ProteomicsData__Metadata$MyoID,
                                                                                              ProteomicsData__Metadata$Timepoint,
                                                                                              sep = "__")
### Non Imputed Data
names(ProteomicsData__Raw.Log2)[grep("^Abundance",
                                     names(ProteomicsData__Raw.Log2))] <- paste(ProteomicsData__Metadata$MyoID,
                                                                                ProteomicsData__Metadata$Timepoint,
                                                                                sep = "__")

### Raw Data
names(ProteomicsData__Raw)[grep("^Abundance",
                                names(ProteomicsData__Raw))] <- paste(ProteomicsData__Metadata$MyoID,
                                                                      ProteomicsData__Metadata$Timepoint,
                                                                      sep = "__")
## Add Rating to Metadata ----
NRE2_Subject_Rating <- read.delim("NRE2_Subject_Rating.txt", header = TRUE,
                                  sep = "\t")
NRE2_Subject_Rating <- NRE2_Subject_Rating %>% 
  mutate(subjectID = as.character(subjectID),
         MyoID = paste0("M", MyoID))

ProteomicsData__Metadata <- left_join(ProteomicsData__Metadata,
                                      NRE2_Subject_Rating,
                                      by = "MyoID")




## Calculating Training FC, preAcute FC and postAcute FC ----
## Data is log-transformed, therefore substraction and not division
table(ProteomicsData__Metadata$Timepoint)

Subjects_preRest <- ProteomicsData__Metadata[grep("pre_resting", ProteomicsData__Metadata$Timepoint), "New_Colname"]
names(Subjects_preRest) <- ProteomicsData__Metadata[grep("pre_resting", ProteomicsData__Metadata$Timepoint), "subjectID"]
Subjects_preAcute <- ProteomicsData__Metadata[grep("pre_acute", ProteomicsData__Metadata$Timepoint), "New_Colname"]
names(Subjects_preAcute) <- ProteomicsData__Metadata[grep("pre_acute", ProteomicsData__Metadata$Timepoint), "subjectID"]
Subjects_postRest <- ProteomicsData__Metadata[grep("post_resting", ProteomicsData__Metadata$Timepoint), "New_Colname"]
names(Subjects_postRest) <- ProteomicsData__Metadata[grep("post_resting", ProteomicsData__Metadata$Timepoint), "subjectID"]
Subjects_postAcute <- ProteomicsData__Metadata[grep("post_acute", ProteomicsData__Metadata$Timepoint), "New_Colname"]
names(Subjects_postAcute) <- ProteomicsData__Metadata[grep("post_acute", ProteomicsData__Metadata$Timepoint), "subjectID"]


ProteomicsData__T2T0 <- ProteomicsData__Raw.Log2.Impute[, Subjects_preAcute] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_preRest]
colnames(ProteomicsData__T2T0) <- sapply(colnames(ProteomicsData__T2T0), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T2T0_FC")
  return(x)
})
row.names(ProteomicsData__T2T0) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


ProteomicsData__T5T0 <- ProteomicsData__Raw.Log2.Impute[, Subjects_postRest] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_preRest]
colnames(ProteomicsData__T5T0) <- sapply(colnames(ProteomicsData__T5T0), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T5T0_FC")
  return(x)
})
row.names(ProteomicsData__T5T0) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


ProteomicsData__T4T5 <- ProteomicsData__Raw.Log2.Impute[, Subjects_postAcute] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_postRest[intersect(names(Subjects_postAcute), names(Subjects_postRest))]]
colnames(ProteomicsData__T4T5) <- sapply(colnames(ProteomicsData__T4T5), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T4T5_FC")
  return(x)
})
row.names(ProteomicsData__T4T5) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


ProteomicsData__T4T2 <- ProteomicsData__Raw.Log2.Impute[, Subjects_postAcute] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_preAcute[intersect(names(Subjects_postAcute), names(Subjects_preAcute))]]
colnames(ProteomicsData__T4T2) <- sapply(colnames(ProteomicsData__T4T2), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T4T2_FC")
  return(x)
})
row.names(ProteomicsData__T4T2) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


ProteomicsData__T5T2 <- ProteomicsData__Raw.Log2.Impute[, Subjects_postRest] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_preAcute]
colnames(ProteomicsData__T5T2) <- sapply(colnames(ProteomicsData__T5T2), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T5T2_FC")
  return(x)
})
row.names(ProteomicsData__T5T2) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


ProteomicsData__T4T0 <- ProteomicsData__Raw.Log2.Impute[, Subjects_postAcute] - 
  ProteomicsData__Raw.Log2.Impute[, Subjects_preRest[intersect(names(Subjects_postAcute), names(Subjects_preRest))]]
colnames(ProteomicsData__T4T0) <- sapply(colnames(ProteomicsData__T4T0), function(x){
  x <- unlist(str_split(x, "_"))[1]
  x <- paste0(x, "__T4T0_FC")
  return(x)
})
row.names(ProteomicsData__T4T0) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID


## Protein Information list ----
Proteomics_Protein_Info <- ProteomicsData__Raw.Log2.Impute[, 1:6]

## Adding more IDs into new tables ----
Proteomics_Protein_Info_IDs <- Proteomics_Protein_Info

## ensembl_gene_id
hsa.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
hsa.attributes <- listAttributes(hsa.mart)


hsa.annotation <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id"),
                        mart = hsa.mart)

table(Proteomics_Protein_Info$Gene_name %in% hsa.annotation$hgnc_symbol)

hsa.annotation <- rename(hsa.annotation, Gene_name = hgnc_symbol)

Proteomics_Protein_Info_IDs <- left_join(Proteomics_Protein_Info_IDs, 
                                         hsa.annotation,
                                         by = "Gene_name")

## refseq_mrna

hsa.annotation <- getBM(attributes = c("hgnc_symbol", "refseq_mrna"),
                        mart = hsa.mart)

table(Proteomics_Protein_Info$Gene_name %in% hsa.annotation$hgnc_symbol[hsa.annotation$refseq_mrna != ""])

hsa.annotation <- rename(hsa.annotation, Gene_name = hgnc_symbol)

Proteomics_Protein_Info_IDs <- left_join(Proteomics_Protein_Info_IDs, 
                                         hsa.annotation,
                                         by = "Gene_name")
Proteomics_Protein_Info_IDs$refseq_mrna[Proteomics_Protein_Info_IDs$refseq_mrna == ""] <- NA
	

# Proteomics_Protein_Info_IDs <- Proteomics_Protein_Info_IDs[complete.cases(Proteomics_Protein_Info_IDs), ]

## New Table: go_id (GO:00...), name_1006 (Term), definition_1006 (definition), namespace_1003 (domain, BP/CC/MF)
## go_linkage_type (Evidence Code)
hsa.annotation <- getBM(attributes = c("uniprot_gn_id", "go_id", "name_1006", "definition_1006", "namespace_1003", "go_linkage_type"),
                        filters = "uniprot_gn_id",
                        values = Proteomics_Protein_Info$Uniprot_ID,
                        mart = hsa.mart)

table(Proteomics_Protein_Info$Uniprot_ID %in% hsa.annotation$uniprot_gn_id)

hsa.annotation <- rename(hsa.annotation, Uniprot_ID = uniprot_gn_id)

Proteomics_Protein_Info_GO <- right_join(Proteomics_Protein_Info, 
                                         hsa.annotation,
                                         by = "Uniprot_ID")

Proteomics_Protein_Info_GO <- rename(Proteomics_Protein_Info_GO,
                                     GO.ID = go_id,
                                     GO.Term = name_1006,
                                     GO.Definition = definition_1006,
                                     GO.Domain = namespace_1003,
                                     GO.Evidence_code = go_linkage_type)

Proteomics_Protein_Info_GO[Proteomics_Protein_Info_GO == ""] <- NA
# Proteomics_Protein_Info_GO <- Proteomics_Protein_Info_GO[complete.cases(Proteomics_Protein_Info_GO), ]




## Remove Meta Data from Data tables ----
row.names(ProteomicsData__Raw) = ProteomicsData__Raw$Uniprot_ID
ProteomicsData__Raw <- ProteomicsData__Raw[, -c(1:5)]

row.names(ProteomicsData__Raw.Log2) <- ProteomicsData__Raw.Log2$Uniprot_ID
ProteomicsData__Raw.Log2 <- ProteomicsData__Raw.Log2[, -c(1:5)]

row.names(ProteomicsData__Raw.Log2.Impute) <- ProteomicsData__Raw.Log2.Impute$Uniprot_ID
ProteomicsData__Raw.Log2.Impute <- ProteomicsData__Raw.Log2.Impute[, -c(1:6)]
# remove which col are imputed (True/False)
ProteomicsData__Raw.Log2.ImputeCols <- ProteomicsData__Raw.Log2.Impute[, grep("^Imputed",
                                                                              names(ProteomicsData__Raw.Log2.Impute))]
ProteomicsData__Raw.Log2.Impute <- ProteomicsData__Raw.Log2.Impute[, -grep("^Imputed",
                                                                           names(ProteomicsData__Raw.Log2.Impute))]
ProteomicsData__Raw.Log2.Impute$KEEP <- NULL

## Save RData of Expression -----
ProteomicsDataNRE2 <- list("ProteomicsData__Raw" = ProteomicsData__Raw,
                           "ProteomicsData__Raw.Log2" = ProteomicsData__Raw.Log2,
                           "ProteomicsData__Raw.Log2.Impute" = ProteomicsData__Raw.Log2.Impute,
                           "ProteomicsData__Metadata" = ProteomicsData__Metadata,
                           "ProteomicsData__Raw.Log2.ImputeCols" = ProteomicsData__Raw.Log2.ImputeCols,
                           "ProteomicsData__T4T5" = ProteomicsData__T4T5,
                           "ProteomicsData__T2T0" = ProteomicsData__T2T0,
                           "ProteomicsData__T5T0" = ProteomicsData__T5T0,
                           "ProteomicsData__T4T2" = ProteomicsData__T4T2,
                           "ProteomicsData__T5T2" = ProteomicsData__T5T2,
                           "ProteomicsData__T4T0" = ProteomicsData__T4T0,
                           "Proteomics_Protein_Info" = Proteomics_Protein_Info,
                           "Proteomics_Protein_Info_IDs" = Proteomics_Protein_Info_IDs,
                           "Proteomics_Protein_Info_GO" = Proteomics_Protein_Info_GO)

save(ProteomicsDataNRE2,
     file = "ProteomicsNRE2__Expression.RData")


