library(openxlsx)

load("05_Korrelation_ProteomicsNRE2_mit_KlinischenDaten.RData")
T4T2 <- KorrelationTable_Results_List[["T4T2"]]

# from new proteomics data set "Proteomics_Protein_Info_GO"
Mito_CI <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB1", "NDUFB10", "NDUFB11", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2", "NDUFC2-KCTD14", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3")
Mito_CII <- c("SDHC", "SDHB", "SDHD", "SDHA")
Mito_CIII <- c("BCS1L", "BCS1L", "BRAWNIN", "CYC1", "LYRM7", "MT-CO1", "MT-CYB", "UQCC1", "UQCC2", "UQCR10", "UQCR11", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRFS1", "UQCRH", "UQCRHL", "UQCRQ")
Mito_CIV <- c("COX4I1", "COX5A", "COX5B", "COX6A2", "COX6B1", "COX6C", "COX7A1", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "MT-CO1", "MT-CO2", "MT-CO3", "NDUFA4", "UQCRC2", "UQCRFS1")
Mito_CV <- c("ATP5F1A", "ATP5F1A", "ATP5F1A", "ATP5F1B", "ATP5F1B", "ATP5F1C", "ATP5F1C", "ATP5F1D", "ATP5F1D", "ATP5F1D", "ATP5F1E", "ATP5F1E", "ATP5ME", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MG", "ATP5MJ", "ATP5MK", "ATP5PB", "ATP5PB", "ATP5PD", "ATP5PD", "ATP5PD", "ATP5PF", "ATP5PF", "ATP5PO", "ATP5PO", "ATPAF1", "ATPAF2", "FMC1", "MT-ATP6", "MT-ATP8", "MT-ATP8", "PPIF")
Mito_TCA <- c("IDH1", "PDHA2", "ACO1", "SDHC", "MDH1", "MDH2", "IDH3B", "OGDH", "PDHB", "IDH3G", "DLAT", "OGDHL", "FH", "PDHA1", "DLST", "SUCLA2", "SUCLG2", "IDH2", "SDHB", "DHTKD1", "SUCLG1", "SDHD", "NNT", "ACO2", "SDHA", "MRPS36", "IDH3A", "CS")

# from old Proteomics data set
Mito_QLink <- c("CHDH", "DHODH", "ETFA", "ETFB", "ETFBKMT", "ETFDH", "ETFRF1", "GPD2", "PRODH", "PRODH2", "SQOR")


table(T4T2$Gene_name %in% unique(Mito_CI))
table(T4T2$Gene_name %in% unique(Mito_CII))
table(T4T2$Gene_name %in% unique(Mito_CIII))
table(T4T2$Gene_name %in% unique(Mito_CIV))
table(T4T2$Gene_name %in% unique(Mito_CV))
table(T4T2$Gene_name %in% unique(Mito_TCA))
table(T4T2$Gene_name %in% unique(Mito_QLink))

## grep all Complexes
Mitolist <- sapply(ls(pattern = "Mito_"), get)

## select List-element with most values
MaxValue <- max(lengths(Mitolist))

## set length of each list-element to max number, filter genes present in data table and put it back into the list
for(i in names(Mitolist)){
  
  MitoElement = unique(sort(Mitolist[[i]]))
  MitoElement = MitoElement[MitoElement %in% T4T2$Gene_name]
  length(MitoElement) = MaxValue
  Mitolist[[i]] = MitoElement
}

MitoDataframe <- do.call(cbind.data.frame, Mitolist)

write.xlsx(MitoDataframe, file = "11_Annotating_Gene_to_Mitochondria_Part2.xlsx")
