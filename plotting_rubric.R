# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R

install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('ggridges')
install.packages('ggrepel')
install.packages('dplyr')
#install.packages('Seurat')
install.packages('plotly')
install.packages('clustree')
install.packages('patchwork')
install.packages('future')
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
BiocManager::install("EnhancedVolcano")
BiocManager::install("DoubletFinder")
BiocManager::install("glmGamPoi")
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("GenomeInfoDb")
BiocManager::install("MeSHDbi")
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("DOSE")
BiocManager::install("dittoSeq")
BiocManager::install("escape")
BiocManager::install("ComplexHeatmap")
BiocManager::install(c("DropletUtils", "Nebulosa"))
BiocManager::install("hdf5r", force = TRUE)
BiocManager::install('multtest')
BiocManager::install("MAST")
BiocManager::install("enrichplot")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
BiocManager::install("EnrichmentBrowser")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("devtools")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")
install.packages('metap')
install.packages('magick')
install.packages("corrplot")
install.packages("RColorBrewer")
install.packages("sunburstR")
install_github("didacs/ggsunburst")
install.packages("ggpubr")


# Run the following code once you have Seurat installed
suppressWarnings(
  {
    library(leiden)
    library(stringr)
    library(hdf5r)
    library(SoupX)
    library(Rcpp)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(purrr)
    library(tidyverse)
    library(data.table)
    library(reticulate)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(Signac)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    library(glmGamPoi)
    library(GOSemSim)
    library(org.Hs.eg.db)
    library(AnnotationHub)
    library(MeSHDbi)
    library(clusterProfiler)
    library(DOSE)
    library(dittoSeq)
    library(escape)
    library(EnrichmentBrowser)
    library(viridisLite)
    library(viridis)
    library(ComplexHeatmap)
    library(circlize)
    #library(scCustomize)
    library(Nebulosa)
    library(DropletUtils)
    library(ggvenn)
    library(ggVennDiagram)
    library(devtools)
    library(R.utils)
    library(qs)
    library(multtest)
    library(metap)
    library(MAST)
    library(magick)
    library(enrichplot)
    library(corrplot)
    library(DESeq2)
    library(RColorBrewer)
    library(sunburstR)
    library(d3r)
    library(ggpubr)
  }
)


# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())
session_info()
sessionInfo()
# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("Seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")
packageVersion("org.Hs.eg.db")
packageVersion("DESeq2")

############################ STAGE  ############################
############################  12a   ############################
#Load dataset
#processed_rna <- qread(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\processed_rna.qs)")
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
#saveRDS(processed_rna, file = r"(E:\2.SexbasedStudyCurrent\RDS files\processed_rna.rds)")

# Organise metadata
# Add metadata
Idents(processed_rna) <- "ancestry"
processed_rna$ancestry_sex <- paste(Idents(processed_rna), processed_rna$'Sex', sep = "_")
table(processed_rna$ancestry_sex)

Idents(processed_rna) <- "celltype_sex_ancestry_disease"
processed_rna$celltype_sex_ancestry_disease_lib <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib)

Idents(processed_rna) <- "celltype_sex_ancestry_disease_lib"
processed_rna$celltype_sex_ancestry_disease_lib_source <- paste(Idents(processed_rna), processed_rna$'tissue_source', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib_source)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', sep = "_")
table(processed_rna$disease_ancestry_lib_sex)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex_source <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'tissue_source', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex_source_celltype <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'tissue_source', processed_rna$'celltype_qadir', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source_celltype)

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M", "ND_black_HP2110001_M", "ND_black_HP2123201_M", "ND_black_HPAP-052_M", "ND_black_HPAP-080_M", #Black M ND
               "ND_black_HP2106201_F", "ND_black_HP2121601_F", "ND_black_HP2132801_F", "ND_black_HP2202101_F", #Black F ND
               
               #Hispanic M ND
               "ND_hispanic_HPAP-099_F", "ND_hispanic_HPAP-101_F", "ND_hispanic_HPAP-105_F", #Hispanic F ND
               
               "ND_white_HP2107001_M", "ND_white_HP2107901_M", "ND_white_HPAP-026_M", "ND_white_HPAP-035_M", "ND_white_HPAP-040_M", "ND_white_HPAP-056_M", "ND_white_HPAP-059_M", "ND_white_HPAP-075_M", "ND_white_HPAP-077_M", "ND_white_HPAP-082_M", "ND_white_SAMN15877725_M", #White M ND
               "ND_white_HP2022801_F", "ND_white_HP2024001_F", "ND_white_HP2105501_F", "ND_white_HP2108601_F", "ND_white_HP2108901_F", "ND_white_HPAP-022_F", "ND_white_HPAP-036_F", "ND_white_HPAP-037_F", "ND_white_HPAP-053_F", "ND_white_HPAP-054_F",  "ND_white_HPAP-063_F", "ND_white_HPAP-074_F", "ND_white_HPAP-103_F", #White F ND  
               
               "T2D_black_HPAP-065_M", "T2D_black_HPAP-070_M", "T2D_black_HPAP-083_M", "T2D_black_HPAP-108_M", #Black M T2D  
               "T2D_black_HPAP-051_F", "T2D_black_HPAP-058_F", "T2D_black_HPAP-061_F", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F", "T2D_hispanic_HPAP-091_F", "T2D_hispanic_HPAP-109_F", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M", "T2D_white_HPAP-100_M", "T2D_white_HPAP-106_M",# White M T2D
               "T2D_white_HPAP-057_F", "T2D_white_HPAP-081_F", "T2D_white_HPAP-085_F") # White F T2D

#Check
table(processed_rna$disease_ancestry_lib_sex)
summary(is.na(table(processed_rna$disease_ancestry_lib_sex)))
is.na(table(processed_rna$disease_ancestry_lib_sex))

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex <- factor(x = processed_rna$disease_ancestry_lib_sex, levels = my_levels)
table(unique((processed_rna$disease_ancestry_lib_sex)))

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M_Tulane", "ND_black_HP2110001_M_Tulane", "ND_black_HP2123201_M_Tulane", "ND_black_HPAP-052_M_UPENN", "ND_black_HPAP-080_M_nPOD", #Black M ND
               "ND_black_HP2106201_F_Tulane", "ND_black_HP2121601_F_Tulane", "ND_black_HP2132801_F_Tulane", "ND_black_HP2202101_F_Tulane", #Black F ND
               
               #Hispanic M ND
               "ND_hispanic_HPAP-099_F_UPENN", "ND_hispanic_HPAP-101_F_nPOD", "ND_hispanic_HPAP-105_F_nPOD", #Hispanic F ND
               
               "ND_white_HP2107001_M_Tulane", "ND_white_HP2107901_M_Tulane", "ND_white_HPAP-026_M_nPOD", "ND_white_HPAP-035_M_UPENN", "ND_white_HPAP-040_M_UPENN", "ND_white_HPAP-056_M_UPENN", "ND_white_HPAP-059_M_UPENN", "ND_white_HPAP-075_M_UPENN", "ND_white_HPAP-077_M_UPENN", "ND_white_HPAP-082_M_nPOD", "ND_white_SAMN15877725_M_Tulane", #White M ND
               "ND_white_HP2022801_F_Tulane", "ND_white_HP2024001_F_Tulane", "ND_white_HP2105501_F_Tulane", "ND_white_HP2108601_F_Tulane", "ND_white_HP2108901_F_Tulane", "ND_white_HPAP-022_F_UPENN", "ND_white_HPAP-036_F_nPOD", "ND_white_HPAP-037_F_UPENN", "ND_white_HPAP-053_F_UPENN", "ND_white_HPAP-054_F_UPENN", "ND_white_HPAP-063_F_nPOD",  "ND_white_HPAP-074_F_UPENN", "ND_white_HPAP-103_F_UPENN", #White F ND  
               
               "T2D_black_HPAP-065_M_nPOD", "T2D_black_HPAP-070_M_UPENN", "T2D_black_HPAP-083_M_UPENN", "T2D_black_HPAP-108_M_nPOD", #Black M T2D  
               "T2D_black_HPAP-051_F_UPENN", "T2D_black_HPAP-058_F_nPOD", "T2D_black_HPAP-061_F_nPOD", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F_nPOD", "T2D_hispanic_HPAP-091_F_nPOD", "T2D_hispanic_HPAP-109_F_nPOD", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M_nPOD", "T2D_white_HPAP-100_M_nPOD", "T2D_white_HPAP-106_M_UPENN",# White M T2D
               "T2D_white_HPAP-057_F_UPENN", "T2D_white_HPAP-081_F_nPOD", "T2D_white_HPAP-085_F_UPENN") # White F T2D

table(processed_rna$disease_ancestry_lib_sex_source)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)
table(unique(processed_rna$disease_ancestry_lib_sex_source))

# Make aggregated pseudobulk for conserved marker analysis
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AggregateExpression(processed_rna, return.seurat = TRUE, slot = 'counts')

# Split Metadata and add columns
{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

# Split dataset into M and F
Idents(combined_processed_rna) <- "disease"
combined_processed_rna <- subset(combined_processed_rna, idents = c("ND"))

Idents(combined_processed_rna) <- "sex"
combined_processed_rna_M <- subset(combined_processed_rna, idents = c("M"))
combined_processed_rna_F <- subset(combined_processed_rna, idents = c("F"))

# Analysis for male
Idents(combined_processed_rna_M) <- "celltype"
Idents(combined_processed_rna_F) <- "celltype"

celltypes <- c("beta", "alpha", "beta+alpha", "beta+delta", "gamma", "delta", "epsilon", "cycling-endo", "ductal", "acinar", 
               "quiescent-stellate", "activated-stellate", "endothelial", "macrophages", "mast", "lymphocyte", "schwann")

ident.2 <- c("acinar", "quiescent-stellate", "ductal", "alpha", "endothelial", "gamma", "delta",
             "activated-stellate", "macrophages", "mast", "lymphocyte", "cycling-endo", "epsilon","schwann")

# Now for males
markers_list <- list()
markers_list <- map(celltypes, function(celltype) {
  FindMarkers(combined_processed_rna_M, assay = "RNA", slot = "data", 
              ident.1 = celltype, ident.2 = setdiff(ident.2, celltype), 
              group.by = "celltype", test.use = "DESeq2", 
              only.pos = FALSE) %>% 
    filter(p_val_adj < 1e-2 & avg_log2FC >= 1) %>% 
    arrange(desc(avg_log2FC))
})
names(markers_list) <- celltypes

# Define the directory path where you want to save the files
save_dir <- "C:\\Users\\mqadir\\Box\\Lab 2301\\1. R_Coding Scripts\\Sex Biology Study\\Data Output\\scRNA\\Conserved markers\\DEtesting\\0.by_sex\\male"

# Create the directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir)
}

# Iterate over the list of markers_list
for (i in seq_along(markers_list)) {
  # Get the current cell type
  celltype <- names(markers_list)[i]
  
  # Construct the file path
  file_path <- file.path(save_dir, paste0(celltype, ".csv"))
  
  # Save the data frame to a CSV file
  write.csv(markers_list[[i]], file_path, row.names = TRUE)
  
  print(paste0("Saved file ", celltype, ".csv to ", save_dir))
}

# Now for females
markers_list <- list()
markers_list <- map(celltypes, function(celltype) {
  FindMarkers(combined_processed_rna_F, assay = "RNA", slot = "data", 
              ident.1 = celltype, ident.2 = setdiff(ident.2, celltype), 
              group.by = "celltype", test.use = "DESeq2", 
              only.pos = FALSE) %>% 
    filter(p_val_adj < 1e-2 & avg_log2FC >= 1) %>% 
    arrange(desc(avg_log2FC))
})
names(markers_list) <- celltypes

# Define the directory path where you want to save the files
save_dir <- "C:\\Users\\mqadir\\Box\\Lab 2301\\1. R_Coding Scripts\\Sex Biology Study\\Data Output\\scRNA\\Conserved markers\\DEtesting\\0.by_sex\\female"

# Create the directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir)
}

# Iterate over the list of markers_list
for (i in seq_along(markers_list)) {
  # Get the current cell type
  celltype <- names(markers_list)[i]
  
  # Construct the file path
  file_path <- file.path(save_dir, paste0(celltype, ".csv"))
  
  # Save the data frame to a CSV file
  write.csv(markers_list[[i]], file_path, row.names = TRUE)
  
  print(paste0("Saved file ", celltype, ".csv to ", save_dir))
}

# Make a side by side dotplot
# Advanced coding for ggplot2
# Select only ND samples clear your GE
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(processed_rna) <- "diabetes_status"
processed_rna_ND <- subset(processed_rna, idents = c("ND"))

# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = processed_rna_ND) <- "celltype_sex"
names(table(processed_rna_ND@meta.data$celltype_sex))

# New metadata column is not paired, so we need to pair
my_levels2 <- c("delta_F", "delta_M", "beta+delta_F", "beta+delta_M", "beta_F", "beta_M",  "beta+alpha_F", "beta+alpha_M", "alpha_F", "alpha_M",  "gamma_F", "gamma_M", "epsilon_F", "epsilon_M", "cycling_endo_F", "cycling_endo_M",
                "ductal_F", "ductal_M", "acinar_F", "acinar_M", "activated_stellate_F", "activated_stellate_M", 
                "quiescent_stellate_F", "quiescent_stellate_M", 
                "endothelial_F", "endothelial_M", "lymphocyte_F", "lymphocyte_M", "macrophages_F", "macrophages_M", "mast_F", "mast_M", "schwann_F", "schwann_M") 

head(processed_rna_ND@meta.data$celltype_sex)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna_ND@meta.data$celltype_sex <- factor(x = processed_rna_ND@meta.data$celltype_sex, levels = my_levels2)
table(processed_rna_ND@meta.data$celltype_sex)

# Re select organized idents
Idents(processed_rna_ND) <- "celltype_sex"
Idents(processed_rna_ND) <- "celltype"
DefaultAssay(object = processed_rna_ND) <- "RNA"

markers.to.plot <- c("SST", "LEPR", "LY6H", "BCHE", "INS", "MAFA", "PDX1", "GLP1R", "ONECUT3", "TMPRSS11B", "GCG", "TM4SF4", "ARX", "TTR", 
                     "PPY", "CARTPT", "GHRL", "ANXA13", "UBE2C", "TOP2A", "MKI67", "CCNB2", "CFTR", "MMP7", "KRT23", "TFF1", "CPA1", "PNLIP", "CELA2A", "AMY2A",
                     "SFRP2", "PTGDS", "LUM", "PDGFRA", "RGS5", "FABP4", "CSRP2", "ESM1", "SOX18", "PECAM1", "VWF",
                     "CD3D", "CD2", "TRAC", "CD7", "C1QA", "SDS", "TYROBP", "FCER1G", "TPSB2", "TPSAB1", "HPGDS", "CDH19", "SOX10", "NGFR", "GDNF")

DotPlot(processed_rna_ND,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Older dotplot configuration, shows only percentage not expression
Idents(pancreas.integrated) <- "celltype"
DotPlot(pancreas.integrated, features = rev(markers.to.plot), 
        cols = c("blue", "red"), 
        dot.scale = 8, 
        split.by = "treatment") + 
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme_light() + 
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =10, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))

## ORA
#Gene Ontology plotting
# Load data
# Make a list of all unique genes
# Set the working directory to the location of your files
setwd("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/Conserved markers/DEtesting/0.by_sex/female")

# Get a list of all files in the directory
files <- list.files(pattern = "*.csv")

# Initialize an empty list to store the rownames
gene.list.F <- list()

# Loop through each file and extract the rownames
for (file in files) {
  # Read in the file
  df <- read.csv(file, row.names = 1)
  
  # Extract the rownames
  rownames_vec <- rownames(df)
  
  # Get the file name without the extension
  file_name <- gsub(".csv", "", file)
  
  # Add the rownames to the list with the file name as the name
  gene.list.F[[file_name]] <- rownames_vec
}

# Take a look at the resulting list
gene.list.F

# Set the working directory to the location of your files
setwd("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/Conserved markers/DEtesting/0.by_sex/male")

# Get a list of all files in the directory
files <- list.files(pattern = "*.csv")

# Initialize an empty list to store the rownames
gene.list.M <- list()

# Loop through each file and extract the rownames
for (file in files) {
  # Read in the file
  df <- read.csv(file, row.names = 1)
  
  # Extract the rownames
  rownames_vec <- rownames(df)
  
  # Get the file name without the extension
  file_name <- gsub(".csv", "", file)
  
  # Add the rownames to the list with the file name as the name
  gene.list.M[[file_name]] <- rownames_vec
}

# Take a look at the resulting list
gene.list.M

# Combine lists, but first add notations _M and _F
# Add _F to the names of all objects in the list
names(gene.list.F) <- paste0(names(gene.list.F), "_F")
names(gene.list.M) <- paste0(names(gene.list.M), "_M")
combined_list <- c(gene.list.F, gene.list.M)
names(combined_list)

# Rearrange the list by indexing
combined_list <- combined_list[c("delta_F", "delta_M", "beta+delta_F", "beta+delta_M", "beta_F", "beta_M",  "beta+alpha_F", "beta+alpha_M", "alpha_F", "alpha_M",  "gamma_F", "gamma_M", "epsilon_F", "epsilon_M", "cycling-endo_F", "cycling-endo_M",
                                 "ductal_F", "ductal_M", "acinar_F", "acinar_M", "activated-stellate_F", "activated-stellate_M", 
                                 "quiescent-stellate_F", "quiescent-stellate_M", 
                                 "endothelial_F", "endothelial_M", "lymphocyte_F", "lymphocyte_M", "macrophages_F", "macrophages_M", "mast_F", "mast_M", "schwann_F", "schwann_M")]

# Split dataset into M and F
Idents(processed_rna) <- "diabetes_status"
processed_rna <- subset(processed_rna, idents = c("ND"))

Idents(processed_rna) <- "Sex"
processed_rna_M <- subset(processed_rna, idents = c("M"))
processed_rna_F <- subset(processed_rna, idents = c("F"))

# Compare
ck <- compareCluster(geneCluster = combined_list, 
                     fun = enrichGO, 
                     #universe = rownames(combined_processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db,
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 0.1, #if not set default is at 0.05
                     readable = TRUE)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck) 
cluster_summary <- data.frame(ck)
ck <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 3)
dotplot(ck, showCategory = 1)
ck.save <- ck@compareClusterResult
write.csv(ck.save, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\ORA\ck.all.save.csv)")
write.csv(ck.save, file = r"(C:\Users\mqadir\Documents\ck.all.save.csv)")

dotplot(ck, showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                             "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                             "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                             "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                             "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                             "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                             "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                             "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                             "extracellular matrix organization", "collagen fibril organization",
                             "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                             "regulation of angiogenesis", "blood vessel endothelial cell migration",
                             "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                             "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                             "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                             "central nervous system myelination", "ensheathment of neurons", "axon development"), font.size=5)

cnetplot(ck)

beta.alpha.delta <- list(
  delta_genes=as.character(delta_genes),
  beta_genes=as.character(beta_genes),
  alpha_genes=as.character(alpha_genes)
)


# Compare
ck.bad <- compareCluster(geneCluster = beta.alpha.delta, 
                         fun = enrichGO, 
                         universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = c("ALL"), 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 0.1, #if not set default is at 0.05
                         readable = TRUE)
ck.bad <- setReadable(ck.bad, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
cnetplot(ck.bad,
         showCategory = c("gamma-aminobutyric acid signaling pathway", "hormone secretion",
                          "peptide transport", "peptide hormone secretion", "calcium-ion regulated exocytosis",
                          "neurotransmitter secretion", "Golgi to endosome transport", "potassium channel complex"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 2,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", "dodgerblue3", "lightseagreen"))

options(ggrepel.max.overlaps = Inf)
cnetplot(ck,
         showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                          "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                          "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                          "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                          "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                          "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                          "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                          "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                          "extracellular matrix organization", "collagen fibril organization",
                          "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                          "regulation of angiogenesis", "blood vessel endothelial cell migration",
                          "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                          "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                          "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                          "central nervous system myelination", "ensheathment of neurons", "axon development"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 10,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", #"delta" = 
                                                            "dodgerblue3", #"beta" = ,
                                                            "turquoise2", #"beta+alpha" =
                                                            "lightseagreen", #"alpha"= 
                                                            "springgreen4", #"gamma" =
                                                            "khaki2", #"epsilon" = 
                                                            "darkseagreen2", #"cycling-endo" = 
                                                            "darkorange2", #"ductal" =
                                                            "salmon3", #"acinar" = 
                                                            "orange", #"activated-stellate" = 
                                                            "salmon", #"quiescent-stellate" = 
                                                            "red", #"endothelial" = 
                                                            "orchid1", #"lymphocyte" = 
                                                            "magenta3", #"macrophages" = 
                                                            "red4", #"mast" = 
                                                            "grey30"#"schwann" = 
         ))

eg <- bitr(as.character(alpha_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
edox <- enrichDGN(as.character(eg$ENTREZID), readable = TRUE)
edox <- setReadable(edox, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
edox <- pairwise_termsim(edox)
emapplot(ck)
treeplot(edox)
mutate(edox, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

############################ STAGE  ############################
############################  12b   ############################
# plotting first load seurat object
# First Plot cell Based clustering
# Add metadata
Idents(processed_rna) <- "ancestry"
processed_rna$ancestry_sex <- paste(Idents(processed_rna), processed_rna$'Sex', sep = "_")
table(processed_rna$ancestry_sex)

Idents(processed_rna) <- "celltype_sex_ancestry_disease"
processed_rna$celltype_sex_ancestry_disease_lib <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib)

Idents(processed_rna) <- "celltype_sex_ancestry_disease_lib"
processed_rna$celltype_sex_ancestry_disease_lib_source <- paste(Idents(processed_rna), processed_rna$'tissue_source', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib_source)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', sep = "_")
table(processed_rna$disease_ancestry_lib_sex)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex_source <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'tissue_source', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source)

Idents(processed_rna) <- "diabetes_status"
processed_rna$disease_ancestry_lib_sex_source_celltype <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'tissue_source', processed_rna$'celltype_qadir', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source_celltype)

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M", "ND_black_HP2110001_M", "ND_black_HP2123201_M", "ND_black_HPAP-052_M", "ND_black_HPAP-080_M", #Black M ND
               "ND_black_HP2106201_F", "ND_black_HP2121601_F", "ND_black_HP2132801_F", "ND_black_HP2202101_F", #Black F ND
               
               #Hispanic M ND
               "ND_hispanic_HPAP-099_F", "ND_hispanic_HPAP-101_F", "ND_hispanic_HPAP-105_F", #Hispanic F ND
               
               "ND_white_HP2107001_M", "ND_white_HP2107901_M", "ND_white_HPAP-026_M", "ND_white_HPAP-035_M", "ND_white_HPAP-040_M", "ND_white_HPAP-056_M", "ND_white_HPAP-059_M", "ND_white_HPAP-075_M", "ND_white_HPAP-077_M", "ND_white_HPAP-082_M", "ND_white_SAMN15877725_M", #White M ND
               "ND_white_HP2022801_F", "ND_white_HP2024001_F", "ND_white_HP2105501_F", "ND_white_HP2108601_F", "ND_white_HP2108901_F", "ND_white_HPAP-022_F", "ND_white_HPAP-036_F", "ND_white_HPAP-037_F", "ND_white_HPAP-053_F", "ND_white_HPAP-054_F",  "ND_white_HPAP-063_F", "ND_white_HPAP-074_F", "ND_white_HPAP-103_F", #White F ND  
               
               "T2D_black_HPAP-065_M", "T2D_black_HPAP-070_M", "T2D_black_HPAP-083_M", "T2D_black_HPAP-108_M", #Black M T2D  
               "T2D_black_HPAP-051_F", "T2D_black_HPAP-058_F", "T2D_black_HPAP-061_F", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F", "T2D_hispanic_HPAP-091_F", "T2D_hispanic_HPAP-109_F", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M", "T2D_white_HPAP-100_M", "T2D_white_HPAP-106_M",# White M T2D
               "T2D_white_HPAP-057_F", "T2D_white_HPAP-081_F", "T2D_white_HPAP-085_F") # White F T2D

table(processed_rna$disease_ancestry_lib_sex)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex <- factor(x = processed_rna$disease_ancestry_lib_sex, levels = my_levels)
table(unique((processed_rna$disease_ancestry_lib_sex)))

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M_Tulane", "ND_black_HP2110001_M_Tulane", "ND_black_HP2123201_M_Tulane", "ND_black_HPAP-052_M_UPENN", "ND_black_HPAP-080_M_nPOD", #Black M ND
               "ND_black_HP2106201_F_Tulane", "ND_black_HP2121601_F_Tulane", "ND_black_HP2132801_F_Tulane", "ND_black_HP2202101_F_Tulane", #Black F ND
               
               #Hispanic M ND
               "ND_hispanic_HPAP-099_F_UPENN", "ND_hispanic_HPAP-101_F_nPOD", "ND_hispanic_HPAP-105_F_nPOD", #Hispanic F ND
               
               "ND_white_HP2107001_M_Tulane", "ND_white_HP2107901_M_Tulane", "ND_white_HPAP-026_M_nPOD", "ND_white_HPAP-035_M_UPENN", "ND_white_HPAP-040_M_UPENN", "ND_white_HPAP-056_M_UPENN", "ND_white_HPAP-059_M_UPENN", "ND_white_HPAP-075_M_UPENN", "ND_white_HPAP-077_M_UPENN", "ND_white_HPAP-082_M_nPOD", "ND_white_SAMN15877725_M_Tulane", #White M ND
               "ND_white_HP2022801_F_Tulane", "ND_white_HP2024001_F_Tulane", "ND_white_HP2105501_F_Tulane", "ND_white_HP2108601_F_Tulane", "ND_white_HP2108901_F_Tulane", "ND_white_HPAP-022_F_UPENN", "ND_white_HPAP-036_F_nPOD", "ND_white_HPAP-037_F_UPENN", "ND_white_HPAP-053_F_UPENN", "ND_white_HPAP-054_F_UPENN", "ND_white_HPAP-063_F_nPOD",  "ND_white_HPAP-074_F_UPENN", "ND_white_HPAP-103_F_UPENN", #White F ND  
               
               "T2D_black_HPAP-065_M_nPOD", "T2D_black_HPAP-070_M_UPENN", "T2D_black_HPAP-083_M_UPENN", "T2D_black_HPAP-108_M_nPOD", #Black M T2D  
               "T2D_black_HPAP-051_F_UPENN", "T2D_black_HPAP-058_F_nPOD", "T2D_black_HPAP-061_F_nPOD", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F_nPOD", "T2D_hispanic_HPAP-091_F_nPOD", "T2D_hispanic_HPAP-109_F_nPOD", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M_nPOD", "T2D_white_HPAP-100_M_nPOD", "T2D_white_HPAP-106_M_UPENN",# White M T2D
               "T2D_white_HPAP-057_F_UPENN", "T2D_white_HPAP-081_F_nPOD", "T2D_white_HPAP-085_F_UPENN") # White F T2D

table(processed_rna$disease_ancestry_lib_sex_source)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)
table(unique(processed_rna$disease_ancestry_lib_sex_source))

# Tulane
Idents(processed_rna) <- 'tissue_source'
#tulane_rna <- subset(processed_rna, idents = "Tulane")
#hpap_rna <- subset(processed_rna, idents = c("nPod", "UPenn"))

#Qc plots add complexity info
processed_rna <- Add_Mito_Ribo_Seurat(seurat_object = processed_rna, species = "Human")
processed_rna <- Add_Cell_Complexity_Seurat(seurat_object = processed_rna)

# Visualize QC metrics as a violin plot
# All functions contain
p1 <- QC_Plots_Genes(seurat_object = processed_rna, pt.size = 0,
                     colors_use = c("red4",
                                    "dodgerblue",
                                    "springgreen4"))
p2 <- QC_Plots_UMIs(seurat_object = processed_rna, pt.size = 0,
                    colors_use = c("red4",
                                   "dodgerblue",
                                   "springgreen4"))
p3 <- QC_Plots_Mito(seurat_object = processed_rna, pt.size = 0,
                    colors_use = c("red4",
                                   "dodgerblue",
                                   "springgreen4"))
p4 <- QC_Plots_Complexity(seurat_object = processed_rna, pt.size = 0,
                          colors_use = c("red4",
                                         "dodgerblue",
                                         "springgreen4"))

wrap_plots(p1, p2, p3, p4, ncol = 2)

#QC all scCustomise
raw_metrics <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Raw Data\metrics_meta_all.csv)", header = TRUE, sep = ',')
raw_metrics$ploting_var <- paste(raw_metrics$origin, raw_metrics$disease, sep = "_")

# Plot
Seq_QC_Plot_Basic_Combined(metrics_dataframe = raw_metrics, plot_by = "origin", significance = T)
Seq_QC_Plot_Alignment_Combined(metrics_dataframe = raw_metrics, plot_by = "origin", significance = T) 

+ theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black"))

#Sunburst plot
meta <- read.csv(file = r"(C:\Users\mqadir\Box\!FAHD\4. Sex and Race Based Study Project\meta.csv)")
tree <- d3_nest(meta, value_cols = "size")
sb1 <- sunburst(tree)

pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\plot.pdf)",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

sunburst(
  tree,
  legend = FALSE,
  width = "100%",
  height = 400
)


DimPlot(processed_rna, #switch here to plot
        reduction = "umap",
        #split.by = "Diabetes Status", 
        group.by = "Sex", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.05,
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange2",      #ductal
                 "salmon3",          #acinar
                 "orange",           #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)

table(processed_rna$disease_ancestry_lib_sex_source)
dittoBarPlot(processed_rna, "celltype_qadir", 
             retain.factor.levels = TRUE,
             scale = "count",
             color.panel = c("dodgerblue3",      #beta
                             "turquoise2",       #beta+alpha
                             "lightseagreen",    #alpha
                             "darkseagreen2",    #cycling_endo
                             "khaki2",           #epsilon 
                             "springgreen4",     #gamma
                             "chartreuse3",      #delta
                             "burlywood3",       #beta+delta
                             "darkorange2",      #ductal
                             "salmon3",          #acinar
                             "orange",           #activated_setallate
                             "salmon",           #quiescent_stellate
                             "red",              #endothelial
                             "magenta3",         #macrophages
                             "orchid1",          #lymphocytes
                             "red4",             #mast
                             "grey30"),          #schwann), group.by = "tx")
                             group.by = "disease_ancestry_lib_sex_source") + coord_flip()


# Umap of Diabetes status
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "diabetes_status", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("dodgerblue", #ND
                 "red2"         #T2D  
        ))
# Umap of sex
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", http://127.0.0.1:42565/graphics/plot_zoom_png?width=1160&height=900
        group.by = "Sex", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("red4", #ND
                 "deepskyblue3"         #T2D  
        ))
# Umap of ancestry_sex
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "ancestry", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("black", #ND
                 "darkorange",
                 "deepskyblue3"
                 ))
# Umap of tissue source
DimPlot(processed_rna, 
        #split.by = "Diabetes Status", 
        group.by = "tissue_source", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("dodgerblue",
                 "springgreen4",         
                 "red4"
        )
)

# Plot genes
FeaturePlot(processed_rna, 
            features = c("INS", "GCG", "SST", "GHRL", "PPY", "MKI67",
                         "CFTR", "PNLIP", "SFRP2", "RGS5", "PECAM1", "CD7",
                         "FCER1G", "TPSB2", "CDH19"), 
            cols = c("lightgrey", "red4"),
            order = FALSE,
            ncol = 5,
            pt.size = 2,
            raster = TRUE,
            raster.dpi = c(1024, 1024))

# Plot genes
# Seperate sex
processed_rna$sex_diabetes <- paste(processed_rna$Sex, processed_rna$diabetes_status, sep = "_")
Idents(processed_rna) <- "sex_diabetes"
female_rna <- subset(processed_rna, idents = c("F_ND"))
male_rna <- subset(processed_rna, idents = c("M_ND"))
FeaturePlot(male_rna, 
            features = c("XIST"), 
            cols = c("lightgrey", "red4"),
            #split.by = "sex_diabetes",
            order = TRUE,
            ncol = 1,
            pt.size = 2,
            max.cutoff = 3,
            raster = TRUE,
            raster.dpi = c(1024, 1024))

#Vln plots
sample_colors <- c("firebrick1", "red4", "dodgerblue", "royalblue4")

# Create Plots
Idents(processed_rna) <- "celltype_qadir"
Stacked_VlnPlot(seurat_object = processed_rna, features = c("XIST"), pt.size = 0.1, x_lab_rotate = TRUE,  plot_legend = TRUE,
                colors_use = sample_colors, split.by = "diabetes_status")

c("XIST", "RPS4Y1", 
  "KDM6A", "KDM5D")

# Heatmap
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')

# subset Tulane
#Idents(processed_rna) <- "Tissue Source"
#tulane_rna <- subset(processed_rna, idents = "Tulane")
#Idents(tulane_rna) <- "disease_ancestry_lib_sex_source_celltype"
#combined_processed_rna <- AverageExpression(tulane_rna, return.seurat = TRUE, slot = 'data')

# Split Metadata and add columns
{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

genes.to.plot <- c('INS', 'GCG')

# All sex chromosome coding genes
x.chrom <- read.csv(r"(C:\Users\mqadir\Box\FMJ lab\4. Chromosome Genes\protcoding_genes_accrosschromosomes\RNA_transcripts_chr_X.csv)", row.names = 1)
y.chrom <- read.csv(r"(C:\Users\mqadir\Box\FMJ lab\4. Chromosome Genes\protcoding_genes_accrosschromosomes\RNA_transcripts_chr_Y.csv)", row.names = 1)
x.chrom.genes <- unique(x.chrom$symbol)
y.chrom.genes <- unique(y.chrom$symbol)

# Load all DE genes
M_ND.vs.F_ND_beta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_alpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_beta_alpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_beta_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)

# WARNING THIS IS FOR T2D DATA
M_ND.vs.F_ND_beta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_alpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_beta_alpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_beta_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
M_ND.vs.F_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t',row.names = 1, header = TRUE)
M_ND.vs.F_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)
#M_ND.vs.F_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1)
M_ND.vs.F_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_T2D.vs.F_T2D.tsv)", sep = '\t', row.names = 1, header = TRUE)

mvsf_beta <- rownames(dplyr::filter(M_ND.vs.F_ND_beta, padj < 0.1))
mvsf_alpha <- rownames(dplyr::filter(M_ND.vs.F_ND_alpha, padj < 0.1))
mvsf_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_delta, padj < 0.1))
mvsf_gamma <- rownames(dplyr::filter(M_ND.vs.F_ND_gamma, padj < 0.1))
mvsf_epsilon <- rownames(dplyr::filter(M_ND.vs.F_ND_epsilon, padj < 0.1))
mvsf_acinar <- rownames(dplyr::filter(M_ND.vs.F_ND_acinar, padj < 0.1))
mvsf_ductal <- rownames(dplyr::filter(M_ND.vs.F_ND_ductal, padj < 0.1))
mvsf_activated_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_activated_stellate, padj < 0.1))
mvsf_quiescent_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_quiescent_stellate, padj < 0.1))
mvsf_beta_alpha <- rownames(dplyr::filter(M_ND.vs.F_ND_beta_alpha, padj < 0.1))
mvsf_beta_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_beta_delta, padj < 0.1))
mvsf_cycling_endo <- rownames(dplyr::filter(M_ND.vs.F_ND_cycling_endo, padj < 0.1))
mvsf_lymphocyte <- rownames(dplyr::filter(M_ND.vs.F_ND_lymphocyte, padj < 0.1))
mvsf_macrophages <- rownames(dplyr::filter(M_ND.vs.F_ND_macrophages, padj < 0.1))
mvsf_mast <- rownames(dplyr::filter(M_ND.vs.F_ND_mast, padj < 0.1))
#mvsf_schwann <- rownames(dplyr::filter(M_ND.vs.F_ND_schwann, padj < 0.1))
mvsf_endothelial <- rownames(dplyr::filter(M_ND.vs.F_ND_endothelial, padj < 0.1))

# T2D
mvsf_beta <- rownames(dplyr::filter(M_ND.vs.F_ND_beta, padj < 0.01))
mvsf_alpha <- rownames(dplyr::filter(M_ND.vs.F_ND_alpha, padj < 0.01))
mvsf_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_delta, padj < 0.01))
mvsf_gamma <- rownames(dplyr::filter(M_ND.vs.F_ND_gamma, padj < 0.01))
mvsf_epsilon <- rownames(dplyr::filter(M_ND.vs.F_ND_epsilon, padj < 0.01))
mvsf_acinar <- rownames(dplyr::filter(M_ND.vs.F_ND_acinar, padj < 0.01))
mvsf_ductal <- rownames(dplyr::filter(M_ND.vs.F_ND_ductal, padj < 0.01))
mvsf_activated_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_activated_stellate, padj < 0.01))
mvsf_quiescent_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_quiescent_stellate, padj < 0.01))
mvsf_beta_alpha <- rownames(dplyr::filter(M_ND.vs.F_ND_beta_alpha, padj < 0.01))
mvsf_beta_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_beta_delta, padj < 0.01))
mvsf_cycling_endo <- rownames(dplyr::filter(M_ND.vs.F_ND_cycling_endo, padj < 0.01))
mvsf_lymphocyte <- rownames(dplyr::filter(M_ND.vs.F_ND_lymphocyte, padj < 0.01))
mvsf_macrophages <- rownames(dplyr::filter(M_ND.vs.F_ND_macrophages, padj < 0.01))
mvsf_mast <- rownames(dplyr::filter(M_ND.vs.F_ND_mast, padj < 0.01))
#mvsf_schwann <- rownames(dplyr::filter(M_ND.vs.F_ND_schwann, padj < 0.1))
mvsf_endothelial <- rownames(dplyr::filter(M_ND.vs.F_ND_endothelial, padj < 0.01))

# Following command gives Unique across all
Reduce(intersect, list(mvsf_beta, mvsf_alpha, mvsf_delta, mvsf_gamma, mvsf_epsilon, mvsf_acinar, mvsf_ductal,
                              mvsf_activated_stellate, mvsf_quiescent_stellate, mvsf_beta_alpha, mvsf_beta_delta, 
                              mvsf_cycling_endo, mvsf_lymphocyte, mvsf_macrophages, mvsf_mast, #mvsf_schwann, #T2D has no schwann
                       mvsf_endothelial))

# Make a list of DE genes across all data
unique.de.genes.all <- unique(Reduce(c, list(mvsf_beta, mvsf_alpha, mvsf_delta, mvsf_gamma, mvsf_epsilon, mvsf_acinar, mvsf_ductal,
                                             mvsf_activated_stellate, mvsf_quiescent_stellate, mvsf_beta_alpha, mvsf_beta_delta, 
                                             mvsf_cycling_endo, mvsf_lymphocyte, mvsf_macrophages, mvsf_mast, #mvsf_schwann, #T2D has no schwann
                                             mvsf_endothelial)))

# Make list of all sex chrom genes
unique.y.chrom.genes <- unique(Reduce(c, list(y.chrom.genes)))
unique.x.chrom.genes <- unique(Reduce(c, list(x.chrom.genes)))

# Intersect against DE genes
sex.chrom.genes.intersected.y <- intersect(unique.y.chrom.genes, unique.de.genes.all)
sex.chrom.genes.intersected.x <- intersect(unique.x.chrom.genes, unique.de.genes.all)
intersect(c(unique.x.chrom.genes, unique.y.chrom.genes), unique.de.genes.all)
auto.genes.intersected <- setdiff(unique.de.genes.all, c(sex.chrom.genes.intersected.x, sex.chrom.genes.intersected.y))
plotgenesnow <- intersect(combined_processed_rna@assays[["RNA"]]@counts@Dimnames[[1]], uniquegenes)

# Intersect against ATAC
plotgenesnow <- intersect(unique.de.genes.all, uniquegenes) #Unique genes comes from snATAC data

genes.to.plot <- c("DDX3Y", "EIF1AY",
                   "KDM5D", "NLGN4Y",
                   "RPS4Y1","USP9Y", 
                   "UTY", "ZFY",
                   "XIST", "TSIX",
                   "ZFX", "KDM5C",
                   "SEPTIN6", "EIF1AX",
                   "KDM6A", "PUDP", "DDX3X")

pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 8,
    height = 6)

Idents(combined_processed_rna) <- "disease"
nd.combined <- subset(combined_processed_rna, idents = c("ND"))
t2d.combined <- subset(combined_processed_rna, idents = c("T2D"))

# BE CAREFUL OF THE OBJECT AND GENE LIST INTERACTIONS
dittoHeatmap(
  t2d.combined,
  genes = unique.de.genes.all, #genes.to.plot.intr,
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("ancestry", "sex", "source", "celltype", "disease"),
  #annot.by = c("lib", "sex", "source"),
  order.by = c("sex", "celltype", "ancestry", "disease"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  breaks=seq(-2, 2, length.out=50),
  scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated-stellate" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "beta+alpha" = "turquoise2",
                                        "beta+delta" = "burlywood3",
                                        "cycling-endo" = "darkseagreen2",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "epsilon" = "khaki2",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophages" = "magenta3",
                                        "mast" = "red4",
                                        "quiescent-stellate" = "salmon",
                                        "schwann" = "grey30"),
                           disease = c("ND" = "dodgerblue",
                                       "T2D" = "red2"),
                           sex = c("F" = "red4",
                                   "M" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black",
                                        "hispanic" = "darkorange"),
                           source = c("nPOD" = "dodgerblue",
                                      "Tulane" = "springgreen4",         
                                      "UPENN" = "red4")),
  # # data.out = FALSE,
  # highlight.features = NULL,
  # highlight.genes = NULL,
  # show_colnames = isBulk(object),
  # show_rownames = TRUE,
  # scale = "row",
  # cluster_cols = TRUE,
  cluster_rows = TRUE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  # breaks = NA,
  # complex = FALSE
  #gaps_col = c(460),
  use_raster = TRUE,
  #raster_quality = 5
)
dev.off()
dev.off()


# Violin plot
Idents(processed_rna) <- "ancestry"
processed_rna$ancestry_sex <- paste(Idents(processed_rna), processed_rna$'Sex', sep = '_')
Idents(processed_rna) <- "ancestry_sex"
processed_rna$ancestry_sex_library <- paste(Idents(processed_rna), processed_rna$'Library', sep = '_')
table(processed_rna@meta.data[["ancestry_sex_library"]])

Idents(processed_rna) <- "ancestry_sex_library"
Stacked_VlnPlot(processed_rna, features = c("DDX3Y", "EIF1AY",
                                            "KDM5D", "NLGN4Y",
                                            "RPS4Y1","USP9Y", 
                                            "UTY", "ZFY"), x_lab_rotate = TRUE)

Stacked_VlnPlot(processed_rna, features = c("XIST", "TSIX",
                                            "ZFX", "KDM5C",
                                            "SEPTIN6", "EIF1AX",
                                            "KDM6A", "PUDP", "DDX3X"), x_lab_rotate = TRUE)

# QC
processed_rna[["percent.mt"]] <- PercentageFeatureSet(processed_rna, pattern = "^MT-")
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
VlnPlot(processed_rna, group.by = "Tissue Source", features = feats, pt.size = 0, ncol = 3) +
  NoLegend()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
Idents(processed_rna) <- "Tissue Source"
plot1 <- FeatureScatter(processed_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(processed_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor = TRUE)
plot1 + plot2

# Volcano plot
# Plotting DE genes ###
Idents(processed_rna) <- "celltype_qadir"
beta.cells <- subset(processed_rna, idents = "beta")
plots <- VlnPlot(beta.cells, features = c("INS", "DDIT3", "MIF", "DEPP1", "PLCG2", "IAPP"), group.by = "Sex", 
                 pt.size = 0, combine = TRUE)
plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CO1", "MT-CYB"), group.by = "Sex", 
                 pt.size = 0, combine = TRUE)
wrap_plots(plots = plots, nrow = 1, ncol = 1)

# Load data BE CAREFUL WHAT YOU CALL
volcanodat <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_ND.vs.F_ND_auto.csv)",
                         header = TRUE, sep = ',', row.names = 1) #check seperator
# volcanodat <- volcanodat %>% 
#   mutate_at(c('pvalue'), ~replace_na(.,0.0000000000000001))
#read.table(file.path(x), sep = '\t', row.names = 1) 
#volcanodat <- epsilon.DHT.response

# Intersect against DE genes
degenes_t2d <- rownames(dplyr::filter(volcanodat, padj < 0.01))
intersect(c(unique.x.chrom.genes, unique.y.chrom.genes), degenes_t2d)

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1 note for T2D padj < 0.01 #keyvals[which(volcanodat$log2FoldChange > 1 & volcanodat$padj < 0.01)]
keyvals[which(volcanodat$log2FoldChange > 0 & volcanodat$padj < 0.1)] <- 'royalblue'
names(keyvals)[which(volcanodat$log2FoldChange > 0 & volcanodat$padj < 0.1)] <- 'high'
  
# modify keyvals for variables with fold change < -1
keyvals[which(volcanodat$log2FoldChange < 0 & volcanodat$padj < 0.1)] <- 'red'
names(keyvals)[which(volcanodat$log2FoldChange < 0 & volcanodat$padj < 0.1)] <- 'low'
    
unique(names(keyvals))
    
unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = TRUE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                #selectLab = c('ATF4', 'DNAJB14', 'EDEM1', 'CREBRF', 'HSPA5', 'AMFR', 'DNAJB14', 'UGCG', 'LDLR',
                #              'COX5A', 'COX6A1', 'NDUFS8', 'NDUFS6', 'NDUFA8', 'NDUFB5', 'NDUFB2', 'CCT5', 'CCT7', 'PTGES3', 'WDR83OS'), # use this for labelling genes on plot
                # selectLab = c('DDX3X', 'KDM6A', 'MAPT', 'BARHL1', 'MAPT-IT1', 'SH3BP5', 'AL672277.1', 'JPX', 'XIST', 'TSIX',
                #              'RPS4Y1', 'EIF1AY', 'USP9Y', 'DDX3Y', 'TTTY14', 'KDM5D', 'UTY', 'ZFY', 'AC244213.1', 'NLGN4Y', 'TMSB4Y', 'AC006157.1',
                #              'SRY', 'LINC00278', 'TTTY10', 'AC010889.2', 'AC010889.1', 'AC011297.1', 'HAR1B', 'CD99', 'BCHE', 'ODF3L1', 'EFEMP2',
                #              'DENND1B', 'SOCS6', 'TWSG1', 'ADGRG7'), # use this for labelling genes on plot
                selectLab = c('DDX3X', 'KIAA0408', 'KDM6A', 'RORB', 'ACCS', 'ADAT1', 'ARSD', 'TMEM196',
                              'L672277.1', 'KDM5C', 'ENPEP', 'JPX', 'FXYD2', 'RBMS3', 'NFIB', 'VWA2', 'XIST', 'TSIX',
                              'RPS4Y1', 'EIF1AY', 'USP9Y', 'TTTY14', 'DDX3Y', 'KDM5D', 'ZFY', 'UTY', 'NLGN4Y', 'AC011297.1',
                              'SRY', 'AC244213.1', 'TMSB4Y', 'AC006157.1', 'PRKY', 'FAM162B', 'SDC3', 'C9orf24', 'CXCL12',
                              'TSPAN8', 'TGFBR3L', 'CCDC141', 'FAM3C', 'TTTY10', 'LINC00278', 'AC010889.1',
                              'TBL1Y', 'AC010889.2', 'RAB38', 'PTPRD-AS1', 'SPESP1', 'C5orf58', 'PAIP1', 'NAE1', 'NAT8L'), # use this for labelling genes on plot
                # selectLab = c('RPS4Y1', 'EIF1AY', 'XIST', 'USP9Y', 'DDX3Y', 'UTY', 'ZFY', 'KDM5D',
                #               'TTTY14', 'AC244213.1', 'TMSB4Y', 'NLGN4Y', 'AC006157.1', 'TSIX',
                #               'SRY', 'AC010889.2', 'LINC00278', 'TTTY10', 'CD99', 'BMP2',
                #               'XDH', 'TNFRSF6B', 'WDR25', 'PCDH9', 'LDHA', 'AC010889.1', 'OAT', 'NANOS3',
                #               'AMY2B', 'TUSC3', 'PSENEN'), # use this for labelling genes on plot
                #selectLab = intersect(c(unique.x.chrom.genes, unique.y.chrom.genes), degenes_t2d), # use this for labelling genes on plot
                #selectLab = c('GSTA1', 'SEPTIN6', 'CTAG2', 'SPTSSB', 'LRATD1', 'GGA1', 'DDIT4', 'CXCL2', 'TNFAIP2', 'CXCL3'), # use this for labelling genes on plot
                #encircle = c('VAMP3'),
                boxedLabels = FALSE,
                #xlim = c(-25,25),
                ylim = c(0,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 0.1,
                FCcutoff = c(0, 0), 
                pointSize = c(ifelse((volcanodat$log2FoldChange > 0 & volcanodat$padj < 0.1) | (volcanodat$log2FoldChange < 0 & volcanodat$padj < 0.1), 3, 2)), #changed from T2D 0.01
                #pointSize = 2,
                labSize = 2,
                labFace = 'bold',
                labCol = 'black',
                shape = c(20, 20, 20, 20),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                lengthConnectors = unit(0.01, "npc"),
                arrowheads = TRUE,
                widthConnectors = 1,
                #directionConnectors = "y",
                #hline = c(10e-8),
                min.segment.length = 0,
                typeConnectors = 'open', 
                max.overlaps = 100,
                #min.segment.length = 20,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black',
                raster = TRUE) + theme(axis.text.x = element_text(colour = "black"),
                                                axis.text.y = element_text(colour = "black"),
                                                axis.title.x = element_text(colour = "black"),
                                                axis.title.y = element_text(colour = "black"))

#Gene Ontology plotting
# Load data NON DIABETIC
# Tulane UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Tulane DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\tulane\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Alldat UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Alldat DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Hpap UP
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\UP\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Hpap DOWN
beta.F_white_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.F_white_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_black_ND.vs.F_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_black_ND.vs.F_black_ND.csv)", sep = ',', header = TRUE)
beta.M_ND.vs.F_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_ND.vs.F_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.F_white_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_white_ND.vs.F_white_ND.csv)", sep = ',', header = TRUE)
beta.M_white_ND.vs.M_black_ND <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\hpap\DOWN\beta.M_white_ND.vs.M_black_ND.csv)", sep = ',', header = TRUE)

# Convert frac to dec
beta.F_white_ND.vs.F_black_ND <- beta.F_white_ND.vs.F_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_black_ND.vs.F_black_ND <- beta.M_black_ND.vs.F_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_ND.vs.F_ND <- beta.M_ND.vs.F_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_white_ND.vs.F_white_ND <- beta.M_white_ND.vs.F_white_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

beta.M_white_ND.vs.M_black_ND <- beta.M_white_ND.vs.M_black_ND %>% 
  mutate(GeneRatio = map_dbl(GeneRatio, ~eval(parse(text = .x))))

# Add col identifier
beta.F_white_ND.vs.F_black_ND["condition"] = 'beta.F_white_ND.vs.F_black_ND'
beta.M_black_ND.vs.F_black_ND["condition"] = 'beta.M_black_ND.vs.F_black_ND'
beta.M_ND.vs.F_ND["condition"] = 'beta.M_ND.vs.F_ND'
beta.M_white_ND.vs.F_white_ND["condition"] = 'beta.M_white_ND.vs.F_white_ND'
beta.M_white_ND.vs.M_black_ND["condition"] = 'beta.M_white_ND.vs.M_black_ND'

# Merge
merged.df <- do.call("rbind", list(beta.F_white_ND.vs.F_black_ND, beta.M_black_ND.vs.F_black_ND, beta.M_ND.vs.F_ND, beta.M_white_ND.vs.F_white_ND, beta.M_white_ND.vs.M_black_ND))
table(merged.df$condition)

# data to plot
data_mod <- merged.df %>%                                     
  arrange(desc(qvalue)) %>%
  group_by(condition) %>%
  slice(1:5)
print(data_mod)


ggplot(data = merged.df,
       aes(x = condition, y = Description, 
       color = `qvalue`, size = GeneRatio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=8, face = "bold"), 
        legend.text=element_text(size=8, face = "bold")) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

#Gene Ontology plotting
# Load data
# delta
M_ND.vs.F_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_delta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_delta, padj < 0.1 & log2FoldChange > 0))
mbvsfb_delta <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_delta, padj < 0.1 & log2FoldChange > 0))
mwvsfw_delta <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_delta, padj < 0.1 & log2FoldChange > 0))
fwvsfb_delta <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_delta, padj < 0.1 & log2FoldChange > 0))
fwvsfh_delta <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_delta, padj < 0.1 & log2FoldChange > 0))
mwvsmb_delta <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_delta, padj < 0.1 & log2FoldChange > 0))
fbvsfh_delta <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_delta, padj < 0.1 & log2FoldChange > 0))
# Extract gene lists DOWN
fvsm_delta <- rownames(dplyr::filter(M_ND.vs.F_ND_delta, padj < 0.1 & log2FoldChange < 0))
fbvsmb_delta <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_delta, padj < 0.1 & log2FoldChange < 0))
fwvsmw_delta <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_delta, padj < 0.1 & log2FoldChange < 0))
fbvsfw_delta <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_delta, padj < 0.1 & log2FoldChange < 0))
fhvsfw_delta <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_delta, padj < 0.1 & log2FoldChange < 0))
mbvsmw_delta <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_delta, padj < 0.1 & log2FoldChange < 0))
fhvsfb_delta <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_delta, padj < 0.1 & log2FoldChange < 0))


# Gamma
M_ND.vs.F_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_gamma <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_gamma <- rownames(dplyr::filter(M_ND.vs.F_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_gamma <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_gamma <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_gamma <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_gamma <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_gamma <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_gamma <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_gamma, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_gamma <- rownames(dplyr::filter(M_ND.vs.F_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_gamma <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_gamma <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_gamma <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_gamma <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_gamma <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_gamma <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_gamma, padj < 0.1 & log2FoldChange < -0.000000000014))


# Epsilon
M_ND.vs.F_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
#F_white_ND.vs.F_hispanic_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
#F_black_ND.vs.F_hispanic_ND_epsilon <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_epsilon <- rownames(dplyr::filter(M_ND.vs.F_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_epsilon <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_epsilon <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_epsilon <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
#fwvsfh_epsilon <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_epsilon <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
#fbvsfh_epsilon <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_epsilon, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_epsilon <- rownames(dplyr::filter(M_ND.vs.F_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_epsilon <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_epsilon <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_epsilon <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
#fhvsfw_epsilon <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_epsilon <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))
#fhvsfb_epsilon <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_epsilon, padj < 0.1 & log2FoldChange < -0.000000000014))


# beta+alpha
M_ND.vs.F_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_betaalpha <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_betaalpha <- rownames(dplyr::filter(M_ND.vs.F_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_betaalpha <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_betaalpha <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_betaalpha <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_betaalpha <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_betaalpha <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_betaalpha <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_betaalpha, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_betaalpha <- rownames(dplyr::filter(M_ND.vs.F_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_betaalpha <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_betaalpha <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_betaalpha <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_betaalpha <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_betaalpha <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_betaalpha <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_betaalpha, padj < 0.1 & log2FoldChange < -0.000000000014))


# beta+delta
M_ND.vs.F_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_betadelta <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_betadelta <- rownames(dplyr::filter(M_ND.vs.F_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_betadelta <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_betadelta <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_betadelta <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_betadelta <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_betadelta <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_betadelta <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_betadelta, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_betadelta <- rownames(dplyr::filter(M_ND.vs.F_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_betadelta <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_betadelta <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_betadelta <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_betadelta <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_betadelta <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_betadelta <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_betadelta, padj < 0.1 & log2FoldChange < -0.000000000014))


# cyc_endo
M_ND.vs.F_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_cycling_endo <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_cycling_endo <- rownames(dplyr::filter(M_ND.vs.F_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_cycling_endo <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_cycling_endo <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_cycling_endo <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_cycling_endo <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_cycling_endo <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_cycling_endo <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_cycling_endo, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_cycling_endo <- rownames(dplyr::filter(M_ND.vs.F_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_cycling_endo <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_cycling_endo <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_cycling_endo <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_cycling_endo <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_cycling_endo <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_cycling_endo <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_cycling_endo, padj < 0.1 & log2FoldChange < -0.000000000014))


# acinar
M_ND.vs.F_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_acinar <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_acinar <- rownames(dplyr::filter(M_ND.vs.F_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_acinar <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_acinar <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_acinar <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_acinar<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_acinar <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_acinar <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_acinar, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_acinar <- rownames(dplyr::filter(M_ND.vs.F_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_acinar <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_acinar <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_acinar <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_acinar <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_acinar <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_acinar <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_acinar, padj < 0.1 & log2FoldChange < -0.000000000014))


# ductal
M_ND.vs.F_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_ductal <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_ductal <- rownames(dplyr::filter(M_ND.vs.F_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_ductal <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_ductal <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_ductal <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_ductal<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_ductal <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_ductal <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_ductal, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_ductal <- rownames(dplyr::filter(M_ND.vs.F_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_ductal <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_ductal <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_ductal <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_ductal <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_ductal <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_ductal <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_ductal, padj < 0.1 & log2FoldChange < -0.000000000014))


# activated_stellate
M_ND.vs.F_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_activated_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_activated_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_activated_stellate <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_activated_stellate <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_activated_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_activated_stellate<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_activated_stellate <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_activated_stellate <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_activated_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_activated_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_activated_stellate <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_activated_stellate <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_activated_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_activated_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_activated_stellate <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_activated_stellate <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_activated_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))


# quiescent_stellate
M_ND.vs.F_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_quiescent_stellate <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_quiescent_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_quiescent_stellate <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_quiescent_stellate <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_quiescent_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_quiescent_stellate<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_quiescent_stellate <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_quiescent_stellate <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_quiescent_stellate, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_quiescent_stellate <- rownames(dplyr::filter(M_ND.vs.F_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_quiescent_stellate <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_quiescent_stellate <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_quiescent_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_quiescent_stellate <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_quiescent_stellate <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_quiescent_stellate <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_quiescent_stellate, padj < 0.1 & log2FoldChange < -0.000000000014))


#endothelial
M_ND.vs.F_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_endothelial <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_endothelial <- rownames(dplyr::filter(M_ND.vs.F_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_endothelial <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_endothelial <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_endothelial <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_endothelial<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_endothelial <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_endothelial <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_endothelial, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_endothelial <- rownames(dplyr::filter(M_ND.vs.F_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_endothelial <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_endothelial <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_endothelial <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_endothelial <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_endothelial <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_endothelial <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_endothelial, padj < 0.1 & log2FoldChange < -0.000000000014))


#lymphocyte
M_ND.vs.F_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_lymphocyte <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_lymphocyte <- rownames(dplyr::filter(M_ND.vs.F_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_lymphocyte <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_lymphocyte <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_lymphocyte <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_lymphocyte<- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_lymphocyte <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_lymphocyte <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_lymphocyte, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_lymphocyte <- rownames(dplyr::filter(M_ND.vs.F_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_lymphocyte <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_lymphocyte <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_lymphocyte <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_lymphocyte <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_lymphocyte <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_lymphocyte <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_lymphocyte, padj < 0.1 & log2FoldChange < -0.000000000014))


#macrophages
M_ND.vs.F_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_macrophages <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_macrophages <- rownames(dplyr::filter(M_ND.vs.F_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_macrophages <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_macrophages <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_macrophages <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_macrophages <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_macrophages <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_macrophages <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_macrophages, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_macrophages <- rownames(dplyr::filter(M_ND.vs.F_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_macrophages <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_macrophages <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_macrophages <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_macrophages <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_macrophages <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_macrophages <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_macrophages, padj < 0.1 & log2FoldChange < -0.000000000014))


#mast
M_ND.vs.F_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_hispanic_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
F_black_ND.vs.F_hispanic_ND_mast <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_mast <- rownames(dplyr::filter(M_ND.vs.F_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_mast <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_mast <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_mast <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfh_mast <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_mast <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
fbvsfh_mast <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_mast, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_mast <- rownames(dplyr::filter(M_ND.vs.F_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_mast <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_mast <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_mast <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfw_mast <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_mast <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))
fhvsfb_mast <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_mast, padj < 0.1 & log2FoldChange < -0.000000000014))


#schwann
M_ND.vs.F_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
M_black_ND.vs.F_black_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.F_white_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
F_white_ND.vs.F_black_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
#F_white_ND.vs.F_hispanic_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
M_white_ND.vs.M_black_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
#F_black_ND.vs.F_hispanic_ND_schwann <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
# Extract gene lists UP
mvsf_schwann <- rownames(dplyr::filter(M_ND.vs.F_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
mbvsfb_schwann <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsfw_schwann <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
fwvsfb_schwann <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
#fwvsfh_schwann <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
mwvsmb_schwann <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
#fbvsfh_schwann <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_schwann, padj < 0.1 & log2FoldChange > 0.000000000014))
# Extract gene lists DOWN
fvsm_schwann <- rownames(dplyr::filter(M_ND.vs.F_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsmb_schwann <- rownames(dplyr::filter(M_black_ND.vs.F_black_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
fwvsmw_schwann <- rownames(dplyr::filter(M_white_ND.vs.F_white_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
fbvsfw_schwann <- rownames(dplyr::filter(F_white_ND.vs.F_black_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
#fhvsfw_schwann <- rownames(dplyr::filter(F_white_ND.vs.F_hispanic_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
mbvsmw_schwann <- rownames(dplyr::filter(M_white_ND.vs.M_black_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))
#fhvsfb_schwann <- rownames(dplyr::filter(F_black_ND.vs.F_hispanic_ND_schwann, padj < 0.1 & log2FoldChange < -0.000000000014))


# Make list    
gene.list.sex <- list("fvsm_delta" = fvsm_delta,     "fvsm_gamma" = fvsm_gamma,     "fvsm_epsilon" = fvsm_epsilon,     "fvsm_betaalpha" = fvsm_betaalpha,     "fvsm_betadelta" = fvsm_betadelta,     "fvsm_cycling_endo" = fvsm_cycling_endo,     "fvsm_acinar" = fvsm_acinar,     "fvsm_ductal" = fvsm_ductal,     "fvsm_activated_stellate" = fvsm_activated_stellate,     "fvsm_quiescent_stellate" = fvsm_quiescent_stellate,     "fvsm_endothelial" = fvsm_endothelial,      "fvsm_lymphocyte" = fvsm_lymphocyte,      "fvsm_macrophages" = fvsm_macrophages,      "fvsm_mast" = fvsm_mast,     "fvsm_schwann" = fvsm_schwann,                                                                                                                                                                          
                      "fbvsmb_delta" = fbvsmb_delta, "fbvsmb_gamma" = fbvsmb_gamma, "fbvsmb_epsilon" = fbvsmb_epsilon, "fbvsmb_betaalpha" = fbvsmb_betaalpha, "fbvsmb_betadelta" = fbvsmb_betadelta, "fbvsmb_cycling_endo" = fbvsmb_cycling_endo, "fbvsmb_acinar" = fbvsmb_acinar, "fbvsmb_ductal" = fbvsmb_ductal, "fbvsmb_activated_stellate" = fbvsmb_activated_stellate, "fbvsmb_quiescent_stellate" = fbvsmb_quiescent_stellate, "fbvsmb_endothelial" = fbvsmb_endothelial,  "fbvsmb_lymphocyte" = fbvsmb_lymphocyte,  "fbvsmb_macrophages" = fbvsmb_macrophages,  "fbvsmb_mast" = fbvsmb_mast, "fbvsmb_schwann" = fbvsmb_schwann,
                      "fwvsmw_delta" = fwvsmw_delta, "fwvsmw_gamma" = fwvsmw_gamma, "fwvsmw_epsilon" = fwvsmw_epsilon, "fwvsmw_betaalpha" = fwvsmw_betaalpha, "fwvsmw_betadelta" = fwvsmw_betadelta, "fwvsmw_cycling_endo" = fwvsmw_cycling_endo, "fwvsmw_acinar" = fwvsmw_acinar, "fwvsmw_ductal" = fwvsmw_ductal, "fwvsmw_activated_stellate" = fwvsmw_activated_stellate, "fwvsmw_quiescent_stellate" = fwvsmw_quiescent_stellate, "fwvsmw_endothelial" = fwvsmw_endothelial,  "fwvsmw_lymphocyte" = fwvsmw_lymphocyte,  "fwvsmw_macrophages" = fwvsmw_macrophages,  "fwvsmw_mast" = fwvsmw_mast, "fwvsmw_schwann" = fwvsmw_schwann,
                  
                      "mvsf_delta" = mvsf_delta,     "mvsf_gamma" = mvsf_gamma,     "mvsf_epsilon" = mvsf_epsilon,     "mvsf_betaalpha" = mvsf_betaalpha,     "mvsf_betadelta" = mvsf_betadelta,     "mvsf_cycling_endo" = mvsf_cycling_endo,     "mvsf_acinar" = mvsf_acinar,     "mvsf_ductal" = mvsf_ductal,     "mvsf_activated_stellate" = mvsf_activated_stellate,     "mvsf_quiescent_stellate" = mvsf_quiescent_stellate,     "mvsf_endothelial" = mvsf_endothelial,      "mvsf_lymphocyte" = mvsf_lymphocyte,      "mvsf_macrophages" = mvsf_macrophages,      "mvsf_mast" = mvsf_mast,     "mvsf_schwann" = mvsf_schwann,
                      "mbvsfb_delta" = mbvsfb_delta, "mbvsfb_gamma" = mbvsfb_gamma, "mbvsfb_epsilon" = mbvsfb_epsilon, "mbvsfb_betaalpha" = mbvsfb_betaalpha, "mbvsfb_betadelta" = mbvsfb_betadelta, "mbvsfb_cycling_endo" = mbvsfb_cycling_endo, "mbvsfb_acinar" = mbvsfb_acinar, "mbvsfb_ductal" = mbvsfb_ductal, "mbvsfb_activated_stellate" = mbvsfb_activated_stellate, "mbvsfb_quiescent_stellate" = mbvsfb_quiescent_stellate, "mbvsfb_endothelial" = mbvsfb_endothelial,  "mbvsfb_lymphocyte" = mbvsfb_lymphocyte,  "mbvsfb_macrophages" = mbvsfb_macrophages,  "mbvsfb_mast" = mbvsfb_mast, "mbvsfb_schwann" = mbvsfb_schwann,
                      "mwvsfw_delta" = mwvsfw_delta, "mwvsfw_gamma" = mwvsfw_gamma, "mwvsfw_epsilon" = mwvsfw_epsilon, "mwvsfw_betaalpha" = mwvsfw_betaalpha, "mwvsfw_betadelta" = mwvsfw_betadelta, "mwvsfw_cycling_endo" = mwvsfw_cycling_endo, "mwvsfw_acinar" = mwvsfw_acinar, "mwvsfw_ductal" = mwvsfw_ductal, "mwvsfw_activated_stellate" = mwvsfw_activated_stellate, "mwvsfw_quiescent_stellate" = mwvsfw_quiescent_stellate, "mwvsfw_endothelial" = mwvsfw_endothelial,  "mwvsfw_lymphocyte" = mwvsfw_lymphocyte,  "mwvsfw_macrophages" = mwvsfw_macrophages,  "mwvsfw_mast" = mwvsfw_mast, "mwvsfw_schwann" = mwvsfw_schwann
)

gene.list.anc <- list("fwvsfb_delta" = fwvsfb_delta, "fwvsfb_gamma" = fwvsfb_gamma, "fwvsfb_epsilon" = fwvsfb_epsilon, "fwvsfb_betaalpha" = fwvsfb_betaalpha, "fwvsfb_betadelta" = fwvsfb_betadelta, "fwvsfb_cycling_endo" = fwvsfb_cycling_endo, "fwvsfb_acinar" = fwvsfb_acinar, "fwvsfb_ductal" = fwvsfb_ductal, "fwvsfb_activated_stellate" = fwvsfb_activated_stellate, "fwvsfb_quiescent_stellate" = fwvsfb_quiescent_stellate, "fwvsfb_endothelial" = fwvsfb_endothelial,  "fwvsfb_lymphocyte" = fwvsfb_lymphocyte,  "fwvsfb_macrophages" = fwvsfb_macrophages,  "fwvsfb_mast" = fwvsfb_mast, "fwvsfb_schwann" = fwvsfb_schwann,
                      "fbvsfw_delta" = fbvsfw_delta, "fbvsfw_gamma" = fbvsfw_gamma, "fbvsfw_epsilon" = fbvsfw_epsilon, "fbvsfw_betaalpha" = fbvsfw_betaalpha, "fbvsfw_betadelta" = fbvsfw_betadelta, "fbvsfw_cycling_endo" = fbvsfw_cycling_endo, "fbvsfw_acinar" = fbvsfw_acinar, "fbvsfw_ductal" = fbvsfw_ductal, "fbvsfw_activated_stellate" = fbvsfw_activated_stellate, "fbvsfw_quiescent_stellate" = fbvsfw_quiescent_stellate, "fbvsfw_endothelial" = fbvsfw_endothelial,  "fbvsfw_lymphocyte" = fbvsfw_lymphocyte,  "fbvsfw_macrophages" = fbvsfw_macrophages,  "fbvsfw_mast" = fbvsfw_mast, "fbvsfw_schwann" = fbvsfw_schwann,
                      "fbvsfh_delta" = fbvsfh_delta, "fbvsfh_gamma" = fbvsfh_gamma,                                    "fbvsfh_betaalpha" = fbvsfh_betaalpha, "fbvsfh_betadelta" = fbvsfh_betadelta, "fbvsfh_cycling_endo" = fbvsfh_cycling_endo, "fbvsfh_acinar" = fbvsfh_acinar, "fbvsfh_ductal" = fbvsfh_ductal, "fbvsfh_activated_stellate" = fbvsfh_activated_stellate, "fbvsfh_quiescent_stellate" = fbvsfh_quiescent_stellate, "fbvsfh_endothelial" = fbvsfh_endothelial,  "fbvsfh_lymphocyte" = fbvsfh_lymphocyte,  "fbvsfh_macrophages" = fbvsfh_macrophages,  "fbvsfh_mast" = fbvsfh_mast,
                      "fhvsfb_delta" = fhvsfb_delta, "fhvsfb_gamma" = fhvsfb_gamma,                                    "fhvsfb_betaalpha" = fhvsfb_betaalpha, "fhvsfb_betadelta" = fhvsfb_betadelta, "fhvsfb_cycling_endo" = fhvsfb_cycling_endo, "fhvsfb_acinar" = fhvsfb_acinar, "fhvsfb_ductal" = fhvsfb_ductal, "fhvsfb_activated_stellate" = fhvsfb_activated_stellate, "fhvsfb_quiescent_stellate" = fhvsfb_quiescent_stellate, "fhvsfb_endothelial" = fhvsfb_endothelial,  "fhvsfb_lymphocyte" = fhvsfb_lymphocyte,  "fhvsfb_macrophages" = fhvsfb_macrophages,  "fhvsfb_mast" = fhvsfb_mast,
                      "fhvsfw_delta" = fhvsfw_delta, "fhvsfw_gamma" = fhvsfw_gamma,                                    "fhvsfw_betaalpha" = fhvsfw_betaalpha, "fhvsfw_betadelta" = fhvsfw_betadelta, "fhvsfw_cycling_endo" = fhvsfw_cycling_endo, "fhvsfw_acinar" = fhvsfw_acinar, "fhvsfw_ductal" = fhvsfw_ductal, "fhvsfw_activated_stellate" = fhvsfw_activated_stellate, "fhvsfw_quiescent_stellate" = fhvsfw_quiescent_stellate, "fhvsfw_endothelial" = fhvsfw_endothelial,  "fhvsfw_lymphocyte" = fhvsfw_lymphocyte,  "fhvsfw_macrophages" = fhvsfw_macrophages,  "fhvsfw_mast" = fhvsfw_mast,
                      "fwvsfh_delta" = fwvsfh_delta, "fwvsfh_gamma" = fwvsfh_gamma,                                    "fwvsfh_betaalpha" = fwvsfh_betaalpha, "fwvsfh_betadelta" = fwvsfh_betadelta, "fwvsfh_cycling_endo" = fwvsfh_cycling_endo, "fwvsfh_acinar" = fwvsfh_acinar, "fwvsfh_ductal" = fwvsfh_ductal, "fwvsfh_activated_stellate" = fwvsfh_activated_stellate, "fwvsfh_quiescent_stellate" = fwvsfh_quiescent_stellate, "fwvsfh_endothelial" = fwvsfh_endothelial,  "fwvsfh_lymphocyte" = fwvsfh_lymphocyte,  "fwvsfh_macrophages" = fwvsfh_macrophages,  "fwvsfh_mast" = fwvsfh_mast,
                  
                      "mwvsmb_delta" = mwvsmb_delta, "mwvsmb_gamma" = mwvsmb_gamma, "mwvsmb_epsilon" = mwvsmb_epsilon, "mwvsmb_betaalpha" = mwvsmb_betaalpha, "mwvsmb_betadelta" = mwvsmb_betadelta, "mwvsmb_cycling_endo" = mwvsmb_cycling_endo, "mwvsmb_acinar" = mwvsmb_acinar, "mwvsmb_ductal" = mwvsmb_ductal, "mwvsmb_activated_stellate" = mwvsmb_activated_stellate, "mwvsmb_quiescent_stellate" = mwvsmb_quiescent_stellate, "mwvsmb_endothelial" = mwvsmb_endothelial,  "mwvsmb_lymphocyte" = mwvsmb_lymphocyte,  "mwvsmb_macrophages" = mwvsmb_macrophages,  "mwvsmb_mast" = mwvsmb_mast, "mwvsmb_schwann" = mwvsmb_schwann,
                      "mbvsmw_delta" = mbvsmw_delta, "mbvsmw_gamma" = mbvsmw_gamma, "mbvsmw_epsilon" = mbvsmw_epsilon, "mbvsmw_betaalpha" = mbvsmw_betaalpha, "mbvsmw_betadelta" = mbvsmw_betadelta, "mbvsmw_cycling_endo" = mbvsmw_cycling_endo, "mbvsmw_acinar" = mbvsmw_acinar, "mbvsmw_ductal" = mbvsmw_ductal, "mbvsmw_activated_stellate" = mbvsmw_activated_stellate, "mbvsmw_quiescent_stellate" = mbvsmw_quiescent_stellate, "mbvsmw_endothelial" = mbvsmw_endothelial,  "mbvsmw_lymphocyte" = mbvsmw_lymphocyte,  "mbvsmw_macrophages" = mbvsmw_macrophages,  "mbvsmw_mast" = mbvsmw_mast, "mbvsmw_schwann" = mbvsmw_schwann
                  )


# For diabetes switch cell type to desire
mt2vsm <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_T2D.vs.M_ND.tsv)", sep = '\t', row.names = 1)
ft2vsf <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.F_T2D.vs.F_ND.tsv)", sep = '\t', row.names = 1)

#mt2vsm <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\T2D_containscopies\alpha.deseq.WaldTest.M_T2D.vs.M_ND.tsv)", sep = '\t', row.names = 1)
#ft2vsf <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\T2D_containscopies\alpha.deseq.WaldTest.F_T2D.vs.F_ND.tsv)", sep = '\t', row.names = 1)

# Extract gene lists UP
mt2vsm_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.1 & log2FoldChange > 0))
ft2vsf_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.1 & log2FoldChange > 0))

# Extract gene lists DOWN
mvsmt2_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.1 & log2FoldChange < 0))
fvsft2_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.1 & log2FoldChange < 0))

#all_genes <- rownames(M_ND.vs.F_ND_delta)
# Compare sex
ck.sex <- compareCluster(geneCluster = gene.list.sex, 
                         fun = enrichGO, 
                         #universe = all_genes, 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 0.1, #if not set default is at 0.05
                         readable = TRUE)
ck.sex <- setReadable(ck.sex, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck.sex) 
cluster_summary_sex <- data.frame(ck.sex)
write.csv(cluster_summary_sex, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\clustersummary_all_sex.csv)")

# Compare ancestry::anc
ck.anc <- compareCluster(geneCluster = gene.list.anc, 
                         fun = enrichGO, 
                         #universe = all_genes, 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 0.1, #if not set default is at 0.05
                         readable = TRUE)
ck.anc <- setReadable(ck.anc, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck.anc) 
cluster_summary_anc <- data.frame(ck.anc)
write.csv(cluster_summary_anc, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\clustersummary_all_anc.csv)")

# Compare sex::t2d
gene.list.diab <- list("mt2vsm_cell" = mt2vsm_cell, "ft2vsf_cell" = ft2vsf_cell,
                       "mvsmt2_cell" = mvsmt2_cell, "fvsft2_cell" = fvsft2_cell)
#gene.list.diab <- list("ft2vsf_cell" = ft2vsf_cell, "mt2vsm_cell" = mt2vsm_cell)

ck.dia <- compareCluster(geneCluster = gene.list.diab, 
                         fun = enrichGO, 
                         #universe = rownames(mt2vsm), 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.2, 
                         qvalueCutoff = 0.2, #if not set default is at 0.05
                         readable = TRUE)
ck.dia <- setReadable(ck.dia, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck.dia) 
cluster_summary_anc <- data.frame(ck.dia)
write.csv(cluster_summary_anc, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\ORA\alldata\clustersummary_all_ck_dia_0.1.csv)")

#ck.sub <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
# Beta
dotplot(ck, showCategory = c("positive regulation of endopeptidase activity", "protein acetylation", #f>m
                             "dopamine uptake", "norepinephrine transport", "dosage compensation by inactivation of X chromosome", #bfvsbm
                             "maintenance of cell polarity", "regulation of collateral sprouting", #wf>wm
                             "histone lysine demethylation", "protein dealkylation", #m>f
                             "histone lysine demethylation", "protein dealkylation", #bmvsbf
                             "histone lysine demethylation", "protein dealkylation", "BMP signaling pathway", "male sex determination", "androgen receptor signaling pathway", #wm>wf
                             #bf>wf
                             #bf>hf
                             "neuropeptide signaling pathway", "protein maturation", "calcium ion homeostasis", #hf>bf
                             "I-kappaB kinase/NF-kappaB signaling", "intracellular lipid transport", "ribosomal small subunit biogenesis", #hf>bf
                             "regulation of histone phosphorylation", "tumor necrosis factor production", "cytokine-mediated signaling pathway" #wf>hf
                             ), font.size=14)

#Alpha
dotplot(ck, showCategory = c("histone lysine demethylation", "protein dealkylation", #fvsm
                             "mitochondrial transcription", "dosage compensation by inactivation of X chromosome", #fbvsmb
                             #fwvsmw
                             "histone demethylase activity",  "cytokine activity", #mvsf
                             "sequestering of actin monomers", "protein dealkylation", #mbvsfb
                             "histone lysine demethylation", "protein dealkylation", "co-SMAD binding", #mwvsfw
                             #fwvsfb
                             "error-prone translesion synthesis", "neurotrophin signaling pathway", "regulation of fibroblast growth factor receptor signaling pathway", #fhvsfw
                             "positive regulation of NF-kappaB transcription factor activity", "toll-like receptor signaling pathway" #fwvsfh
), font.size=14)

#all sex
ck.sex <- simplify(ck.sex, cutoff = 0.6, measure = "Wang")
dotplot(ck.sex, showCategory = c("positive regulation of translation", "positive regulation of cellular amide metabolic process", "primary miRNA processing", "toll-like receptor 7 signaling pathway",
                                 "dosage compensation by inactivation of X chromosome", "positive regulation of neurotransmitter secretion", "bicarbonate transport", "gluconeogenesis", "response to electrical stimulus",
                                 "response to peptide hormone", "regulation of translation in response to endoplasmic reticulum stress", "lipid catabolic process", "calcium-mediated signaling using intracellular calcium source",
                                 "positive regulation of T cell migration", "translational initiation", "positive regulation of ERK1 and ERK2 cascade", "extracellular matrix disassembly", "insulin metabolic process",
                                 "regulation of insulin secretion", "iron ion transport", "alcohol metabolic process", "pancreatic juice secretion", "leukotriene metabolic process", "antigen processing and presentation of exogenous peptide antigen via MHC class II",
                                 "mesenchymal stem cell differentiation", "negative regulation of oxidative stress-induced neuron death", "humoral immune response", "positive regulation of protein tyrosine kinase activity",
                                 "positive regulation of ERAD pathway", "endoplasmic reticulum calcium ion homeostasis", "positive regulation of response to endoplasmic reticulum stress", 
                                 "retinoid metabolic process", "carbohydrate catabolic process", "endothelial cell morphogenesis", "vascular associated smooth muscle cell development", "histone demethylation", "protein dealkylation", "cellular response to BMP stimulus", "regulation of interleukin-6-mediated signaling pathway",
                                 "response to transforming growth factor beta", "negative regulation of mast cell activation involved in immune response", "protein O-linked glycosylation via serine", "cellular response to fibroblast growth factor stimulus"
), font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

#all anc
ck.anc <- simplify(ck.anc, cutoff = 0.6, measure = "Wang")
dotplot(ck.anc, showCategory = c("positive regulation of receptor-mediated endocytosis", "iron ion transport", "insulin metabolic process", "regulation of peptide transport", "replicative senescence",
                                 "negative regulation of cell adhesion mediated by integrin", "neuropeptide signaling pathway", "protein secretion", "fibroblast growth factor receptor signaling pathway", "macromolecule glycosylation",
                                 "positive regulation of ERAD pathway", "endoplasmic reticulum calcium ion homeostasis", "positive regulation of response to endoplasmic reticulum stress", "leukotriene metabolic process",
                                 "maintenance of gastrointestinal epithelium", "negative regulation of endopeptidase activity", "transforming growth factor beta production", "cellular response to amino acid stimulus", 
                                 "regulation of peptidase activity", "developmental growth involved in morphogenesis", "response to chemokine", "negative regulation of natural killer cell mediated cytotoxicity",
                                 "negative regulation of cell killing", "mesenchymal stem cell differentiation", "Wnt signaling pathway involved in midbrain dopaminergic neuron differentiation", "positive regulation of DNA-binding transcription factor activity",
                                 "lipid catabolic process", "regulation of NMDA receptor activity", "regulation of AMPA receptor activity", "cellular response to interleukin-1", "cellular response to interferon-gamma",
                                 "positive regulation of triglyceride lipase activity", "retinol metabolic process", "positive regulation of protein kinase A signaling", "negative regulation of voltage-gated calcium channel activity",
                                 "negative regulation of mitochondrion organization", "myelination", "positive regulation of response to oxidative stress", "negative regulation of protein polyubiquitination", "regulation of interleukin-6-mediated signaling pathway",
                                 "positive regulation of tumor necrosis factor-mediated signaling pathway", "regulation of protein K63-linked ubiquitination", "positive regulation of NF-kappaB transcription factor activity", "attachment of spindle microtubules to kinetochore", 
                                 "non-motile cilium assembly", "lipid transport", "positive regulation of protein kinase activity", "positive regulation of cytokine production", "negative regulation of intracellular estrogen receptor signaling pathway",
                                 "calcium ion transmembrane transport", "glycolipid catabolic process", "response to hepatocyte growth factor", "glycosaminoglycan biosynthetic process", "response to histamine", "cellular response to copper ion", "extracellular matrix organization",
                                 "response to hypoxia", "positive regulation of epidermal growth factor-activated receptor activity", "branch elongation of an epithelium", "pancreatic juice secretion", "cellular calcium ion homeostasis", "interleukin-5 production", 
                                 "interleukin-13 production", "response to fibroblast growth factor", "natural killer cell chemotaxis", "T cell migration", "TRAIL-activated apoptotic signaling pathway", "hormone catabolic process", "iron coordination entity transport",
                                 "digestion", "organ or tissue specific immune response"
),
font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

#all diab
#ck.dia <- simplify(ck.dia, cutoff = 0.6, measure = "Wang")
# Beta
dotplot(ck.dia, showCategory = c("positive regulation of receptor recycling", "regulation of morphogenesis of an epithelium", "regulation of steroid hormone biosynthetic process",
                                 "negative regulation of transport", "negative regulation of secretion by cell", "positive regulation of calcium ion import", "innate immune response activating cell surface receptor signaling pathway",
                                 "positive regulation of mRNA metabolic process", "positive regulation of RNA splicing", "ER-nucleus signaling pathway", "mRNA processing", "post-transcriptional regulation of gene expression", "histone H3 acetylation", "histone methylation"), 
        size = "Count",
        #color = "qvalue",
font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("peptide hormone secretion", "insulin secretion", "amino acid import", "glucose homeostasis",
                                 "cellular ion homeostasis", "cholesterol transport", "triglyceride homeostasis", "response to glucose", "dopamine metabolic process", "endocrine system development",
                                 "cellular amino acid metabolic process", "mitochondrial gene expression", "generation of precursor metabolites and energy", "cellular modified amino acid metabolic process",
                                 "mitochondrial translation", "aerobic respiration", "electron transport chain", "carbohydrate metabolic process", "mitochondrial RNA metabolic process"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Alpha
dotplot(ck.dia, showCategory = c("response to chemokine", "cytokine-mediated signaling pathway", "negative regulation of ATP metabolic process", "negative regulation of carbohydrate metabolic process", "cellular glucuronidation",
                                 "peptidyl-lysine trimethylation", "SREBP signaling pathway", "ER-nucleus signaling pathway", "cellular response to sterol depletion"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("linoleic acid metabolic process", "sphingosine biosynthetic process", "protein localization to ciliary membrane", "Golgi to plasma membrane protein transport", "membrane lipid metabolic process", "ceramide metabolic process", "cilium organization",
                                 "mitochondrial respiratory chain complex assembly", "NADH dehydrogenase complex assembly", "aerobic respiration", "protein secretion", "telomerase RNA localization", "pentose-phosphate shunt", "acetyl-CoA metabolic process", "fatty acid beta-oxidation", "protein folding"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Delta
dotplot(ck.dia, showCategory = c("establishment of mitotic spindle orientation", "G protein-coupled receptor signaling pathway", "centrosome localization", "positive regulation of cell cycle process", "Ras protein signal transduction", "organelle localization"), 
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)


dotplot(ck.dia, showCategory = c("chaperone-mediated protein folding", "regulation of protein stability", "NADH dehydrogenase complex assembly", "mitochondrial respiratory chain complex assembly",
                                 "response to retinoic acid", "protein folding", "regulation of establishment of protein localization to telomere", "epithelial cell differentiation", "regulation of telomere maintenance", "electron transport chain", "pentose-phosphate shunt", "mitochondrial electron transport, NADH to ubiquinone", "response to lipid"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Gamma
dotplot(ck.dia, showCategory = c("peptidyl-lysine monomethylation", "protein methylation", "protein alkylation", "peptidyl-lysine modification",
                                 "response to topologically incorrect protein", "ERAD pathway", "PERK-mediated unfolded protein response", "ER-nucleus signaling pathway", "regulation of triglyceride biosynthetic process"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("retinol metabolic process", "cholesterol transport", "cellular hormone metabolic process", "electron transport chain", "ATP synthesis coupled electron transport", "mitochondrial respiratory chain complex I assembly",
                                 "regulation of protein stability", "chaperone-mediated protein folding", "protein stabilization", "cellular amino acid metabolic process", "telomerase RNA localization",
                                 "mitochondrial translation", "mitochondrial gene expression", "positive regulation of telomere maintenance via telomerase", "hormone metabolic process"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Epsilon
dotplot(ck.dia, showCategory = c("multicellular organism growth", "stem cell proliferation", "extracellular matrix disassembly", "negative regulation of stem cell differentiation", "developmental growth"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("negative regulation of amyloid fibril formation", "protein kinase A signaling", "adenylate cyclase-activating G protein-coupled receptor signaling pathway",
                                 "tissue homeostasis", "positive regulation of cytosolic calcium ion concentration", "positive regulation of MAPK cascade", "mitochondrion organization",
                                 "regulation of cell differentiation", "iron ion homeostasis", "cation homeostasis", "metal ion transport", "cellular homeostasis"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Ductal
dotplot(ck.dia, showCategory = c("negative regulation of signaling receptor activity", "negative regulation of execution phase of apoptosis",
                                 "mRNA processing", "sphingolipid catabolic process", "nuclear export", "protein modification by small protein removal", "proteasome-mediated ubiquitin-dependent protein catabolic process", 
                                 "RNA splicing", "ribonucleoprotein complex biogenesis", "response to endoplasmic reticulum stress", "cellular response to topologically incorrect protein", "cellular response to unfolded protein", "ERAD pathway", "positive regulation of endoplasmic reticulum unfolded protein response"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("biomineralization", "organic acid catabolic process", "pentose-phosphate shunt", "mitochondrial respiratory chain complex assembly", "ATP metabolic process",
                                 "protein secretion", "regulation of proton transport", "glucose 6-phosphate metabolic process", "fatty acid metabolic process", "protein maturation",
                                 "electron transport chain", "pyruvate metabolic process", "cell redox homeostasis"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Acinar
dotplot(ck.dia, showCategory = c("mRNA processing", "post-translational protein targeting to endoplasmic reticulum membrane", "ribonucleoprotein complex biogenesis",
                                 "transcription by RNA polymerase III", "post-translational protein targeting to endoplasmic reticulum membrane", "protein dephosphorylation",
                                 "regulation of BMP signaling pathway", "mRNA polyadenylation", "response to endoplasmic reticulum stress", "protein targeting to ER"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("biomineralization", "organic acid catabolic process", "pentose-phosphate shunt", "mitochondrial respiratory chain complex assembly", "ATP metabolic process",
                                 "protein secretion", "regulation of proton transport", "glucose 6-phosphate metabolic process", "fatty acid metabolic process", "protein maturation",
                                 "electron transport chain", "pyruvate metabolic process", "cell redox homeostasis"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Qstell
dotplot(ck.dia, showCategory = c("response to insulin", "regulation of peptide transport", "regulation of homotypic cell-cell adhesion",
                                 "mRNA processing", "RNA splicing", "response to endoplasmic reticulum stress", "regulation of translation", "signal transduction in response to DNA damage", "post-translational protein targeting to endoplasmic reticulum membrane"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("positive regulation of interleukin-6 production", "actin filament organization", "cellular oxidant detoxification", "ATP metabolic process",
                                 "substrate adhesion-dependent cell spreading", "regulation of actin cytoskeleton organization", "collagen metabolic process", "muscle contraction", "Arp2/3 complex-mediated actin nucleation",
                                 "regulation of cell-substrate adhesion", "fibroblast migration", "actomyosin structure organization", "integrin-mediated signaling pathway"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Astell
dotplot(ck.dia, showCategory = c("transcription initiation from RNA polymerase I promoter", "cell fate determination", "positive regulation of epithelial cell differentiation", "negative regulation of smoothened signaling pathway",
                                 "mRNA processing", "macromolecule methylation", "RNA splicing", "RNA polyadenylation", "DNA replication", "hippo signaling", "positive regulation of cell cycle"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("peptidyl-proline hydroxylation to 4-hydroxy-L-proline", "positive regulation of transforming growth factor beta receptor signaling pathway",
                                 "N-glycan processing", "cell redox homeostasis", "regulation of epithelial to mesenchymal transition", "collagen fibril organization",
                                 "ATP metabolic process", "extracellular matrix organization", "aerobic respiration", "electron transport chain", "mitochondrial respiratory chain complex I assembly", "mitotic cytokinesis",
                                 "mitochondrial translation", "protein maturation", "regulation of cell morphogenesis", "protein processing", "androgen biosynthetic process"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Endothelial
dotplot(ck.dia, showCategory = c("mRNA processing", "mRNA destabilization", "negative regulation of translation", "response to starvation", "mitochondrion disassembly",
                                 "autophagy of nucleus", "signal transduction in response to DNA damage", "organelle disassembly", "stem cell population maintenance"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("regulation of actin cytoskeleton organization", "regulation of telomerase RNA localization to Cajal body", "ATP metabolic process", "aerobic respiration", "mitochondrial gene expression",
                                 "cellular oxidant detoxification", "protein folding", "oxidative phosphorylation", "cell redox homeostasis", "mitochondrial respiratory chain complex I assembly", "telomere maintenance via telomerase",
                                 "SMAD protein complex assembly", "contractile actin filament bundle assembly"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Macro
dotplot(ck.dia, showCategory = c("humoral immune response",
                                 "cellular response to interleukin-1", "granulocyte activation", 
                                 "regulation of cellular response to transforming growth factor beta stimulus", "response to amino acid starvation", "positive regulation by host of viral transcription",
                                 "response to oxygen levels"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("ATP metabolic process", "aerobic respiration", "oxidative phosphorylation", "mitochondrial respiratory chain complex I assembly",
                                 "cellular detoxification", "interleukin-12 production", "interleukin-6 production", "ATP biosynthetic process"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Lymphocyes
dotplot(ck.dia, showCategory = c("negative regulation of execution phase of apoptosis", "cytoplasmic translation", "positive regulation of p38MAPK cascade", "positive regulation of JNK cascade",
                                 "positive regulation of stress-activated MAPK cascade"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("tissue homeostasis", "positive regulation of protein kinase A signaling", "positive regulation of protein kinase B signaling",
                                 "B cell receptor signaling pathway", "regulation of mitochondrion organization", "adenylate cyclase-activating G protein-coupled receptor signaling pathway",
                                 "immunoglobulin production", "calcium ion homeostasis", "regulation of phagocytosis", "acute inflammatory response",
                                 "natural killer cell mediated immunity", "cellular detoxification", "response to interleukin-4", "regulation of protein stability", "negative regulation of response to endoplasmic reticulum stress", "protein stabilization",
                                 "glucose homeostasis", "electron transport chain"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

# Mast
dotplot(ck.dia, showCategory = c("cytoplasmic translation", "negative regulation of ubiquitin protein ligase activity", "negative regulation of miRNA-mediated gene silencing",
                                 "lipid digestion", "retinol metabolic process", "N-glycan processing", "protein localization to cytoskeleton", "protein secretion"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)

dotplot(ck.dia, showCategory = c("regulation of mitochondrial membrane potential", "cellular oxidant detoxification", "icosanoid metabolic process", "proton motive force-driven mitochondrial ATP synthesis",
                                 "regulation of mitochondrial depolarization", "Fc-epsilon receptor signaling pathway", "cell redox homeostasis", "regulation of mitochondrion organization"), 
        size = "Count",
        #color = "qvalue",
        font.size=10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size_area(max_size = 5)


# Heatmaps
# For diabetes switch cell type to desire
mt2vsm <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_T2D.vs.M_ND.tsv)", sep = '\t', row.names = 1)
ft2vsf <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.F_T2D.vs.F_ND.tsv)", sep = '\t', row.names = 1)

# Extract gene lists UP
mt2vsm_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.01 & log2FoldChange > 0))
ft2vsf_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.01 & log2FoldChange > 0))

# Extract gene lists DOWN
mvsmt2_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.01 & log2FoldChange < 0))
fvsft2_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.01 & log2FoldChange < 0))

# Heatmap
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')

# subset male and female
Idents(processed_rna) <- "celltype_qadir"
male_rna <- subset(processed_rna, idents = c("beta"))
female_rna <- subset(processed_rna, idents = "beta_F")
Idents(male_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(male_rna, return.seurat = TRUE, slot = 'data')

# Split Metadata and add columns
{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

# Be careful when running this, it was run by individually running nd vs t2d for each celltype
# genes.to.plot.beta.m <- c(mvsmt2_cell, mt2vsm_cell)
# genes.to.plot.beta.f <- c(fvsft2_cell, ft2vsf_cell)
# 
# genes.to.plot.alpha.m <- c(mvsmt2_cell, mt2vsm_cell)
# genes.to.plot.alpha.f <- c(fvsft2_cell, ft2vsf_cell)
# 
# genes.to.plot.delta.m <- c(mvsmt2_cell, mt2vsm_cell)
# genes.to.plot.delta.f <- c(fvsft2_cell, ft2vsf_cell)
# 
# genes.to.plot.gamma.m <- c(mvsmt2_cell, mt2vsm_cell)
# genes.to.plot.gamma.f <- c(fvsft2_cell, ft2vsf_cell)
# 
# genes.to.plot.epsilon.m <- c(mvsmt2_cell, mt2vsm_cell)
# genes.to.plot.epsilon.f <- c(fvsft2_cell, ft2vsf_cell)

genes.to.plot <- c(genes.to.plot.beta.f)

# Male beta
label_genes <- c('RIPOR2', 'ADAMTS5', 'RAMP3', 'VAMP3', 'TBX', 'HGF', 
                 'GJA1', 'DKK3', 'POR', 'LGALS3', 'FRMD4A', 'GEM', 'PLCG2', 'FCN1',
                 'HNF1A','GLUL', 'RFX6', 'KCNA5', 'CLTRN', 'HNF4A', 'FFAR4',
                 'CLTRN', 'SLC9A3R1', 'RGS4')

# Female beta
label_genes <- c('ERN1', 'SLC38A2', 'RBMX', 'TRA2A', 'DYRK1A', 'TRA2B', 'PRPF19', 'DAZAP1',
                 'PAX5', 'SETD1B', 'DOT1L', 'ZNF304', 'PHF1', 'MLLT6', 'NELFA', 'SETD1A', 'RLF', 'PRDM2', 
                 'KMT5A', 'MTHFR', 'WDR82', 'BRD4', 'KDM3A', 'COA3', 'GFM2')

pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 8,
    height = 6)
dittoHeatmap(
  object = combined_processed_rna,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = genes.to.plot, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("ancestry", "sex", "source", "disease", "celltype"),
  #annot.by = c("lib", "sex", "source"),
  order.by = c("disease", "celltype"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  # scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated-stellate" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "beta+alpha" = "turquoise2",
                                        "beta+delta" = "burlywood3",
                                        "cycling-endo" = "darkseagreen2",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "epsilon" = "khaki2",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophages" = "magenta3",
                                        "mast" = "red4",
                                        "quiescent-stellate" = "salmon",
                                        "schwann" = "grey30"),
                           disease = c("ND" = "dodgerblue",
                                       "T2D" = "red2"),
                           sex = c("F" = "red4",
                                   "M" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black",
                                        "hispanic" = "darkorange"),
                           source = c("nPOD" = "dodgerblue",
                                      "Tulane" = "springgreen4",
                                      "UPENN" = "red4")),
  # # data.out = FALSE,
  # highlight.features = c("INS", "MAFA", "IAPP", "ENTPD3", "NKX6-1", "PDX1", #beta
  #                        "GCG", "TTR", "IRX2", "ARX", "TM4SF4", "PCSK1N", #alpha
  #                        "SST", "RBP4", "HHEX", "LY6H", "F5", "BHLHE41", #delta
  #                        "PPY", #gamma
  #                        "GHRL", #epsilon
  #                        "TOP2A", "CCNB2", "HMGB2", "CDKN3", "MKI67", "CENPF", #cycendo
  #                        "SPP1", "TFPI2", "KRT19", "ONECUT1", "TM4SF1", #ductal
  #                        "CTRB1", "CTRB2", "PRSS2", "PRSS1", "PNLIP", "CELA2A", #acinar
  #                        "SFRP2", "VIM", "DCN", "COL1A1", "LUM", "PTGDS", #activated
  #                        "GADD45B", "HMGB1", "PDGFRB", "PRDX1", "PTMA", "RGS5", #quiescent
  #                        "PECAM1", "VWF", "SOX18", "FCN3", "CD59", "ESM1", #endo
  #                        "CCL5", "NKG7", "CD3E", "IL32", "TRAC", "HLA-B", #lymphocyte
  #                        "TPSAB1", "TPSB2", #mast
  #                        "CRYAB", "SOX10", "NGFR", "RUNX2", "BTC", "CDH19", #schwann
  #                        "SDS", "C1QB", "CD68", "APOE", "VMO1", "MS4A7"), #macrophages
  #right_annotations = rowAnnotation(foo = anno_mark(at = c(1), labels = c("HHEX"))),
  # show_colnames = isBulk(object),
  # show_rownames = TRUE,
  # scale = "row",
  cluster_row = FALSE,
  # cluster_cols = FALSE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  breaks=seq(-2, 2, length.out=50),
  complex = TRUE,
  #column_km = 1,
  use_raster = TRUE,
  raster_quality = 5,
  #column_split = combined_processed_rna$celltype,
  #border_color = "black",
  #gaps_col = c(51, 103, 155, 207, 259, 311, 343, 385, 437, 489, 541, 593, 645, 696, 748, 799),
  #gaps_row = c(59, 159, 165, 265, 269, 272, 372, 472, 572, 672, 755, 855, 909, 1008, 1046)
) + rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_rna[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1
)
)
dev.off()
dev.off()



# Plotting venn diagrams
#Beta cells ND
deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\beta.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
beta <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\alpha.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
alpha <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\delta.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
delta <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\gamma.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
gamma <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\epsilon.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
epsilon <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\betaalpha.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
betaalpha <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\betadelta.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
betadelta <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\cycling_endo.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
cycendo <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\acinar.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
acinar <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\ductal.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
ductal <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\activated_stellate.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 5e-2 & avg_log2FC > 0.000000000014) #pos
  #deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
activated_stellate <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\quiescent_stellate.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 5e-2 & avg_log2FC > 0.000000000014) #pos
  #deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
quiescent_stellate <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\endothelial.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
endothelial <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\lymphocyte.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
lymphocyte <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\macrophages.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
macrophages <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\mast.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
mast <- deseq2_dge_genes

deseq2_dge <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\deseq2\schwann.csv)", sep = ',', row.names = 1)
if (exists("deseq2_dge")) {
  deseq2_dge$p_val_adj[deseq2_dge$p_val_adj == 0] <- 2e-302
  deseq2_dge_genes <- dplyr::filter(deseq2_dge, p_val_adj < 1e-5 & avg_log2FC > 0.000000000014) #pos
  deseq2_dge_genes <- deseq2_dge_genes %>% slice_max(avg_log2FC, n = 100)
  deseq2_dge_genes <- rownames(deseq2_dge_genes) } else {deseq2_dge <- character()}
schwann <- deseq2_dge_genes

# Plotting venn diagrams
#Beta cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Alpha cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Delta cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\delta.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Gamma cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\gamma.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Epsilon cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\epsilon.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Beta+alpha cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+alpha.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Beta+delta cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta+delta.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)


#Cycling Endo cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\cycling_endo.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Acinar cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\acinar.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Ductal cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\ductal.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Activated stellate cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\activated_stellate.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Quiescent stellate cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\quiescent_stellate.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Endothelial cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\endothelial.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#Lymphocyte cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\lymphocyte.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#macrophages cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\macrophages.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#mast cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\mast.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

#schwann cells ND
wfvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_white_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
wfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_white_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)
bmvsbf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_black_ND.vs.F_black_ND.tsv)", sep = '\t', row.names = 1)
mvsf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_ND.vs.F_ND.tsv)", sep = '\t', row.names = 1)
wmvswf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_white_ND.vs.F_white_ND.tsv)", sep = '\t', row.names = 1)
wmvsbm <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.M_white_ND.vs.M_black_ND.tsv)", sep = '\t', row.names = 1)
bfvshf <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\schwann.deseq.WaldTest.F_black_ND.vs.F_hispanic_ND.tsv)", sep = '\t', row.names = 1)

# Sometimes enough comparisons are not present, for example we dont have enough hispanic males
# if else conditional formatting allows us to eliminate 
# First compare across Sex
#UP
if (exists("bmvsbf")) {
  bmvsbf_genes <- dplyr::filter(bmvsbf, padj < 0.1 & log2FoldChange > 0) # >1.2x
  bmvsbf_genes <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}

if (exists("mvsf")) {
  mvsf_genes <- dplyr::filter(mvsf, padj < 0.1 & log2FoldChange > 0) # >1.2x
  mvsf_genes <- rownames(mvsf_genes)} else {mvsf_genes <- character()}

if (exists("wmvswf")) {
  wmvswf_genes <- dplyr::filter(wmvswf, padj < 0.1 & log2FoldChange > 0) # >1.2x
  wmvswf_genes <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}

#DOWN
if (exists("bmvsbf")) {
  bmvsbf_genes <- dplyr::filter(bmvsbf, padj < 0.1 & log2FoldChange < 0) # >1.2x
  bmvsbf_genes <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}

if (exists("mvsf")) {
  mvsf_genes <- dplyr::filter(mvsf, padj < 0.1 & log2FoldChange < 0) # >1.2x
  mvsf_genes <- rownames(mvsf_genes)} else {mvsf_genes <- character()}

if (exists("wmvswf")) {
  wmvswf_genes <- dplyr::filter(wmvswf, padj < 0.1 & log2FoldChange < 0) # >1.2x
  wmvswf_genes <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}

x <- list(
  bmvsbf = bmvsbf_genes,
  mvsf = mvsf_genes,
  wmvswf = wmvswf_genes
)

# Note how I am switching x depending on whether I want a venn of T2 data or ND
# The gene.list.diab comes from code above used to calculate GO 
# ggVennDiagram https://gaospecial.github.io/ggVennDiagram/articles/fully-customed.html
# For diabetes switch cell type to desire
mt2vsm <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.M_T2D.vs.M_ND.tsv)", sep = '\t', row.names = 1)
ft2vsf <- read.table(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\alpha.deseq.WaldTest.F_T2D.vs.F_ND.tsv)", sep = '\t', row.names = 1)

# Extract gene lists UP
mt2vsm_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.01 & log2FoldChange > 0))
ft2vsf_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.01 & log2FoldChange > 0))

# Extract gene lists DOWN
mvsmt2_cell <- rownames(dplyr::filter(mt2vsm, padj < 0.01 & log2FoldChange < 0))
fvsft2_cell <- rownames(dplyr::filter(ft2vsf, padj < 0.01 & log2FoldChange < 0))

# Compare sex::t2d
gene.list.diab <- list("mt2vsm_cell" = mt2vsm_cell, "ft2vsf_cell" = ft2vsf_cell,
                       "mvsmt2_cell" = mvsmt2_cell, "fvsft2_cell" = fvsft2_cell)

x <- gene.list.diab
x <- list(
  beta_m = beta_m,
  beta_f = beta_f,
  alpha_m = alpha_m,
  alpha_f = alpha_f
)

beta_m <-  as.character(cgenes_beta_male$gene_name)
beta_f <- as.character(cgenes_beta_female$gene_name)
alpha_m <- as.character(cgenes_alpha_male$gene_name)
alpha_f <-  as.character(cgenes_alpha_female$gene_name)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
               data = venn_regionedge(data)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data)) +
  # 4. region label layer
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data)) +
  scale_fill_gradient(low = "white", high = "lightseagreen") + # change color based on celltype lightseagreen dodgerblue3
  # scale_color_manual(values = c("bmvsbf" = "black",
  #                               "mvsf" ="black", 
  #                               "wmvswf" = 'black'),
  # scale_color_manual(values = c("beta_m" = "black",
  #                               "beta_f" ="black",
  #                               "alpha_m" = 'black',
  #                               "alpha_f" = 'black')) +
  scale_color_manual(values = c("mt2vsm_cell" = "black",
                                "ft2vsf_cell" ="black",
                                "mvsmt2_cell" = 'black',
                                "fvsft2_cell" = 'black'), labels = c('D' = 'D = bdiv_human')) +
  coord_equal() +
  theme_void()


# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist


# Second compare across Ancestry
#UP
if (exists("wfvsbf")) {
  wfvsbf_genes <- dplyr::filter(wfvsbf, padj < 0.1 & log2FoldChange > 0.000000000014) # >1.2x
  wfvsbf_genes <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}

if (exists("wfvshf")) {
  wfvshf_genes <- dplyr::filter(wfvshf, padj < 0.1 & log2FoldChange > 0.000000000014) # >1.2x
  wfvshf_genes <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}

if (exists("wmvsbm")) {
  wmvsbm_genes <- dplyr::filter(wmvsbm, padj < 0.1 & log2FoldChange > 0.000000000014) # >1.2x
  wmvsbm_genes <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}

if (exists("bfvshf")) {
  bfvshf_genes <- dplyr::filter(bfvshf, padj < 0.1 & log2FoldChange > 0.000000000014) # >1.2x
  bfvshf_genes <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}

#DOWN
if (exists("wfvsbf")) {
  wfvsbf_genes <- dplyr::filter(wfvsbf, padj < 0.1 & log2FoldChange < -0.000000000014) # >1.2x
  wfvsbf_genes <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}

if (exists("wfvshf")) {
  wfvshf_genes <- dplyr::filter(wfvshf, padj < 0.1 & log2FoldChange < -0.000000000014) # >1.2x
  wfvshf_genes <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}

if (exists("wmvsbm")) {
  wmvsbm_genes <- dplyr::filter(wmvsbm, padj < 0.1 & log2FoldChange < -0.000000000014) # >1.2x
  wmvsbm_genes <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}

if (exists("bfvshf")) {
  bfvshf_genes <- dplyr::filter(bfvshf, padj < 0.1 & log2FoldChange < -0.000000000014) # >1.2x
  bfvshf_genes <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}

x <- list(
  wfvsbf = wfvsbf_genes,
  wfvshf = wfvshf_genes,
  wmvsbm = wmvsbm_genes,
  bfvshf = bfvshf_genes
)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, "(", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "grey30") + # change color based on celltype
  scale_color_manual(values = c("beta.bmvsbf" = "black",
                                "beta.mvsf" ="black", 
                                "beta.wmvswf" = 'black'),
                     labels = c('D' = 'D = bdiv_human')) +
  theme_void()

# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist

# Save gene sets
{
write.csv(beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\beta_100.csv)", row.names = FALSE)
write.csv(alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\alpha_100.csv)", row.names = FALSE)
write.csv(delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\delta_100.csv)", row.names = FALSE)
write.csv(gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\gamma_100.csv)", row.names = FALSE)
write.csv(epsilon, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\epsilon_100.csv)", row.names = FALSE)
write.csv(betaalpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\betaalpha_100.csv)", row.names = FALSE)
write.csv(betadelta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\betadelta_100.csv)", row.names = FALSE)
write.csv(cycendo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\cycendo_100.csv)", row.names = FALSE)
write.csv(acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\acinar_100.csv)", row.names = FALSE)
write.csv(ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\ductal_100.csv)", row.names = FALSE)
write.csv(activated_stellate, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\activated_stellate_100.csv)", row.names = FALSE)
write.csv(quiescent_stellate, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\quiescent_stellate_100.csv)", row.names = FALSE)
write.csv(endothelial, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\endothelial_100.csv)", row.names = FALSE)
write.csv(lymphocyte, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\lymphocyte_100.csv)", row.names = FALSE)
write.csv(macrophages, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\macrophages_100.csv)", row.names = FALSE)
write.csv(mast, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\mast_100.csv)", row.names = FALSE)
write.csv(schwann, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\schwann_100.csv)", row.names = FALSE)
}

{
beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\beta_100.csv)", sep = ',')
alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\alpha_100.csv)", sep = ',')
delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\delta_100.csv)", sep = ',')
gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\gamma_100.csv)", sep = ',')
epsilon <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\epsilon_100.csv)", sep = ',')
betaalpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\betaalpha_100.csv)", sep = ',')
betadelta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\betadelta_100.csv)", sep = ',')
cycendo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\cycendo_100.csv)", sep = ',')
acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\acinar_100.csv)", sep = ',')
ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\ductal_100.csv)", sep = ',')
activated_stellate <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\activated_stellate_100.csv)", sep = ',')
quiescent_stellate <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\quiescent_stellate_100.csv)", sep = ',')
endothelial <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\endothelial_100.csv)", sep = ',')
lymphocyte <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\lymphocyte_100.csv)", sep = ',')
macrophages <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\macrophages_100.csv)", sep = ',')
mast <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\mast_100.csv)", sep = ',')
schwann <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\1. conserved_genes_deseq2\schwann_100.csv)", sep = ',')
}

x <- list(
  betadelta_genes=as.character(betadelta$x),
  delta_genes=as.character(delta$x),
  beta_genes=as.character(beta$x),
  betaalpha_genes=as.character(betaalpha$x),
  alpha_genes=as.character(alpha$x),
  gamma_genes=as.character(gamma$x),
  epsilon_genes=as.character(epsilon$x),
  cycendo_genes=as.character(cycendo$x),
  ductal_genes=as.character(ductal$x),
  acinar_genes=as.character(acinar$x),
  activated_stellate_genes=as.character(activated_stellate$x),
  quiescent_stellate_genes=as.character(quiescent_stellate$x),
  endothelial_genes=as.character(endothelial$x),
  lymphocyte_genes=as.character(lymphocyte$x),
  macrophages_genes=as.character(macrophages$x),
  mast_genes=as.character(mast$x),
  schwann_genes=as.character(schwann$x)
)

# Concatenate and remove dupliates
Reduce(intersect, x)
allgenes <- unlist(x, use.names = FALSE)
allgenes_unique <- unique(allgenes)

# Find markers to make lists of genes different across cells
Idents(processed_rna) <- "Diabetes Status"
nd.pancreas <- subset(processed_rna, idents = c("ND"))

# Heatmap
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"

#DefaultAssay(processed_rna) <- "RNA"
#DefaultAssay(processed_rna) <- "SCT"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')

{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("delta", "beta+delta", "beta", "beta+alpha", "alpha", "gamma", "epsilon", "cycling-endo",
               "ductal", "acinar",
               "activated-stellate", "quiescent-stellate", "endothelial",
               "lymphocyte", "macrophages", "mast", "schwann") 

table(combined_processed_rna$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_rna$celltype <- factor(x = combined_processed_rna$celltype, levels = my_levels)
table(combined_processed_rna$celltype)
Idents(combined_processed_rna) <- "celltype"

genes.to.plot <- as.character(as.character(allgenes_unique))
write.csv(genes.to.plot, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\genes.csv)")
# genes.to.plot <- c("INS", "MAFA", "IAPP", "ENTPD3", "NKX6-1", "PDX1", #beta
#                    "GCG", "TTR", "IRX2", "ARX", "TM4SF4", "PCSK1N", #alpha
#                    "SST", "RBP4", "HHEX", "LY6H", "F5", "BHLHE41", #delta
#                    "PPY", #gamma
#                    "GHRL", #epsilon
#                    "TOP2A", "CCNB2", "HMGB2", "CDKN3", "MKI67", "CENPF", #cycendo
#                    "SPP1", "TFPI2", "KRT19", "ONECUT1", "TM4SF1", #ductal
#                    "CTRB1", "CTRB2", "PRSS2", "PRSS1", "PNLIP", "CELA2A", #acinar
#                    "SFRP2", "VIM", "DCN", "COL1A1", "LUM", "PTGDS", #activated
#                    "GADD45B", "HMGB1", "PDGFRB", "PRDX1", "PTMA", "RGS5", #quiescent
#                    "PECAM1", "VWF", "SOX18", "FCN3", "CD59", "ESM1", #endo
#                    "CCL5", "NKG7", "CD3E", "IL32", "TRAC", "HLA-B", #lymphocyte
#                    "TPSAB1", "TPSB2", #mast
#                    "CRYAB", "SOX10", "NGFR", "RUNX2", "BTC", "CDH19", #schwann
#                    "SDS", "C1QB", "CD68", "APOE", "VMO1", "MS4A7")

# Get the genes we want to label.
label_genes <- c("INS", "MAFA", "GLP1R", "PDX1", #beta
                 "ONECUT3", "TMPRSS11B", #betaalpha
                 "GCG", "TTR", "ARX", "TM4SF4", #alpha
                 "SST", "BCHE", "LEPR", "LY6H", #delta
                 "PPY", "THSD7A", "CARTPT", #gamma
                 "GHRL", "ANXA13", #epsilon
                 "TOP2A", "UBE2C", "CCNB2", "MKI67", #cycendo
                 "CPA1", "PNLIP", "CELA2A", "AMY2A", #acinar
                 "CFTR", "MMP7", "KRT23", "TFF1", #ductal
                 "SFRP2", "PDGFRA", "LUM", "PTGDS", #activated
                 "RGS5", "FABP4", "CSRP2", "SCGA", #quiescent
                 "PECAM1", "VWF", "SOX18", "ESM1", #endo
                 "TRAC", "CD3D", "CD7", "CD2", #lymphocyte
                 "TPSAB1", "TPSB2", "HPGDS", #mast
                 "GDNF", "SOX10", "NGFR", "CDH19", #schwann
                 "SDS", "TYROBP", "FCER1G", "C1QA") #macrophages

label_genes <- c("HHEX", 
                 "BMP5",
                 "AC132217.2",
                 "ABLIM1",
                 "PPY",
                 "GHRL",
                 "ANLN",
                 "ABCC3",
                 "AKR1C3",
                 "A2M",
                 "ADAMTS4",
                 "ADGRF5",
                 "CKLF",
                 "ABCA1",
                 "TPSAB1",
                 "ANK3")

# Get the indices of said genes in the original object
#locs <- match(label_genes, rownames(combined_processed_rna[genes,]))
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 16,
    height = 6)
dittoHeatmap(
  object = combined_processed_rna,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = genes.to.plot, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("ancestry", "sex", "source", "disease", "celltype"),
  #annot.by = c("lib", "sex", "source"),
  order.by = c("celltype", "sex", "ancestry"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  # scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("acinar" = "salmon3",
                                        "activated-stellate" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "beta+alpha" = "turquoise2",
                                        "beta+delta" = "burlywood3",
                                        "cycling-endo" = "darkseagreen2",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "epsilon" = "khaki2",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophages" = "magenta3",
                                        "mast" = "red4",
                                        "quiescent-stellate" = "salmon",
                                        "schwann" = "grey30"),
                           disease = c("ND" = "dodgerblue",
                                       "T2D" = "red2"),
                           sex = c("F" = "red4",
                                   "M" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black",
                                        "hispanic" = "darkorange"),
                           source = c("nPOD" = "dodgerblue",
                                      "Tulane" = "springgreen4",
                                      "UPENN" = "red4")),
  # # data.out = FALSE,
  # highlight.features = c("INS", "MAFA", "IAPP", "ENTPD3", "NKX6-1", "PDX1", #beta
  #                        "GCG", "TTR", "IRX2", "ARX", "TM4SF4", "PCSK1N", #alpha
  #                        "SST", "RBP4", "HHEX", "LY6H", "F5", "BHLHE41", #delta
  #                        "PPY", #gamma
  #                        "GHRL", #epsilon
  #                        "TOP2A", "CCNB2", "HMGB2", "CDKN3", "MKI67", "CENPF", #cycendo
  #                        "SPP1", "TFPI2", "KRT19", "ONECUT1", "TM4SF1", #ductal
  #                        "CTRB1", "CTRB2", "PRSS2", "PRSS1", "PNLIP", "CELA2A", #acinar
  #                        "SFRP2", "VIM", "DCN", "COL1A1", "LUM", "PTGDS", #activated
  #                        "GADD45B", "HMGB1", "PDGFRB", "PRDX1", "PTMA", "RGS5", #quiescent
  #                        "PECAM1", "VWF", "SOX18", "FCN3", "CD59", "ESM1", #endo
  #                        "CCL5", "NKG7", "CD3E", "IL32", "TRAC", "HLA-B", #lymphocyte
  #                        "TPSAB1", "TPSB2", #mast
  #                        "CRYAB", "SOX10", "NGFR", "RUNX2", "BTC", "CDH19", #schwann
  #                        "SDS", "C1QB", "CD68", "APOE", "VMO1", "MS4A7"), #macrophages
  #right_annotations = rowAnnotation(foo = anno_mark(at = c(1), labels = c("HHEX"))),
  # show_colnames = isBulk(object),
  # show_rownames = TRUE,
  # scale = "row",
  cluster_row = FALSE,
  # cluster_cols = FALSE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  breaks=seq(-2, 2, length.out=50),
  complex = TRUE,
  #column_km = 1,
  use_raster = TRUE,
  raster_quality = 5,
  #column_split = combined_processed_rna$celltype,
  #border_color = "black",
  gaps_col = c(51, 103, 155, 207, 259, 311, 343, 385, 437, 489, 541, 593, 645, 696, 748, 799),
  gaps_row = c(59, 159, 165, 265, 269, 272, 372, 472, 572, 672, 755, 855, 909, 1008, 1046)
) + rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_rna[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1
                                   )
                                   )

dev.off()
dev.off()


allseuratgenes <- as.character(rownames(combined_processed_rna@assays[["RNA"]]))
intersect_genes <- intersect(allseuratgenes, allgenes_unique)
write.csv(allseuratgenes, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\allseuratgenes.csv)")
x <- list(
  allseuratgenes = allseuratgenes,
  allgenes_unique = allgenes_unique
  )

venn <- Venn(x)
overlap(venn) # https://cran.r-project.org/web/packages/RVenn/vignettes/vignette.html
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, "(", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "chartreuse3") + # change color based on celltype
  scale_color_manual(values = c("beta.bmvsbf" = "black",
                                "beta.mvsf" ="black", 
                                "beta.wmvswf" = 'black'),
                     labels = c('D' = 'D = bdiv_human')) +
  theme_void()

# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist

#Gene Ontology plotting
# Load data
# Make a list of all unique genes
beta_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\beta.csv)", sep = ',')
alpha_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\alpha.csv)", sep = ',')
delta_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\delta.csv)", sep = ',')
gamma_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\gamma.csv)", sep = ',')
epsilon_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\epsilon.csv)", sep = ',')
betaalpha_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\betaalpha.csv)", sep = ',')
betadelta_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\betadelta.csv)", sep = ',')
cycendo_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\cycling_endo.csv)", sep = ',')
acinar_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\acinar.csv)", sep = ',')
ductal_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\ductal.csv)", sep = ',')
activated_stellate_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\activated_stellate.csv)", sep = ',')
quiescent_stellate_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\quiescent_stellate.csv)", sep = ',')
endothelial_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\endothelial.csv)", sep = ',')
lymphocyte_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\lymphocyte.csv)", sep = ',')
macrophages_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\macrophages.csv)", sep = ',')
mast_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\mast.csv)", sep = ',')
schwann_genes <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\schwann.csv)", sep = ',')

beta_genes <- (beta_genes$X)
alpha_genes <- (alpha_genes$X)
delta_genes <- (delta_genes$X)
gamma_genes <- (gamma_genes$X)
epsilon_genes <- (epsilon_genes$X)
betaalpha_genes <- (betaalpha_genes$X)
betadelta_genes <- (betadelta_genes$X)
cycendo_genes <- (cycendo_genes$X)
ductal_genes <- (ductal_genes$X)
acinar_genes <- (acinar_genes$X)
activated_stellate_genes <- (activated_stellate_genes$X)
quiescent_stellate_genes <- (quiescent_stellate_genes$X)
endothelial_genes <- (endothelial_genes$X)
lymphocyte_genes <- (lymphocyte_genes$X)
mast_genes <- (mast_genes$X)
schwann_genes <- (schwann_genes$X)
macrophages_genes <- (macrophages_genes$X)

gene.list <- list(
  delta_genes=as.character(delta_genes),
  betadelta_genes=as.character(betadelta_genes),
  beta_genes=as.character(beta_genes),
  betaalpha_genes=as.character(betaalpha_genes),
  alpha_genes=as.character(alpha_genes),
  gamma_genes=as.character(gamma_genes),
  epsilon_genes=as.character(epsilon_genes),
  cycendo_genes=as.character(cycendo_genes),
  ductal_genes=as.character(ductal_genes),
  acinar_genes=as.character(acinar_genes),
  activatedstellate_genes=as.character(activated_stellate_genes),
  quiescentstellate_genes=as.character(quiescent_stellate_genes),
  endothelial_genes=as.character(endothelial_genes),
  lymphocyte_genes=as.character(lymphocyte_genes),
  macro_genes=as.character(macrophages_genes),
  mast_genes=as.character(mast_genes),
  schwann_genes=as.character(schwann_genes)
)

# Compare
ck <- compareCluster(geneCluster = gene.list, 
                     fun = enrichGO, 
                     universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 0.1, #if not set default is at 0.05
                     readable = TRUE)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck) 
cluster_summary <- data.frame(ck)
ck <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 3)
dotplot(ck, showCategory = 1)
ck.save <- ck@compareClusterResult
write.csv(ck.save, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\ORA\ck.save.csv)")

dotplot(ck, showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                             "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                             "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                             "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                             "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                             "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                             "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                             "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                             "extracellular matrix organization", "collagen fibril organization",
                             "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                             "regulation of angiogenesis", "blood vessel endothelial cell migration",
                             "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                             "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                             "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                             "central nervous system myelination", "ensheathment of neurons", "axon development"), font.size=14)

cnetplot(ck)

beta.alpha.delta <- list(
    delta_genes=as.character(delta_genes),
    beta_genes=as.character(beta_genes),
    alpha_genes=as.character(alpha_genes)
    )


# Compare
ck.bad <- compareCluster(geneCluster = beta.alpha.delta, 
                     fun = enrichGO, 
                     universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 0.1, #if not set default is at 0.05
                     readable = TRUE)
ck.bad <- setReadable(ck.bad, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
cnetplot(ck.bad,
         showCategory = c("gamma-aminobutyric acid signaling pathway", "hormone secretion",
                          "peptide transport", "peptide hormone secretion", "calcium-ion regulated exocytosis",
                          "neurotransmitter secretion", "Golgi to endosome transport", "potassium channel complex"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 2,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", "dodgerblue3", "lightseagreen"))

options(ggrepel.max.overlaps = Inf)
cnetplot(ck,
         showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                          "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                          "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                          "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                          "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                          "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                          "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                          "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                          "extracellular matrix organization", "collagen fibril organization",
                          "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                          "regulation of angiogenesis", "blood vessel endothelial cell migration",
                          "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                          "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                          "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                          "central nervous system myelination", "ensheathment of neurons", "axon development"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 10,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", #"delta" = 
                                                            "dodgerblue3", #"beta" = ,
                                                            "turquoise2", #"beta+alpha" =
                                                            "lightseagreen", #"alpha"= 
                                                            "springgreen4", #"gamma" =
                                                            "khaki2", #"epsilon" = 
                                                            "darkseagreen2", #"cycling-endo" = 
                                                            "darkorange2", #"ductal" =
                                                            "salmon3", #"acinar" = 
                                                            "orange", #"activated-stellate" = 
                                                            "salmon", #"quiescent-stellate" = 
                                                            "red", #"endothelial" = 
                                                            "orchid1", #"lymphocyte" = 
                                                            "magenta3", #"macrophages" = 
                                                            "red4", #"mast" = 
                                                            "grey30"#"schwann" = 
                                                            ))

eg <- bitr(as.character(alpha_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
edox <- enrichDGN(as.character(eg$ENTREZID), readable = TRUE)
edox <- setReadable(edox, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
edox <- pairwise_termsim(edox)
emapplot(ck)
treeplot(edox)
mutate(edox, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

# Corr plot
# you must have generated a combined pseudobulk seurat obj first (called: combined_processed_rna)
# Idents(combined_processed_rna) <- "celltype"
# cells_tocorr <- subset(combined_processed_rna, idents = c("beta", "endothelial"))
# av.exp <- cor(GetAssayData(cells_tocorr, slot = combined_processed_rna@assays[["RNA"]]@counts))
# corrplot(av.exp, method = "circle",
#          tl.cex = 0.1)
# cor.exp$x <- rownames(cor.exp)
# cor.df <- tidyr::gather(data = cor.exp, y, correlation, c('0', '1', '2'))
# ggplot(cor.df, aes(x, y, fill = correlation)) +
#   geom_tile()

# PCA
# you must have generated a combined pseudobulk seurat obj first (called: combined_processed_rna)
test_rna <- combined_processed_rna
Idents(test_rna) <- "celltype"
test_rna$celltype_sex <- paste(Idents(test_rna), test_rna$'sex', sep = "_")
Idents(test_rna) <- "celltype_sex"
test_rna$celltype_sex_ancestry <- paste(Idents(test_rna), test_rna$'ancestry', sep = "_")
table(test_rna$celltype_sex_ancestry)

Idents(test_rna) <- "celltype"
testing_rna <- subset(test_rna, idents = c("beta"))

testing_rna <- FindVariableFeatures(testing_rna, selection.method = "vst", nfeatures = 2000)
testing_rna <- RunPCA(testing_rna, features = VariableFeatures(object = combined_processed_rna))
testing_rna <- FindNeighbors(testing_rna, dims = 1:10)
testing_rna <- FindClusters(testing_rna, resolution = 0.5)
testing_rna <- RunUMAP(testing_rna, dims = 1:10)

Idents(testing_rna) <- "sex"
DimPlot(testing_rna, reduction = "umap")

#Create pseudobulk matrix from ALL CELL TYPES
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'counts')

# Split Metadata and add columns, this is because the pseudobulk seurat obj loses metadata
{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

# Choice 1 (choose this or Choice 2)
# Make meta df and select samples
meta <- processed_rna@meta.data[,c('Library', 'Sex', 'tissue_source', 'Chemistry', 'ancestry', 'diabetes_status', 'celltype_qadir')]
rownames(meta) <- NULL
meta <- meta[!duplicated(meta),]
meta <- meta %>% 
  rename("tissue_source" = "tissue_source",
         "diabetes" = "diabetes_status")
meta$diabetes_ancestry_library_sex_source_celltype <- paste0(meta$'diabetes', '_', meta$ancestry, '_', meta$Library, '_', meta$Sex, '_', meta$'tissue_source', '_', meta$celltype_qadir)
meta$ancestry_sex <- paste0(meta$ancestry, '_', meta$Sex)
meta$disease_sex <- paste0(meta$'diabetes', '_',meta$Sex)
samples <- meta$diabetes_ancestry_library_sex_source_celltype

# Create bulk matrix, calc var genes for use later (3K)
DefaultAssay(combined_processed_rna) <- "RNA"
combined_processed_rna <- FindVariableFeatures(combined_processed_rna, selection.method = "vst", nfeatures = 3000)
bulk <- as.data.frame(GetAssayData(object = combined_processed_rna, slot = "counts", assay = "RNA"))

# adding +1 to adjust for 0s and selecting var genes only
bulk <- (bulk + 1)
bulk <- dplyr::filter(bulk, row.names(bulk) %in% combined_processed_rna@assays[["RNA"]]@var.features)
bulk <- as.matrix(round(bulk))
length(rownames(bulk))
length(combined_processed_rna@assays[["RNA"]]@var.features)
"INS" %in% combined_processed_rna@assays[["RNA"]]@var.features

# Run deseq2
dds_bulk <- DESeqDataSetFromMatrix(
  countData = bulk,
  meta,
  design= ~Sex + Chemistry + tissue_source) #Sex_ancestry_diabetes
dds_bulk <- estimateSizeFactors(dds_bulk)
dds_bulk <- estimateDispersions(dds_bulk)    
vsd_bulk <- varianceStabilizingTransformation(dds_bulk)

#Batch corr using limma
mat <- assay(vsd_bulk)
mm <- model.matrix(~Sex, colData(vsd_bulk))
mat <- limma::removeBatchEffect(mat, batch=vsd_bulk$tissue_source, design=mm)
assay(vsd_bulk) <- mat

#Make a PCA
options(repr.plot.height = 7, repr.plot.width = 14)
pcaData <- plotPCA(vsd_bulk, intgroup=c('ancestry_sex', 'celltype_qadir'), returnData=TRUE, ntop=5000)
pcaData <- plotPCA(vsd_bulk, intgroup=c('disease_sex', 'celltype_qadir'), returnData=TRUE, ntop=5000)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
ggplot(pcaData, aes(PC1, PC2, color=celltype_qadir, shape=disease_sex)) +
  geom_point(size=3) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
  coord_fixed() + theme_minimal() + scale_color_manual(values = c("dodgerblue3",
                                                                  "turquoise2",
                                                                  "lightseagreen",
                                                                  "darkseagreen2",
                                                                  "khaki2",
                                                                  "springgreen4",
                                                                  "chartreuse3",
                                                                  "burlywood3",
                                                                  "darkorange2",
                                                                  "salmon3",
                                                                  "orange",
                                                                  "salmon",
                                                                  "red",
                                                                  "magenta3",
                                                                  "orchid1",
                                                                  "red4",
                                                                  "grey30"))


#Choice 2
# Make average seurat object only SAMPLES
processed_rna$ancestry_sex_lib <- paste0(processed_rna$ancestry, '_', processed_rna$Sex, '_', processed_rna$Library)
Idents(processed_rna) <- "ancestry_sex_lib"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')
meta <- processed_rna@meta.data[,c('Library', 'Sex', 'tissue_source', 'Chemistry', 'ancestry', 'diabetes_status', 'ancestry_sex_lib')]
rownames(meta) <- NULL
meta <- meta[!duplicated(meta),]
meta <- meta %>% 
  rename("tissue_source" = "tissue_source",
         "diabetes" = "diabetes_status")
meta$Sex_ancestry_diabetes <- paste0(meta$Sex, '_', meta$ancestry, '_',  meta$'diabetes')
meta$ancestry_sex <- paste0(meta$ancestry, '_', meta$Sex)
meta$disease_sex <- paste0(meta$'diabetes', '_',meta$Sex)
samples <- meta$ancestry_sex_lib
DefaultAssay(combined_processed_rna) <- "RNA"
combined_processed_rna <- FindVariableFeatures(combined_processed_rna, selection.method = "vst", nfeatures = 3000)
bulk <- as.data.frame(GetAssayData(object = combined_processed_rna, slot = "counts", assay = "RNA"))

# Some kung fu
bulk <- dplyr::filter(bulk, row.names(bulk) %in% combined_processed_rna@assays[["RNA"]]@var.features)
bulk <- as.matrix(round(bulk))
length(rownames(bulk))
length(combined_processed_rna@assays[["RNA"]]@var.features)

"INS" %in% combined_processed_rna@assays[["RNA"]]@var.features

# Deseq2
dds_bulk <- DESeqDataSetFromMatrix(
  countData = bulk,
  meta,
  design= ~Sex + Chemistry + tissue_source) #Sex_ancestry_diabetes

dds_bulk <- estimateSizeFactors(dds_bulk)
dds_bulk <- estimateDispersions(dds_bulk)    
vsd_bulk <- varianceStabilizingTransformation(dds_bulk)

#Batch corr using limma
mat <- assay(vsd_bulk)
mm <- model.matrix(~Sex, colData(vsd_bulk))
mat <- limma::removeBatchEffect(mat, batch=vsd_bulk$tissue_source, design=mm)
assay(vsd_bulk) <- mat

#Make a PCA
options(repr.plot.height = 7, repr.plot.width = 14)
pcaData <- plotPCA(vsd_bulk, intgroup=c('ancestry_sex', 'tissue_source'), returnData=TRUE, ntop=5000)
pcaData <- plotPCA(vsd_bulk, intgroup=c('ancestry_sex', 'diabetes'), returnData=TRUE, ntop=5000)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
ggplot(pcaData, aes(PC1, PC2, color=tissue_source, shape=ancestry_sex)) +
  geom_point(size=3) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
  coord_fixed() + theme_minimal() + scale_color_manual(values = c("dodgerblue",
                                                                  "springgreen4",         
                                                                  "red4")) + theme(aspect.ratio=1)

ggplot(pcaData, aes(PC1, PC2, color=diabetes, shape=ancestry_sex)) +
  geom_point(size=3) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
  coord_fixed() + theme_minimal() + scale_color_manual(values = c("dodgerblue",
                                                                  "red2")) + theme(aspect.ratio=1)

mat <- assay(vsd_bulk)
mm <- model.matrix(~Sex, colData(vsd_bulk))
mat <- limma::removeBatchEffect(mat, batch=vsd_bulk$tissue_source, design=mm)
assay(vsd_bulk) <- mat
plotPCA(vsd_bulk, intgroup=c('diabetes', 'Sex'))
pcaData <- plotPCA(vsd_bulk, intgroup=c('Sex', 'diabetes'), returnData=TRUE, ntop=5000)

percentVar <- round(100 * attr(pcaData, 'percentVar'))
ggplot(pcaData, aes(PC1, PC2, color=diabetes, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
  coord_fixed() + theme_minimal() + scale_color_manual(values = c("dodgerblue3",
                                                                  "turquoise2",
                                                                  "lightseagreen",
                                                                  "darkseagreen2",
                                                                  "khaki2",
                                                                  "springgreen4",
                                                                  "chartreuse3",
                                                                  "burlywood3",
                                                                  "darkorange2",
                                                                  "salmon3",
                                                                  "orange",
                                                                  "salmon",
                                                                  "red",
                                                                  "magenta3",
                                                                  "orchid1",
                                                                  "red4",
                                                                  "grey30"))

#Make a heatmap
options(repr.plot.height = 15, repr.plot.width = 20)
sampleDists <- dist(t(assay(vsd_bulk)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd_bulk$Library
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# prep sample annotations
samp_annot <- data.frame(row.names = rownames(d), cyl = factor(d[,2]))

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# Corelation plots
# RNA
install.packages("ggcorrplot")
library(ggcorrplot)
seurobj <- processed_rna
cluster.id <- "beta"
gene.list <- c("INS", "GCG", "PDX1")
clusterCorPlot <- function(seurObj, cluster.id, gene.list, assay='RNA', slot='data') {
  gene.dat <- GetAssayData(
    subset(seurObj, idents = cluster.id), 
    assay=assay, 
    slot=slot)[gene.list,]
  ggcorrplot(as.data.frame(t(gene.dat)))
}

processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(processed_rna) <- "celltype_sex_ancestry_disease"
av.exp <- AverageExpression(processed_rna)$RNA
av.exp <- av.exp[rownames(av.exp) %in% processed_rna@assays[["RNA"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, unique(processed_rna@meta.data[["celltype_sex_ancestry_disease"]]))
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")

# Venn diagrams from scRNAseq cell by cell analysis
# Plotting venn diagrams
#Beta cells ND
bfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bfvwf.csv)", sep = ',', row.names = 1)
hfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.hfvwf.csv)", sep = ',', row.names = 1)
bfvbm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bfvbm.csv)", sep = ',', row.names = 1)
fvm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.fvm.csv)", sep = ',', row.names = 1)
wfvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.wfvwm.csv)", sep = ',', row.names = 1)
bmvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bmvwm.csv)", sep = ',', row.names = 1)
hfvbf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.hfvbf.csv)", sep = ',', row.names = 1)

#Alpha cells ND
bfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bfvwf.csv)", sep = ',', row.names = 1)
hfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.hfvwf.csv)", sep = ',', row.names = 1)
bfvbm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bfvbm.csv)", sep = ',', row.names = 1)
fvm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.fvm.csv)", sep = ',', row.names = 1)
wfvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.wfvwm.csv)", sep = ',', row.names = 1)
bmvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bmvwm.csv)", sep = ',', row.names = 1)
hfvbf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.hfvbf.csv)", sep = ',', row.names = 1)

# Sometimes enough comparisons are not present, for example we dont have enough hispanic males
# if else conditional formatting allows us to eliminate 
# First compare across Sex
#UP
if (exists("bfvbm")) {
  bmvsbf_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  bmvsbf_genes <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}

if (exists("fvm")) {
  mvsf_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  mvsf_genes <- rownames(mvsf_genes)} else {mvsf_genes <- character()}

if (exists("wfvwm")) {
  wmvswf_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  wmvswf_genes <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}

#DOWN
if (exists("bfvbm")) {
  bmvsbf_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  bmvsbf_genes <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}

if (exists("fvm")) {
  mvsf_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  mvsf_genes <- rownames(mvsf_genes)} else {mvsf_genes <- character()}

if (exists("wfvwm")) {
  wmvswf_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wmvswf_genes <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}

x <- list(
  bfvbm = bmvsbf_genes,
  fvm = mvsf_genes,
  wfvwm = wmvswf_genes
)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, "(", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "lightseagreen") + # change color based on celltype
  # scale_color_manual(values = c("bmvsbf" = "black",
  #                               "mvsf" ="black", 
  #                               "wmvswf" = 'black'),
  # scale_color_manual(values = c("beta_m" = "black",
  #                               "beta_f" ="black",
  #                               "alpha_m" = 'black',
  #                               "alpha_f" = 'black')) +
  scale_color_manual(values = c("mt2vsm_cell" = "black",
                                "ft2vsf_cell" ="black",
                                "mvsmt2_cell" = 'black',
                                "fvsft2_cell" = 'black'),
                     labels = c('D' = 'D = bdiv_human')) +
  theme_void()

# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist


###
# Second compare across Ancestry
#UP
if (exists("bfvwf")) {
  wfvsbf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  wfvsbf_genes <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}

if (exists("hfvwf")) {
  wfvshf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  wfvshf_genes <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}

if (exists("bmvwm")) {
  wmvsbm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  wmvsbm_genes <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}

if (exists("hfvbf")) {
  bfvshf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  bfvshf_genes <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}

#DOWN
if (exists("bfvwf")) {
  wfvsbf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wfvsbf_genes <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}

if (exists("hfvwf")) {
  wfvshf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wfvshf_genes <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}

if (exists("bmvwm")) {
  wmvsbm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wmvsbm_genes <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}

if (exists("hfvbf")) {
  bfvshf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  bfvshf_genes <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}

x <- list(
  bfvwf = wfvsbf_genes,
  hfvwf = wfvshf_genes,
  bmvwm = wmvsbm_genes,
  hfvbf = bfvshf_genes
)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, "(", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "lightseagreen") + # change color based on celltype
  scale_color_manual(values = c("beta.bmvsbf" = "black",
                                "beta.mvsf" ="black", 
                                "beta.wmvswf" = 'black'),
                     labels = c('D' = 'D = bdiv_human')) +
  theme_void()

# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist

#Gene Ontology plotting
# Load data
# Make a list of all unique genes
#Beta cells ND
{
bfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bfvwf.csv)", sep = ',', row.names = 1)
hfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.hfvwf.csv)", sep = ',', row.names = 1)
bfvbm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bfvbm.csv)", sep = ',', row.names = 1)
fvm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.fvm.csv)", sep = ',', row.names = 1)
wfvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.wfvwm.csv)", sep = ',', row.names = 1)
bmvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.bmvwm.csv)", sep = ',', row.names = 1)
hfvbf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\beta.hfvbf.csv)", sep = ',', row.names = 1)

# Sometimes enough comparisons are not present, for example we dont have enough hispanic males
# if else conditional formatting allows us to eliminate errors
#UP
if (exists("bfvbm")) {
  bfvbm_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  bfvbm_genes.b <- rownames(bfvbm_genes) } else {bfvbm_genes <- character()}

if (exists("fvm")) {
  fvm_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  fvm_genes.b <- rownames(fvm_genes)} else {fvm_genes <- character()}

if (exists("wfvwm")) {
  wfvwm_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  wfvwm_genes.b <- rownames(wfvwm_genes) } else {wfvwm_genes <- character()}

if (exists("bfvwf")) {
  bfvwf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  bfvwf_genes.b <- rownames(bfvwf_genes) } else {bfvwf_genes <- character()}

if (exists("hfvwf")) {
  hfvwf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  hfvwf_genes.b <- rownames(hfvwf_genes)} else {hfvwf_genes <- character()}

if (exists("bmvwm")) {
  bmvwm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  bmvwm_genes.b <- rownames(bmvwm_genes) } else {bmvwm_genes <- character()}

if (exists("hfvbf")) {
  hfvbf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
  hfvbf_genes.b <- rownames(hfvbf_genes) } else {hfvbf_genes <- character()}

#DOWN

if (exists("bfvbm")) {
  bmvsbf_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  bmvsbf_genes.b <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}

if (exists("fvm")) {
  mvsf_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  mvsf_genes.b <- rownames(mvsf_genes)} else {mvsf_genes <- character()}

if (exists("wfvwm")) {
  wmvswf_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wmvswf_genes.b <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}

if (exists("bfvwf")) {
  wfvsbf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wfvsbf_genes.b <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}

if (exists("hfvwf")) {
  wfvshf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wfvshf_genes.b <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}

if (exists("bmvwm")) {
  wmvsbm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  wmvsbm_genes.b <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}

if (exists("hfvbf")) {
  bfvshf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
  bfvshf_genes.b <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}
}

#Alpha cells ND
{
bfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bfvwf.csv)", sep = ',', row.names = 1)
hfvwf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.hfvwf.csv)", sep = ',', row.names = 1)
bfvbm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bfvbm.csv)", sep = ',', row.names = 1)
fvm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.fvm.csv)", sep = ',', row.names = 1)
wfvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.wfvwm.csv)", sep = ',', row.names = 1)
bmvwm <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.bmvwm.csv)", sep = ',', row.names = 1)
hfvbf <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\scAnalysis_mast\alpha.hfvbf.csv)", sep = ',', row.names = 1)

# Sometimes enough comparisons are not present, for example we dont have enough hispanic males
# if else conditional formatting allows us to eliminate errors
#UP
  if (exists("bfvbm")) {
    bfvbm_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    bfvbm_genes.a <- rownames(bfvbm_genes) } else {bfvbm_genes <- character()}
  
  if (exists("fvm")) {
    fvm_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    fvm_genes.a <- rownames(fvm_genes)} else {fvm_genes <- character()}
  
  if (exists("wfvwm")) {
    wfvwm_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    wfvwm_genes.a <- rownames(wfvwm_genes) } else {wfvwm_genes <- character()}
  
  if (exists("bfvwf")) {
    bfvwf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    bfvwf_genes.a <- rownames(bfvwf_genes) } else {bfvwf_genes <- character()}
  
  if (exists("hfvwf")) {
    hfvwf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    hfvwf_genes.a <- rownames(hfvwf_genes)} else {hfvwf_genes <- character()}
  
  if (exists("bmvwm")) {
    bmvwm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    bmvwm_genes.a <- rownames(bmvwm_genes) } else {bmvwm_genes <- character()}
  
  if (exists("hfvbf")) {
    hfvbf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC > 0.000000000014) # >1.2x
    hfvbf_genes.a <- rownames(hfvbf_genes) } else {hfvbf_genes <- character()}

#DOWN
  if (exists("bfvbm")) {
    bmvsbf_genes <- dplyr::filter(bfvbm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    bmvsbf_genes.a <- rownames(bmvsbf_genes) } else {bmvsbf_genes <- character()}
  
  if (exists("fvm")) {
    mvsf_genes <- dplyr::filter(fvm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    mvsf_genes.a <- rownames(mvsf_genes)} else {mvsf_genes <- character()}
  
  if (exists("wfvwm")) {
    wmvswf_genes <- dplyr::filter(wfvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    wmvswf_genes.a <- rownames(wmvswf_genes) } else {wmvswf_genes <- character()}
  
  if (exists("bfvwf")) {
    wfvsbf_genes <- dplyr::filter(bfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    wfvsbf_genes.a <- rownames(wfvsbf_genes) } else {wfvsbf_genes <- character()}
  
  if (exists("hfvwf")) {
    wfvshf_genes <- dplyr::filter(hfvwf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    wfvshf_genes.a <- rownames(wfvshf_genes)} else {wfvshf_genes <- character()}
  
  if (exists("bmvwm")) {
    wmvsbm_genes <- dplyr::filter(bmvwm, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    wmvsbm_genes.a <- rownames(wmvsbm_genes) } else {wmvsbm_genes <- character()}
  
  if (exists("hfvbf")) {
    bfvshf_genes <- dplyr::filter(hfvbf, p_val_adj < 0.05 & avg_log2FC < -0.000000000014) # >1.2x
    bfvshf_genes.a <- rownames(bfvshf_genes) } else {bfvshf_genes <- character()}
}

gene.list <- list(
  #beta
  fvm_genes.b=as.character(fvm_genes.b),
  bfvbm_genes.b=as.character(bfvbm_genes.b),
  wfvwm_genes.b=as.character(wfvwm_genes.b),
 
  bfvwf_genes.b=as.character(bfvwf_genes.b),
  hfvwf_genes.b=as.character(hfvwf_genes.b),
  hfvbf_genes.b=as.character(hfvbf_genes.b),
  
  mvsf_genes.b=as.character(mvsf_genes.b),
  bmvsbf_genes.b=as.character(bmvsbf_genes.b),
  wmvswf_genes.b=as.character(wmvswf_genes.b),
  
  wfvsbf_genes.b=as.character(wfvsbf_genes.b),
  wfvshf_genes.b=as.character(wfvshf_genes.b),
  bfvshf_genes.b=as.character(bfvshf_genes.b),
  
  wmvsbm_genes.b=as.character(wmvsbm_genes.b),
  bmvwm_genes.b=as.character(bmvwm_genes.b),
  
  #alpha
  fvm_genes.a=as.character(fvm_genes.a),
  bfvbm_genes.a=as.character(bfvbm_genes.a),
  wfvwm_genes.a=as.character(wfvwm_genes.a),
  
  bfvwf_genes.a=as.character(bfvwf_genes.a),
  hfvwf_genes.a=as.character(hfvwf_genes.a),
  hfvbf_genes.a=as.character(hfvbf_genes.a),
  
  mvsf_genes.a=as.character(mvsf_genes.a),
  bmvsbf_genes.a=as.character(bmvsbf_genes.a),
  wmvswf_genes.a=as.character(wmvswf_genes.a),
  
  wfvsbf_genes.a=as.character(wfvsbf_genes.a),
  wfvshf_genes.a=as.character(wfvshf_genes.a),
  bfvshf_genes.a=as.character(bfvshf_genes.a),
  
  wmvsbm_genes.a=as.character(wmvsbm_genes.a),
  bmvwm_genes.a=as.character(bmvwm_genes.a)
)

# Compare
ck <- compareCluster(geneCluster = gene.list, 
                     fun = enrichGO, 
                     universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 0.1, #if not set default is at 0.05
                     readable = TRUE)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
head(ck) 
cluster_summary <- data.frame(ck)
ck <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 3)
dotplot(ck, showCategory = 1)
ck.save <- ck@compareClusterResult
write.csv(ck.save, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\ORA\ck.scrna.save.csv)")

dotplot(ck, showCategory = c("chromatin remodeling", "dosage compensation by inactivation of X chromosome",
                             "negative regulation of protein modification process", "positive regulation of cytokine production", "stress response to copper ion", "response to unfolded protein", "response to topologically incorrect protein", "negative regulation of mitochondrion organization", "cellular response to oxidative stress", "PERK-mediated unfolded protein response", "protein folding chaperone", "integrated stress response signaling",
                             "mitochondrial matrix", "mitochondrial protein-containing complex", "proteasomal protein catabolic process", "ATPase complex", "regulation of protein stability", "protein folding", "aerobic respiration",
                             "mRNA processing", "RNA splicing", "positive regulation of translation", "cellular response to insulin stimulus",
                             "positive regulation of protein maturation", "antioxidant activity", "response to calcium ion", "hormone transport", "regulation of protein processing", "hormone secretion",
                             "regulation of translation", "dosage compensation by inactivation of X chromosome",
                             "hormone transport", "glucose homeostasis",
                             "cellular respiration", "aerobic respiration", "ATP metabolic process", "oxidative phosphorylation", "respiratory electron transport chain", "mRNA processing", "chromatin remodeling",
                             "chromatin remodeling", "mRNA processing", "glucose metabolic process", "pancreas development", "SMAD binding", 
                             "positive regulation of protein transport", "regulation of establishment of protein localization to telomere", "telomerase RNA localization", 
                             "mRNA processing", "histone modification", "chromatin remodeling", 
                             "mRNA processing", "aerobic respiration", "oxidative phosphorylation", "protein folding", "endoplasmic reticulum protein-containing complex", "mitochondrial respiratory chain complex assembly", "electron transport chain"
                             ), font.size=2)

cnetplot(ck)

beta.alpha.delta <- list(
  delta_genes=as.character(delta_genes),
  beta_genes=as.character(beta_genes),
  alpha_genes=as.character(alpha_genes)
)


# Compare
ck.bad <- compareCluster(geneCluster = beta.alpha.delta, 
                         fun = enrichGO, 
                         universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                         keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                         OrgDb = org.Hs.eg.db, 
                         ont = c("ALL"), 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 0.1, #if not set default is at 0.05
                         readable = TRUE)
ck.bad <- setReadable(ck.bad, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
cnetplot(ck.bad,
         showCategory = c("gamma-aminobutyric acid signaling pathway", "hormone secretion",
                          "peptide transport", "peptide hormone secretion", "calcium-ion regulated exocytosis",
                          "neurotransmitter secretion", "Golgi to endosome transport", "potassium channel complex"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 2,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", "dodgerblue3", "lightseagreen"))

options(ggrepel.max.overlaps = Inf)
cnetplot(ck,
         showCategory = c("synapse organization", "gamma-aminobutyric acid signaling pathway",
                          "insulin secretion", "cilium assembly", "peptide hormone secretion", 
                          "cellular response to glucose starvation", "neurotransmitter secretion", "amide transport",
                          "neuropeptide signaling pathway", "protein secretion", "glucagon secretion",
                          "glucocorticoid secretion", "growth hormone secretion", "positive regulation of feeding behavior",
                          "nuclear division", "mitotic cell cycle phase transition", "organelle fission",
                          "epithelial cell proliferation", "digestive tract development", "water homeostasis", "organic anion transport", "SMAD protein signal transduction",
                          "digestion", "morphogenesis of a branching structure", "primary alcohol metabolic process",
                          "extracellular matrix organization", "collagen fibril organization",
                          "muscle contraction", "muscle cell differentiation", "regulation of systemic arterial blood pressure by hormone",
                          "regulation of angiogenesis", "blood vessel endothelial cell migration",
                          "T cell activation", "lymphocyte mediated immunity", "T cell selection",
                          "myeloid leukocyte activation", "antigen processing and presentation", "cell chemotaxis",
                          "immune response-regulating cell surface receptor signaling pathway", "mast cell activation", "activation of immune response",
                          "central nervous system myelination", "ensheathment of neurons", "axon development"),
         foldChange = NULL,
         layout = "kk",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "category",
         cex_category = 10,
         cex_gene = 0.5,
         node_label_size = NULL,
         cex_label_category = 1,
         cex_label_gene = 1) + scale_fill_manual(values = c("chartreuse3", #"delta" = 
                                                            "dodgerblue3", #"beta" = ,
                                                            "turquoise2", #"beta+alpha" =
                                                            "lightseagreen", #"alpha"= 
                                                            "springgreen4", #"gamma" =
                                                            "khaki2", #"epsilon" = 
                                                            "darkseagreen2", #"cycling-endo" = 
                                                            "darkorange2", #"ductal" =
                                                            "salmon3", #"acinar" = 
                                                            "orange", #"activated-stellate" = 
                                                            "salmon", #"quiescent-stellate" = 
                                                            "red", #"endothelial" = 
                                                            "orchid1", #"lymphocyte" = 
                                                            "magenta3", #"macrophages" = 
                                                            "red4", #"mast" = 
                                                            "grey30"#"schwann" = 
         ))

eg <- bitr(as.character(alpha_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
edox <- enrichDGN(as.character(eg$ENTREZID), readable = TRUE)
edox <- setReadable(edox, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
edox <- pairwise_termsim(edox)
emapplot(ck)
treeplot(edox)
mutate(edox, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

############################ END ############################
############################ END ############################

# MOGONET

# Statistical testing
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

# Add metadata
Idents(processed_rna) <- "ancestry"
processed_rna$ancestry_sex <- paste(Idents(processed_rna), processed_rna$'Sex', sep = "_")
table(processed_rna$ancestry_sex)

Idents(processed_rna) <- "celltype_sex_ancestry_disease"
processed_rna$celltype_sex_ancestry_disease_lib <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib)

Idents(processed_rna) <- "celltype_sex_ancestry_disease_lib"
processed_rna$celltype_sex_ancestry_disease_lib_source <- paste(Idents(processed_rna), processed_rna$'Tissue Source', sep = "_")
table(processed_rna$celltype_sex_ancestry_disease_lib_source)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', sep = "_")
table(processed_rna$disease_ancestry_lib_sex)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex_source <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'Tissue Source', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source)

Idents(processed_rna) <- "Diabetes Status"
processed_rna$disease_ancestry_lib_sex_source_celltype <- paste(Idents(processed_rna), processed_rna$'ancestry', processed_rna$'Library', processed_rna$'Sex', processed_rna$'Tissue Source', processed_rna$'celltype_qadir', sep = "_")
table(processed_rna$disease_ancestry_lib_sex_source_celltype)

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M", "ND_black_HP2110001_M", "ND_black_HP2123201_M", "ND_black_HPAP-052_M", #Black M ND
               "ND_black_HP2106201_F", "ND_black_HP2121601_F", "ND_black_HP2132801_F", "ND_black_HP2202101_F", #Black F ND
               
               "ND_hispanic_HPAP-080_M", #Hispanic M ND
               "ND_hispanic_HPAP-099_F", "ND_hispanic_HPAP-101_F", "ND_hispanic_HPAP-105_F", #Hispanic F ND
               
               "ND_white_HP2107001_M", "ND_white_HP2107901_M", "ND_white_HPAP-026_M", "ND_white_HPAP-035_M", "ND_white_HPAP-040_M", "ND_white_HPAP-056_M", "ND_white_HPAP-059_M", "ND_white_HPAP-075_M", "ND_white_HPAP-077_M", "ND_white_HPAP-082_M", "ND_white_SAMN15877725_M", #White M ND
               "ND_white_HP2022801_F", "ND_white_HP2024001_F", "ND_white_HP2105501_F", "ND_white_HP2108601_F", "ND_white_HP2108901_F", "ND_white_HPAP-022_F", "ND_white_HPAP-036_F", "ND_white_HPAP-037_F", "ND_white_HPAP-053_F", "ND_white_HPAP-054_F", "ND_white_HPAP-063_F",  "ND_white_HPAP-074_F", "ND_white_HPAP-103_F", #White F ND  
               
               "T2D_black_HPAP-070_M", "T2D_black_HPAP-083_M", "T2D_black_HPAP-108_M", #Black M T2D  
               "T2D_black_HPAP-051_F", "T2D_black_HPAP-058_F", "T2D_black_HPAP-061_F", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F", "T2D_hispanic_HPAP-091_F", "T2D_hispanic_HPAP-109_F", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M", "T2D_white_HPAP-100_M", "T2D_white_HPAP-106_M", "T2D_white_HPAP-065_M_nPod",# White M T2D
               "T2D_white_HPAP-057_F", "T2D_white_HPAP-081_F", "T2D_white_HPAP-085_F") # White F T2D

table(processed_rna$disease_ancestry_lib_sex)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex <- factor(x = processed_rna$disease_ancestry_lib_sex, levels = my_levels)
table(unique((processed_rna$disease_ancestry_lib_sex)))

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("ND_black_HP2031401_M_Tulane", "ND_black_HP2110001_M_Tulane", "ND_black_HP2123201_M_Tulane", "ND_black_HPAP-052_M_UPenn", #Black M ND
               "ND_black_HP2106201_F_Tulane", "ND_black_HP2121601_F_Tulane", "ND_black_HP2132801_F_Tulane", "ND_black_HP2202101_F_Tulane", #Black F ND
               
               "ND_hispanic_HPAP-080_M_nPod", #Hispanic M ND
               "ND_hispanic_HPAP-099_F_UPenn", "ND_hispanic_HPAP-101_F_nPod", "ND_hispanic_HPAP-105_F_nPod", #Hispanic F ND
               
               "ND_white_HP2107001_M_Tulane", "ND_white_HP2107901_M_Tulane", "ND_white_HPAP-026_M_nPod", "ND_white_HPAP-035_M_UPenn", "ND_white_HPAP-040_M_UPenn", "ND_white_HPAP-056_M_UPenn", "ND_white_HPAP-059_M_UPenn", "ND_white_HPAP-075_M_UPenn", "ND_white_HPAP-077_M_UPenn", "ND_white_HPAP-082_M_nPod", "ND_white_SAMN15877725_M_Tulane", #White M ND
               "ND_white_HP2022801_F_Tulane", "ND_white_HP2024001_F_Tulane", "ND_white_HP2105501_F_Tulane", "ND_white_HP2108601_F_Tulane", "ND_white_HP2108901_F_Tulane", "ND_white_HPAP-022_F_UPenn", "ND_white_HPAP-036_F_nPod", "ND_white_HPAP-037_F_UPenn", "ND_white_HPAP-053_F_UPenn", "ND_white_HPAP-054_F_UPenn", "ND_white_HPAP-063_F_UPenn",  "ND_white_HPAP-074_F_UPenn", "ND_white_HPAP-103_F_UPenn", #White F ND  
               
               "T2D_black_HPAP-070_M_nPod", "T2D_black_HPAP-083_M_UPenn", "T2D_black_HPAP-108_M_nPod", #Black M T2D  
               "T2D_black_HPAP-051_F_UPenn", "T2D_black_HPAP-058_F_nPod", "T2D_black_HPAP-061_F_UPenn", #Black F T2D
               
               "T2D_hispanic_HPAP-079_F_nPod", "T2D_hispanic_HPAP-091_F_nPod", "T2D_hispanic_HPAP-109_F_nPod", #Hispanic F T2D) #Hispanic F T2D
               
               "T2D_white_HPAP-088_M_nPod", "T2D_white_HPAP-100_M_nPod", "T2D_white_HPAP-106_M_UPenn", "T2D_white_HPAP-065_M_nPod",# White M T2D
               "T2D_white_HPAP-057_F_UPenn", "T2D_white_HPAP-081_F_nPod", "T2D_white_HPAP-085_F_UPenn") # White F T2D

table(processed_rna$disease_ancestry_lib_sex_source)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)
table(processed_rna$disease_ancestry_lib_sex_source)

# Sort out ND only 
Idents(processed_rna) <- "Diabetes Status"
nd.pancreas <- subset(processed_rna, idents = c("ND"))

Idents(nd.pancreas) <- "disease_ancestry_lib_sex_source_celltype"

DefaultAssay(nd.pancreas) <- "RNA"
#DefaultAssay(nd.pancreas) <- "SCT"
#combined_processed_rna <- AverageExpression(nd.pancreas, return.seurat = TRUE, slot = 'data')
combined_processed_rna <- Seurat:::PseudobulkExpression(object = nd.pancreas, 
                                                        pb.method = 'aggregate', 
                                                        return.seurat = TRUE,
                                                        slot = 'counts')

{
  combined_processed_rna$disease_ancestry_lib_sex_source_celltype <- combined_processed_rna@active.ident
  Idents(combined_processed_rna) <- 'disease_ancestry_lib_sex_source_celltype'
  combined_processed_rna$disease <- combined_processed_rna$orig.ident
  metadat <- combined_processed_rna@meta.data
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "quiescent_stellate", "quiescent-stellate"))
  metadat <- metadat %>% 
    mutate(disease_ancestry_lib_sex_source_celltype = str_replace(disease_ancestry_lib_sex_source_celltype, "cycling_endo", "cycling-endo"))
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, "_", -5)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -4)
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -3)
  metadat$source <- metadat[c('source')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -2)
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$disease_ancestry_lib_sex_source_celltype, '_', -1)
  combined_processed_rna@meta.data = metadat
}

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("delta", "beta+delta", "beta", "beta+alpha", "alpha", "gamma", "epsilon", "cycling-endo",
               "ductal", "acinar",
               "activated-stellate", "quiescent-stellate", "endothelial",
               "lymphocyte", "macrophages", "mast", "schwann") 

table(combined_processed_rna$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_rna$celltype <- factor(x = combined_processed_rna$celltype, levels = my_levels)
table(combined_processed_rna$celltype)
Idents(combined_processed_rna) <- "celltype"
    
# DE testing to determine celltype specificity
DefaultAssay(nd.pancreas) <- "RNA"
Idents(combined_processed_rna) <- "celltype"
combined_processed_rna[["RNA"]]@counts<-as.matrix(combined_processed_rna[["RNA"]]@counts)+1

plan(strategy = "multicore", workers = 80)
beta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 1, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "beta", 
                                      group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                      only.pos = TRUE)
beta.conserved.markers <- dplyr::filter(beta.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 1)


alpha.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "alpha", 
                                       group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
alpha.conserved.markers <- dplyr::filter(alpha.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


delta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "delta", 
                                       group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
delta.conserved.markers <- dplyr::filter(delta.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


gamma.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "gamma", 
                                       group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
gamma.conserved.markers <- dplyr::filter(gamma.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


epsilon.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "epsilon", 
                                         group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                         only.pos = TRUE)
epsilon.conserved.markers <- dplyr::filter(epsilon.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


betaalpha.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "betaalpha", 
                                           group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                           only.pos = TRUE)
betaalpha.conserved.markers <- dplyr::filter(betaalpha.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


betadelta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "betadelta", 
                                           group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                           only.pos = TRUE)
betadelta.conserved.markers <- dplyr::filter(betadelta.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


cycling_endo.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "cycling_endo", 
                                              group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                              only.pos = TRUE)
cycling_endo.conserved.markers <- dplyr::filter(cycling_endo.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


acinar.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "acinar", 
                                        group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                        only.pos = TRUE)
acinar.conserved.markers <- dplyr::filter(acinar.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


ductal.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "ductal", 
                                        group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                        only.pos = TRUE)
ductal.conserved.markers <- dplyr::filter(ductal.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


activated_stellate.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "activated_stellate", 
                                                    group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                                    only.pos = TRUE)
activated_stellate.conserved.markers <- dplyr::filter(activated_stellate.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


quiescent_stellate.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "quiescent_stellate", 
                                                    group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                                    only.pos = TRUE)
quiescent_stellate.conserved.markers <- dplyr::filter(quiescent_stellate.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


endothelial.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "endothelial", 
                                             group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                             only.pos = TRUE)
endothelial.conserved.markers <- dplyr::filter(endothelial.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


lymphocyte.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "lymphocyte", 
                                            group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                            only.pos = TRUE)
lymphocyte.conserved.markers <- dplyr::filter(lymphocyte.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


macrophages.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "macrophages", 
                                             group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                             only.pos = TRUE)
macrophages.conserved.markers <- dplyr::filter(macrophages.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


mast.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "mast", 
                                      group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                      only.pos = TRUE)
mast.conserved.markers <- dplyr::filter(mast.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 


schwann.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "schwann", 
                                         group.by = "celltype", min.cells.group = 1, pseudocount.use = 1, test.use = "DESeq2", 
                                         only.pos = TRUE)
schwann.conserved.markers <- dplyr::filter(schwann.conserved.markers, p_val_adj < 0.001 & avg_log2FC >= 1) 










Idents(processed_rna) <- "ancestry_sex"
    processed_rna$ancestry_sex_diabetes <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes"]])
    
    Idents(processed_rna) <- "ancestry_sex_diabetes"
    processed_rna$ancestry_sex_diabetes_sample <- paste(Idents(processed_rna), processed_rna$'Library', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes_sample"]])
    
    Idents(processed_rna) <- "ancestry_sex_diabetes_sample"
    processed_rna$ancestry_sex_diabetes_sample_origin <- paste(Idents(processed_rna), processed_rna$'Tissue Source', sep = "_")
    table(processed_rna@meta.data[["ancestry_sex_diabetes_sample_origin"]])
    
    
    # Dotplots
    Idents(processed_rna) <- "Diabetes Status"
    ND <- subset(processed_rna, idents = c('ND'))
    T2D <- subset(processed_rna, idents = c('T2D'))
    ND<- NULL
    bulk <- NULL
    genes.to.plot <- c('XIST', 'SRY')
    
    DefaultAssay(ND) <- 'SCT'
    Idents(ND) <- 'celltype_sex'
    table(ND$celltype_sex)
    # Define an order of cluster identities remember after this step-
    # cluster re-assignment occurs, which re-assigns clustering in my_levels
    my_levels <- c("alpha_F", 
                   "beta+alpha_F",
                   "beta_F", 
                   "beta+delta_F",
                   "delta_F",
                   "gamma_F",
                   "cycling_endo_F",
                   "epsilon_F",
                   "ductal_F",
                   "acinar_F",
                   "quiescent_stellate_F",
                   "activated_stellate_F",
                   "schwann_F",
                   "mast_F",
                   "endothelial_F",
                   "lymphocyte_F",
                   "macrophages_F",
                   "alpha_M", 
                   "beta+alpha_M",
                   "beta_M", 
                   "beta+delta_M", 
                   "delta_M", 
                   "gamma_M", 
                   "cycling_endo_M",
                   "epsilon_M",
                   "ductal_M",
                   "acinar_M",
                   "quiescent_stellate_M",
                   "activated_stellate_M",
                   "schwann_M",
                   "mast_M",
                   "endothelial_M",
                   "lymphocyte_M",
                   "macrophages_M")
    head(ND@meta.data$celltype_sex)
    
    # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
    ND@meta.data$celltype_sex <- factor(x = ND@meta.data$celltype_sex, levels = my_levels)
    Idents(ND) <- "celltype_sex"
    
    # Selected genes
    markers.to.plot <- c("DDX3Y", "EIF1AY",
                         "KDM5D", "NLGN4Y",
                         "RPS4Y1","USP9Y", 
                         "UTY", "ZFY",
                         "XIST", "TSIX",
                         "ZFX", "KDM5C",
                         "SEPTIN6", "EIF1AX",
                         "KDM6A", "PUDP", "DDX3X")
    
    # Dotplot
    DotPlot(ND,  
            dot.scale = 8,
            col.min = -1, #minimum level
            col.max = 1,  #maximum level
            features = rev(markers.to.plot),
            scale = TRUE) + 
      geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
      theme_light() +
      #facet_wrap(~??? what metadata should be here??)
      #coord_flip() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.title=element_text(size=12, face = "bold"), 
            legend.text=element_text(size=12, face = "bold")) +
      scale_colour_gradient2(low =c("dodgerblue"), high =c("red3")) +
      guides(color = guide_colorbar(title = 'Average Expression')) + coord_flip() 
    
    # Dotplot
    DotPlot(ND,  
            dot.scale = 1,
            col.min = -1, #minimum level
            col.max = 1,  #maximum level
            features = rev(x.chrom),
            scale = TRUE) + 
      geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
      theme_light() +
      #facet_wrap(~??? what metadata should be here??)
      #coord_flip() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.title=element_text(size=12, face = "bold"), 
            legend.text=element_text(size=12, face = "bold")) +
      scale_colour_gradient2(low =c("dodgerblue"), high =c("red3")) +
      guides(color = guide_colorbar(title = 'Average Expression')) + coord_flip() 
    # 
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    x.chrom <- c("ABCB7",
                 "ABCD1",
                 "ACE2",
                 "ACOT9",
                 "ACSL4",
                 "ACTRT1",
                 "ADGRG2",
                 "ADGRG4",
                 "AFF2",
                 "AGTR2",
                 "AIFM1",
                 "AKAP4",
                 "AKAP14",
                 "AKAP17A",
                 "ALAS2",
                 "ALG13",
                 "AMELX",
                 "AMER1",
                 "AMMECR1",
                 "AMOT",
                 "ANOS1",
                 "AP1S2",
                 "APEX2",
                 "APLN",
                 "APOO",
                 "APOOL",
                 "AR",
                 "ARAF",
                 "ARHGAP4",
                 "ARHGAP6",
                 "ARHGAP36",
                 "ARHGEF6",
                 "ARHGEF9",
                 "ARL13A",
                 "ARMCX1",
                 "ARMCX2",
                 "ARMCX3",
                 "ARMCX4",
                 "ARMCX5",
                 "ARMCX6",
                 "ARR3",
                 "ARSD",
                 "ARSF",
                 "ARSH",
                 "ARSL",
                 "ARX",
                 "ASB9",
                 "ASB11",
                 "ASB12",
                 "ASMT",
                 "ASMTL",
                 "ATG4A",
                 "ATP1B4",
                 "ATP2B3",
                 "ATP6AP1",
                 "ATP6AP2",
                 "ATP7A",
                 "ATP11C",
                 "ATRX",
                 "ATXN3L",
                 "AVPR2",
                 "AWAT1",
                 "AWAT2",
                 "BCAP31",
                 "BCLAF3",
                 "BCOR",
                 "BCORL1",
                 "BEND2",
                 "BEX1",
                 "BEX2",
                 "BEX3",
                 "BEX4",
                 "BEX5",
                 "BGN",
                 "BMP15",
                 "BMX",
                 "BRCC3",
                 "BRS3",
                 "BRWD3",
                 "BTK",
                 "C1GALT1C1",
                 "CA5B",
                 "CACNA1F",
                 "CAPN6",
                 "CASK",
                 "CBLL2",
                 "CCDC22",
                 "CCDC120",
                 "CCDC160",
                 "CCNB3",
                 "CCNQ",
                 "CD40LG",
                 "CD99",
                 "CD99L2",
                 "CDK16",
                 "CDKL5",
                 "CDX4",
                 "CENPI",
                 "CENPVL1",
                 "CENPVL2",
                 "CENPVL3",
                 "CETN2",
                 "CFAP47",
                 "CFP",
                 "CHIC1",
                 "CHM",
                 "CHRDL1",
                 "CHST7",
                 "CITED1",
                 "CLCN4",
                 "CLCN5",
                 "CLDN2",
                 "CLDN34",
                 "CLIC2",
                 "CLTRN",
                 "CMC4",
                 "CNGA2",
                 "CNKSR2",
                 "COL4A5",
                 "COL4A6",
                 "COX7B",
                 "CPXCR1",
                 "CRLF2",
                 "CSAG1",
                 "CSAG2",
                 "CSAG3",
                 "CSF2RA",
                 "CSTF2",
                 "CT45A1",
                 "CT45A2",
                 "CT45A3",
                 "CT45A5",
                 "CT45A6",
                 "CT45A7",
                 "CT45A8",
                 "CT45A9",
                 "CT45A10",
                 "CT47A1",
                 "CT47A2",
                 "CT47A3",
                 "CT47A4",
                 "CT47A5",
                 "CT47A6",
                 "CT47A7",
                 "CT47A8",
                 "CT47A9",
                 "CT47A10",
                 "CT47A11",
                 "CT47A12",
                 "CT47B1",
                 "CT47C1",
                 "CT55",
                 "CT83",
                 "CTAG1A",
                 "CTAG1B",
                 "CTAG2",
                 "CTPS2",
                 "CUL4B",
                 "CXCR3",
                 "CXorf38",
                 "CXorf49",
                 "CXorf49B",
                 "CXorf51A",
                 "CXorf51B",
                 "CXorf58",
                 "CXorf65",
                 "CXorf66",
                 "CYBB",
                 "CYLC1",
                 "CYSLTR1",
                 "DACH2",
                 "DCAF8L1",
                 "DCAF8L2",
                 "DCAF12L1",
                 "DCAF12L2",
                 "DCX",
                 "DDX3X",
                 "DDX53",
                 "DGAT2L6",
                 "DGKK",
                 "DHRSX",
                 "DIAPH2",
                 "DIPK2B",
                 "DKC1",
                 "DLG3",
                 "DMD",
                 "DMRTC1",
                 "DMRTC1B",
                 "DNAAF6",
                 "DNASE1L1",
                 "DOCK11",
                 "DRP2",
                 "DUSP9",
                 "DUSP21",
                 "DYNLT3",
                 "EBP",
                 "EDA",
                 "EDA2R",
                 "EFHC2",
                 "EFNB1",
                 "EGFL6",
                 "EIF1AX",
                 "EIF2S3",
                 "ELF4",
                 "ELK1",
                 "EMD",
                 "ENOX2",
                 "EOLA1",
                 "EOLA2",
                 "ERAS",
                 "ERCC6L",
                 "ESX1",
                 "ETDA",
                 "ETDB",
                 "ETDC",
                 "EZHIP",
                 "F8",
                 "F8A1",
                 "F8A2",
                 "F8A3",
                 "F9",
                 "FAAH2",
                 "FAM3A",
                 "FAM9A",
                 "FAM9B",
                 "FAM9C",
                 "FAM47A",
                 "FAM47B",
                 "FAM47C",
                 "FAM50A",
                 "FAM104B",
                 "FAM120C",
                 "FAM133A",
                 "FAM156A",
                 "FAM156B",
                 "FAM199X",
                 "FAM236A",
                 "FAM236B",
                 "FAM236C",
                 "FAM236D",
                 "FANCB",
                 "FATE1",
                 "FGD1",
                 "FGF13",
                 "FGF16",
                 "FHL1",
                 "FLNA",
                 "FMR1",
                 "FMR1NB",
                 "FOXO4",
                 "FOXP3",
                 "FOXR2",
                 "FRMD7",
                 "FRMPD3",
                 "FRMPD4",
                 "FTHL17",
                 "FTSJ1",
                 "FUNDC1",
                 "FUNDC2",
                 "G6PD",
                 "GAB3",
                 "GABRA3",
                 "GABRE",
                 "GABRQ",
                 "GAGE1",
                 "GAGE2A",
                 "GAGE2B",
                 "GAGE2C",
                 "GAGE2D",
                 "GAGE2E",
                 "GAGE4",
                 "GAGE5",
                 "GAGE6",
                 "GAGE7",
                 "GAGE8",
                 "GAGE10",
                 "GAGE12B",
                 "GAGE12C",
                 "GAGE12D",
                 "GAGE12E",
                 "GAGE12F",
                 "GAGE12G",
                 "GAGE12H",
                 "GAGE12I",
                 "GAGE12J",
                 "GAGE13",
                 "GATA1",
                 "GCNA",
                 "GDI1",
                 "GDPD2",
                 "GEMIN8",
                 "GJB1",
                 "GK",
                 "GLA",
                 "GLOD5",
                 "GLRA2",
                 "GLUD2",
                 "GNG5B",
                 "GNL3L",
                 "GPC3",
                 "GPC4",
                 "GPKOW",
                 "GPM6B",
                 "GPR34",
                 "GPR50",
                 "GPR82",
                 "GPR101",
                 "GPR119",
                 "GPR143",
                 "GPR173",
                 "GPR174",
                 "GPRASP1",
                 "GPRASP2",
                 "GPRASP3",
                 "GRIA3",
                 "GRIPAP1",
                 "GRPR",
                 "GSPT2",
                 "GTPBP6",
                 "GUCY2F",
                 "GYG2",
                 "H2AB1",
                 "H2AB2",
                 "H2AB3",
                 "H2AL3",
                 "H2AP",
                 "H2BW1",
                 "H2BW2",
                 "HAPSTR2",
                 "HAUS7",
                 "HCCS",
                 "HCFC1",
                 "HDAC6",
                 "HDAC8",
                 "HDX",
                 "HEPH",
                 "HMGB3",
                 "HMGN5",
                 "HNRNPH2",
                 "HPRT1",
                 "HS6ST2",
                 "HSD17B10",
                 "HSFX1",
                 "HSFX2",
                 "HSFX3",
                 "HSFX4",
                 "HTATSF1",
                 "HTR2C",
                 "HUWE1",
                 "IDH3G",
                 "IDS",
                 "IGBP1",
                 "IGSF1",
                 "IKBKG",
                 "IL1RAPL1",
                 "IL1RAPL2",
                 "IL2RG",
                 "IL3RA",
                 "IL9R",
                 "IL13RA1",
                 "IL13RA2",
                 "INTS6L",
                 "IQSEC2",
                 "IRAK1",
                 "IRS4",
                 "ITGB1BP2",
                 "ITIH6",
                 "ITM2A",
                 "JADE3",
                 "KANTR",
                 "KCND1",
                 "KCNE5",
                 "KDM5C",
                 "KDM6A",
                 "KIAA1210",
                 "KIF4A",
                 "KLF8",
                 "KLHL4",
                 "KLHL13",
                 "KLHL15",
                 "KLHL34",
                 "KRBOX4",
                 "L1CAM",
                 "LAGE3",
                 "LAMP2",
                 "LANCL3",
                 "LAS1L",
                 "LDOC1",
                 "LHFPL1",
                 "LONRF3",
                 "LPAR4",
                 "LRCH2",
                 "LUZP4",
                 "MAGEA1",
                 "MAGEA2",
                 "MAGEA2B",
                 "MAGEA3",
                 "MAGEA4",
                 "MAGEA6",
                 "MAGEA8",
                 "MAGEA9",
                 "MAGEA9B",
                 "MAGEA10",
                 "MAGEA11",
                 "MAGEA12",
                 "MAGEB1",
                 "MAGEB2",
                 "MAGEB3",
                 "MAGEB4",
                 "MAGEB5",
                 "MAGEB6",
                 "MAGEB6B",
                 "MAGEB10",
                 "MAGEB16",
                 "MAGEB17",
                 "MAGEB18",
                 "MAGEC1",
                 "MAGEC2",
                 "MAGEC3",
                 "MAGED1",
                 "MAGED2",
                 "MAGED4",
                 "MAGED4B",
                 "MAGEE1",
                 "MAGEE2",
                 "MAGEH1",
                 "MAGIX",
                 "MAGT1",
                 "MAMLD1",
                 "MAOA",
                 "MAOB",
                 "MAP3K15",
                 "MAP7D2",
                 "MAP7D3",
                 "MBNL3",
                 "MBTPS2",
                 "MCF2",
                 "MCTS1",
                 "MECP2",
                 "MED12",
                 "MED14",
                 "MED14OS",
                 "MID1",
                 "MID1IP1",
                 "MID2",
                 "MMGT1",
                 "MORC4",
                 "MORF4L2",
                 "MOSPD1",
                 "MOSPD2",
                 "MPC1L",
                 "MPP1",
                 "MSL3",
                 "MSN",
                 "MTCP1",
                 "MTM1",
                 "MTMR1",
                 "MTMR8",
                 "MXRA5",
                 "NAA10",
                 "NALF2",
                 "NAP1L2",
                 "NAP1L3",
                 "NBDY",
                 "NCBP2L",
                 "NDP",
                 "NDUFA1",
                 "NDUFB11",
                 "NEXMIF",
                 "NHS",
                 "NHSL2",
                 "NKAP",
                 "NKRF",
                 "NLGN3",
                 "NLGN4X",
                 "NLRP2B",
                 "NONO",
                 "NOX1",
                 "NR0B1",
                 "NRK",
                 "NSDHL",
                 "NUDT10",
                 "NUDT11",
                 "NUP62CL",
                 "NXF2",
                 "NXF2B",
                 "NXF3",
                 "NXF5",
                 "NXT2",
                 "NYX",
                 "OCRL",
                 "OFD1",
                 "OGT",
                 "OPHN1",
                 "OPN1LW",
                 "OPN1MW",
                 "OPN1MW2",
                 "OPN1MW3",
                 "OR13H1",
                 "OTC",
                 "OTUD5",
                 "OTUD6A",
                 "P2RY4",
                 "P2RY8",
                 "P2RY10",
                 "PABIR2",
                 "PABIR3",
                 "PABPC1L2A",
                 "PABPC1L2B",
                 "PABPC5",
                 "PAGE1",
                 "PAGE2",
                 "PAGE2B",
                 "PAGE3",
                 "PAGE4",
                 "PAGE5",
                 "PAK3",
                 "PASD1",
                 "PBDC1",
                 "PCDH11X",
                 "PCDH19",
                 "PCSK1N",
                 "PCYT1B",
                 "PDHA1",
                 "PDK3",
                 "PDZD4",
                 "PDZD11",
                 "PFKFB1",
                 "PGAM4",
                 "PGK1",
                 "PGRMC1",
                 "PHEX",
                 "PHF6",
                 "PHF8",
                 "PHKA1",
                 "PHKA2",
                 "PIGA",
                 "PIM2",
                 "PIN4",
                 "PIR",
                 "PJA1",
                 "PLAC1",
                 "PLCXD1",
                 "PLP1",
                 "PLP2",
                 "PLS3",
                 "PLXNA3",
                 "PLXNB3",
                 "PNCK",
                 "PNMA3",
                 "PNMA5",
                 "PNMA6A",
                 "PNMA6E",
                 "PNMA6F",
                 "PNPLA4",
                 "POF1B",
                 "POLA1",
                 "PORCN",
                 "POU3F4",
                 "PPEF1",
                 "PPP1R2C",
                 "PPP1R3F",
                 "PPP2R3B",
                 "PPP4R3C",
                 "PQBP1",
                 "PRAF2",
                 "PRDX4",
                 "PRICKLE3",
                 "PRKX",
                 "PRPS1",
                 "PRPS2",
                 "PRR32",
                 "PRRG1",
                 "PRRG3",
                 "PSMD10",
                 "PTCHD1",
                 "PUDP",
                 "PWWP3B",
                 "PWWP4",
                 "RAB9A",
                 "RAB9B",
                 "RAB33A",
                 "RAB39B",
                 "RAB40A",
                 "RAB40AL",
                 "RAB41",
                 "RADX",
                 "RAI2",
                 "RAP2C",
                 "RBBP7",
                 "RBM3",
                 "RBM10",
                 "RBM41",
                 "RBMX",
                 "RBMX2",
                 "RBMXL3",
                 "RENBP",
                 "REPS2",
                 "RGN",
                 "RHOXF1",
                 "RHOXF2",
                 "RHOXF2B",
                 "RIBC1",
                 "RIPPLY1",
                 "RLIM",
                 "RNF113A",
                 "RNF128",
                 "RP2",
                 "RPA4",
                 "RPGR",
                 "RPL10",
                 "RPL36A",
                 "RPL39",
                 "RPS4X",
                 "RPS6KA3",
                 "RPS6KA6",
                 "RRAGB",
                 "RS1",
                 "RTL3",
                 "RTL4",
                 "RTL5",
                 "RTL8A",
                 "RTL8B",
                 "RTL8C",
                 "RTL9",
                 "S100G",
                 "SAGE1",
                 "SASH3",
                 "SAT1",
                 "SATL1",
                 "SCML1",
                 "SCML2",
                 "SEPTIN6",
                 "SERPINA7",
                 "SERTM2",
                 "SH2D1A",
                 "SH3BGRL",
                 "SH3KBP1",
                 "SHOX",
                 "SHROOM2",
                 "SHROOM4",
                 "SLC6A8",
                 "SLC6A14",
                 "SLC7A3",
                 "SLC9A6",
                 "SLC9A7",
                 "SLC10A3",
                 "SLC16A2",
                 "SLC25A5",
                 "SLC25A6",
                 "SLC25A14",
                 "SLC25A43",
                 "SLC25A53",
                 "SLC35A2",
                 "SLC38A5",
                 "SLITRK2",
                 "SLITRK4",
                 "SMARCA1",
                 "SMC1A",
                 "SMIM9",
                 "SMIM10",
                 "SMIM10L2A",
                 "SMIM10L2B",
                 "SMPX",
                 "SMS",
                 "SNX12",
                 "SOWAHD",
                 "SOX3",
                 "SPACA5",
                 "SPACA5B",
                 "SPANXA1",
                 "SPANXA2",
                 "SPANXB1",
                 "SPANXC",
                 "SPANXD",
                 "SPANXN1",
                 "SPANXN2",
                 "SPANXN3",
                 "SPANXN4",
                 "SPANXN5",
                 "SPIN2A",
                 "SPIN2B",
                 "SPIN3",
                 "SPIN4",
                 "SPRY3",
                 "SRPK3",
                 "SRPX",
                 "SRPX2",
                 "SSR4",
                 "SSX1",
                 "SSX2",
                 "SSX2B",
                 "SSX3",
                 "SSX4",
                 "SSX4B",
                 "SSX5",
                 "SSX7",
                 "STAG2",
                 "STARD8",
                 "STEEP1",
                 "STK26",
                 "STS",
                 "SUPT20HL1",
                 "SUPT20HL2",
                 "SUV39H1",
                 "SYAP1",
                 "SYN1",
                 "SYP",
                 "SYTL4",
                 "SYTL5",
                 "TAB3",
                 "TAF1",
                 "TAF7L",
                 "TAF9B",
                 "TAFAZZIN",
                 "TASL",
                 "TBC1D8B",
                 "TBC1D25",
                 "TBL1X",
                 "TBX22",
                 "TCEAL1",
                 "TCEAL2",
                 "TCEAL3",
                 "TCEAL4",
                 "TCEAL5",
                 "TCEAL6",
                 "TCEAL7",
                 "TCEAL8",
                 "TCEAL9",
                 "TCEANC",
                 "TCP11X1",
                 "TCP11X2",
                 "TENM1",
                 "TENT5D",
                 "TEX11",
                 "TEX13A",
                 "TEX13B",
                 "TEX13C",
                 "TEX13D",
                 "TEX28",
                 "TFDP3",
                 "TFE3",
                 "TGIF2LX",
                 "THOC2",
                 "TIMM8A",
                 "TIMM17B",
                 "TIMP1",
                 "TKTL1",
                 "TLR7",
                 "TLR8",
                 "TMEM31",
                 "TMEM35A",
                 "TMEM47",
                 "TMEM164",
                 "TMEM185A",
                 "TMEM187",
                 "TMEM255A",
                 "TMLHE",
                 "TMSB4X",
                 "TMSB15A",
                 "TMSB15B",
                 "TMSB15C",
                 "TNMD",
                 "TRAPPC2",
                 "TREX2",
                 "TRMT2B",
                 "TRO",
                 "TRPC5",
                 "TRPC5OS",
                 "TSC22D3",
                 "TSPAN6",
                 "TSPAN7",
                 "TSPYL2",
                 "TSR2",
                 "TXLNG",
                 "UBA1",
                 "UBE2A",
                 "UBE2NL",
                 "UBL4A",
                 "UBQLN2",
                 "UPF3B",
                 "UPRT",
                 "USP9X",
                 "USP11",
                 "USP26",
                 "USP27X",
                 "USP51",
                 "UTP14A",
                 "UXT",
                 "VAMP7",
                 "VBP1",
                 "VCX",
                 "VCX2",
                 "VCX3A",
                 "VCX3B",
                 "VEGFD",
                 "VGLL1",
                 "VMA21",
                 "VSIG1",
                 "VSIG4",
                 "WAS",
                 "WDR13",
                 "WDR44",
                 "WDR45",
                 "WNK3",
                 "WWC3",
                 "XAGE1A",
                 "XAGE1B",
                 "XAGE2",
                 "XAGE3",
                 "XAGE5",
                 "XG",
                 "XIAP",
                 "XK",
                 "XKRX",
                 "XPNPEP2",
                 "YIPF6",
                 "YY2",
                 "ZBED1",
                 "ZBTB33",
                 "ZC3H12B",
                 "ZC4H2",
                 "ZCCHC12",
                 "ZCCHC13",
                 "ZCCHC18",
                 "ZDHHC9",
                 "ZDHHC15",
                 "ZFP92",
                 "ZFX",
                 "ZIC3",
                 "ZMAT1",
                 "ZMYM3",
                 "ZNF41",
                 "ZNF75D",
                 "ZNF81",
                 "ZNF157",
                 "ZNF182",
                 "ZNF185",
                 "ZNF275",
                 "ZNF280C",
                 "ZNF449",
                 "ZNF630",
                 "ZNF674",
                 "ZNF711",
                 "ZRSR2",
                 "ZXDA",
                 "ZXDB")
    
    BiocManager::install("scuttle", force = TRUE)
    library(scuttle)
    remotes::install_github("omnideconv/SimBu")
    library(SimBu)
    # View clustering
    DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
    pancreas_rna <- FindClusters(object = pancreas_rna, algorithm=4, resolution = 0, method = 'igraph')
    
    
    Idents(pancreas_rna) <- "Chemistry"
    
    # QC for Ruth split on chemistry
    v2 <- subset(pancreas_rna, idents = c("10Xv2"))
    v3 <- subset(pancreas_rna, idents = c("10Xv3"))
    
    Idents(pancreas_rna) <- "Chemistry"
    DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
    
    p1 <- DimPlot(v2, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
    p2 <- DimPlot(v3, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
    p1|p2
    
    RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE, reduction.name = 'umap')
    DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
    DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
