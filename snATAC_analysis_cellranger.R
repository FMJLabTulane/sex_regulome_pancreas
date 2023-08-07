# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R
# make sure a Rtools version is installed that supports your current versin of R
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
BiocManager::install(version = "3.17")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), force = TRUE)

BiocManager::install(c("EnhancedVolcano", "DoubletFinder", "glmGamPoi",
                       "GOSemSim", "org.Hs.eg.db", "AnnotationHub",
                       "GenomeInfoDb", "MeSHDbi", "clusterProfiler",
                       "dittoSeq", "escape", "ComplexHeatmap", "DropletUtils", 
                       "Nebulosa", "hdf5r", "scDblFinder", "JASPAR2020",
                       "TFBSTools", "motifmatchr", "GreenleafLab/chromVAR",
                       "EnrichmentBrowser"),
                     dependencies = T, force = TRUE)
BiocManager::install("BiocParallel")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

# Run the following code once you have Seurat installed
suppressWarnings({
    library(ggplot2)
    library(ggrepel)
    library(leiden)
    library(stringr)
    library(hdf5r)
    library(SoupX)
    library(Rcpp)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(dplyr)
    library(tidyverse)
    library(data.table)
    library(reticulate)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(Signac)
    library(scDblFinder)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(plotly)
    library(clustree)
    library(patchwork)
    #library(future)
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
    library(scCustomize)
    library(Nebulosa)
    library(DropletUtils)
    library(ggvenn)
    library(ggVennDiagram)
    library(devtools)
    library(R.utils)
    library(qs)
    library(JASPAR2020)
    library(TFBSTools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(motifmatchr)
    library(chromVAR)
    library(SeuratWrappers)
    library(cicero)
    library(BiocParallel)
})

# Set global environment parameter par-proc
# options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())


# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")
packageVersion("monocle3")
packageVersion("cicero")

############################ STAGE ############################
############################   1   ############################
# Calculation of Doublets https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scATAC.html
# Provide GRanges of repeat elements for exclusion:
suppressPackageStartupMessages(library(GenomicRanges))
repeats <- GRanges("chr6", IRanges(1000,2000))

# Combine with mitochondrial and sex chromosomes
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))

# Combining them:
toExclude <- suppressWarnings(c(repeats, otherChroms))

# Running amulet method
{
  fragfile.HP2022801 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.SAMN15877725 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2024001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2031401 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2105501 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2106201 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2107001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2107901 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2108601 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2108901 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2110001 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2121601 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2123201 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2132801 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)", regionsToExclude=toExclude)
  fragfile.HP2202101 <- amulet(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)", regionsToExclude=toExclude)
}

############################ STAGE ############################
############################   2   ############################
# SAMPLE PRE_PROCESSING
# Adding sampleID to barcode nomenclature
head(combined_atac) # shows that config of cell nomenclature is SAMPLEID_barcode [yeah I know, this shouldnt but hey live dangerous.]
{
  rownames(fragfile.HP2022801) <- paste0('HP2022801_', rownames(fragfile.HP2022801))
  rownames(fragfile.SAMN15877725) <- paste0('SAMN15877725_', rownames(fragfile.SAMN15877725))
  rownames(fragfile.HP2024001) <- paste0('HP2024001_', rownames(fragfile.HP2024001))
  rownames(fragfile.HP2031401) <- paste0('HP2031401_', rownames(fragfile.HP2031401))
  rownames(fragfile.HP2105501) <- paste0('HP2105501_', rownames(fragfile.HP2105501))
  rownames(fragfile.HP2106201) <- paste0('HP2106201_', rownames(fragfile.HP2106201))
  rownames(fragfile.HP2107001) <- paste0('HP2107001_', rownames(fragfile.HP2107001))
  rownames(fragfile.HP2107901) <- paste0('HP2107901_', rownames(fragfile.HP2107901))
  rownames(fragfile.HP2108601) <- paste0('HP2108601_', rownames(fragfile.HP2108601))
  rownames(fragfile.HP2108901) <- paste0('HP2108901_', rownames(fragfile.HP2108901))
  rownames(fragfile.HP2110001) <- paste0('HP2110001_', rownames(fragfile.HP2110001))
  rownames(fragfile.HP2121601) <- paste0('HP2121601_', rownames(fragfile.HP2121601))
  rownames(fragfile.HP2123201) <- paste0('HP2123201_', rownames(fragfile.HP2123201))
  rownames(fragfile.HP2132801) <- paste0('HP2132801_', rownames(fragfile.HP2132801))
  rownames(fragfile.HP2202101) <- paste0('HP2202101_', rownames(fragfile.HP2202101))
}

combined_atac_doublet <- do.call('rbind', list(fragfile.HP2022801, fragfile.SAMN15877725, fragfile.HP2024001, fragfile.HP2031401,
                                               fragfile.HP2105501, fragfile.HP2106201, fragfile.HP2107001, fragfile.HP2107901,
                                               fragfile.HP2108601, fragfile.HP2108901, fragfile.HP2110001, fragfile.HP2121601,
                                               fragfile.HP2123201, fragfile.HP2132801, fragfile.HP2202101))

#Save file
# saveRDS(fragfile.HP2022801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2022801.rds)")
# saveRDS(fragfile.SAMN15877725, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.SAMN15877725.rds)")
# saveRDS(fragfile.HP2024001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2024001.rds)")
# saveRDS(fragfile.HP2031401, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2031401.rds)")
# saveRDS(fragfile.HP2105501, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2105501.rds)")
# saveRDS(fragfile.HP2106201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2106201.rds)")
# saveRDS(fragfile.HP2107001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2107001.rds)")
# saveRDS(fragfile.HP2107901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2107901.rds)")
# saveRDS(fragfile.HP2108601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2108601.rds)")
# saveRDS(fragfile.HP2108901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2108901.rds)")
# saveRDS(fragfile.HP2110001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2110001.rds)")
# saveRDS(fragfile.HP2121601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2121601.rds)")
# saveRDS(fragfile.HP2123201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2123201.rds)")
# saveRDS(fragfile.HP2132801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2132801.rds)")
# saveRDS(fragfile.HP2202101, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2202101.rds)")
# saveRDS(combined_atac_doublet, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\combined_atac_doublet.rds)")

# Load
# fragfile.HP2022801 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2022801.rds)")
# fragfile.SAMN15877725 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.SAMN15877725.rds)")
# fragfile.HP2024001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2024001.rds)")
# fragfile.HP2031401 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2031401.rds)")
# fragfile.HP2105501 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2105501.rds)")
# fragfile.HP2106201 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2106201.rds)")
# fragfile.HP2107001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2107001.rds)")
# fragfile.HP2107901 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2107901.rds)")
# fragfile.HP2108601 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2108601.rds)")
# fragfile.HP2108901 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2108901.rds)")
# fragfile.HP2110001 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2110001.rds)")
# fragfile.HP2121601 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2121601.rds)")
# fragfile.HP2123201 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2123201.rds)")
# fragfile.HP2132801 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\snATACseq\fragfile.HP2132801.rds)")
# fragfile.HP2202101 <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\fragfile.HP2202101.rds)")
combined_atac_doublet <- readRDS(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

# OBJECT SETUP AND NORMALIZATION #
# STEP 1: Load 10X data #### https://stuartlab.org/signac/articles/merging.html
{
  # read in peak sets
  peaks.HP2022801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.SAMN15877725 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2024001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2031401 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2105501 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2106201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2107001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2107901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2108601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2108901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2110001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2121601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2123201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2132801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
  
  peaks.HP2202101 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\peaks.bed)",
    col.names = c("chr", "start", "end")
  )
}

# Conversion of peaks to genomic Ranges
{
  gr.HP2022801 <- makeGRangesFromDataFrame(peaks.HP2022801)
  gr.SAMN15877725 <- makeGRangesFromDataFrame(peaks.SAMN15877725)
  gr.HP2024001 <- makeGRangesFromDataFrame(peaks.HP2024001)
  gr.HP2031401 <- makeGRangesFromDataFrame(peaks.HP2031401)
  gr.HP2105501 <- makeGRangesFromDataFrame(peaks.HP2105501)
  gr.HP2106201 <- makeGRangesFromDataFrame(peaks.HP2106201)
  gr.HP2107001 <- makeGRangesFromDataFrame(peaks.HP2107001)
  gr.HP2107901 <- makeGRangesFromDataFrame(peaks.HP2107901)
  gr.HP2108601 <- makeGRangesFromDataFrame(peaks.HP2108601)
  gr.HP2108901 <- makeGRangesFromDataFrame(peaks.HP2108901)
  gr.HP2110001 <- makeGRangesFromDataFrame(peaks.HP2110001)
  gr.HP2121601 <- makeGRangesFromDataFrame(peaks.HP2121601)
  gr.HP2123201 <- makeGRangesFromDataFrame(peaks.HP2123201)
  gr.HP2132801 <- makeGRangesFromDataFrame(peaks.HP2132801)
  gr.HP2202101 <- makeGRangesFromDataFrame(peaks.HP2202101)
  
  # Create a unified set of peaks to quantify in each dataset
  combined.peaks <- reduce(x = c(gr.HP2022801, gr.SAMN15877725, gr.HP2024001, gr.HP2031401,
                                 gr.HP2105501, gr.HP2106201, gr.HP2107001, gr.HP2107901,
                                 gr.HP2108601, gr.HP2108901, gr.HP2110001, gr.HP2121601,
                                 gr.HP2123201, gr.HP2132801, gr.HP2202101))
  
  # Filter out bad peaks based on length
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  combined.peaks
}

# Create Fragment objects
# Load metadata
{
  md.HP2022801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.SAMN15877725 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2024001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2031401 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2105501 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2106201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2107001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2107901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2108601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2108901 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2110001 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2121601 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2123201 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2132801 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  md.HP2202101 <- read.table(
    file = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\singlecell.csv)",
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
}

# perform an initial filtering of low count cells de-trash the data
md.HP2022801 <- md.HP2022801[md.HP2022801$passed_filters > 500, ]
md.SAMN15877725 <- md.SAMN15877725[md.SAMN15877725$passed_filters > 500, ]
md.HP2024001 <- md.HP2024001[md.HP2024001$passed_filters > 500, ]
md.HP2031401 <- md.HP2031401[md.HP2031401$passed_filters > 500, ]
md.HP2105501 <- md.HP2105501[md.HP2105501$passed_filters > 500, ]
md.HP2106201 <- md.HP2106201[md.HP2106201$passed_filters > 500, ]
md.HP2107001 <- md.HP2107001[md.HP2107001$passed_filters > 500, ]
md.HP2107901 <- md.HP2107901[md.HP2107901$passed_filters > 500, ]
md.HP2108601 <- md.HP2108601[md.HP2108601$passed_filters > 500, ]
md.HP2108901 <- md.HP2108901[md.HP2108901$passed_filters > 500, ]
md.HP2110001 <- md.HP2110001[md.HP2110001$passed_filters > 500, ]
md.HP2121601 <- md.HP2121601[md.HP2121601$passed_filters > 500, ]
md.HP2123201 <- md.HP2123201[md.HP2123201$passed_filters > 500, ]
md.HP2132801 <- md.HP2132801[md.HP2132801$passed_filters > 500, ]
md.HP2202101 <- md.HP2202101[md.HP2202101$passed_filters > 500, ]

# create fragment objects
{
  frags.HP2022801 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)",
    cells = rownames(md.HP2022801)
  )
  
  frags.SAMN15877725 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)",
    cells = rownames(md.SAMN15877725)
  )
  
  frags.HP2024001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)",
    cells = rownames(md.HP2024001)
  )
  
  frags.HP2031401 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)",
    cells = rownames(md.HP2031401)
  )
  
  frags.HP2105501 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)",
    cells = rownames(md.HP2105501)
  )
  
  frags.HP2106201 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)",
    cells = rownames(md.HP2106201)
  )
  
  frags.HP2107001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)",
    cells = rownames(md.HP2107001)
  )
  
  frags.HP2107901 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)",
    cells = rownames(md.HP2107901)
  )
  
  frags.HP2108601 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)",
    cells = rownames(md.HP2108601)
  )
  
  frags.HP2108901 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)",
    cells = rownames(md.HP2108901)
  )
  
  frags.HP2110001 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)",
    cells = rownames(md.HP2110001)
  )
  
  frags.HP2121601 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)",
    cells = rownames(md.HP2121601)
  )
  
  frags.HP2123201 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)",
    cells = rownames(md.HP2123201)
  )
  
  frags.HP2132801 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)",
    cells = rownames(md.HP2132801)
  )
  
  frags.HP2202101 <- CreateFragmentObject(
    path = r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)",
    cells = rownames(md.HP2202101)
  )
}

# Quantify Peaks
{
  HP2022801.counts <- FeatureMatrix(
    fragments = frags.HP2022801,
    features = combined.peaks,
    cells = rownames(md.HP2022801)
  )
  
  SAMN15877725.counts <- FeatureMatrix(
    fragments = frags.SAMN15877725,
    features = combined.peaks,
    cells = rownames(md.SAMN15877725)
  )
  
  HP2024001.counts <- FeatureMatrix(
    fragments = frags.HP2024001,
    features = combined.peaks,
    cells = rownames(md.HP2024001)
  )
  
  HP2031401.counts <- FeatureMatrix(
    fragments = frags.HP2031401,
    features = combined.peaks,
    cells = rownames(md.HP2031401)
  )
  
  HP2105501.counts <- FeatureMatrix(
    fragments = frags.HP2105501,
    features = combined.peaks,
    cells = rownames(md.HP2105501)
  )
  
  HP2106201.counts <- FeatureMatrix(
    fragments = frags.HP2106201,
    features = combined.peaks,
    cells = rownames(md.HP2106201)
  )
  
  HP2107001.counts <- FeatureMatrix(
    fragments = frags.HP2107001,
    features = combined.peaks,
    cells = rownames(md.HP2107001)
  )
  
  HP2107901.counts <- FeatureMatrix(
    fragments = frags.HP2107901,
    features = combined.peaks,
    cells = rownames(md.HP2107901)
  )
  
  HP2108601.counts <- FeatureMatrix(
    fragments = frags.HP2108601,
    features = combined.peaks,
    cells = rownames(md.HP2108601)
  )
  
  HP2108901.counts <- FeatureMatrix(
    fragments = frags.HP2108901,
    features = combined.peaks,
    cells = rownames(md.HP2108901)
  )
  
  HP2110001.counts <- FeatureMatrix(
    fragments = frags.HP2110001,
    features = combined.peaks,
    cells = rownames(md.HP2110001)
  )
  
  HP2121601.counts <- FeatureMatrix(
    fragments = frags.HP2121601,
    features = combined.peaks,
    cells = rownames(md.HP2121601)
  )
  
  HP2123201.counts <- FeatureMatrix(
    fragments = frags.HP2123201,
    features = combined.peaks,
    cells = rownames(md.HP2123201)
  )
  
  HP2132801.counts <- FeatureMatrix(
    fragments = frags.HP2132801,
    features = combined.peaks,
    cells = rownames(md.HP2132801)
  )
  
  HP2202101.counts <- FeatureMatrix(
    fragments = frags.HP2202101,
    features = combined.peaks,
    cells = rownames(md.HP2202101)
  )
}

# STEP 2: Create Seurat objects ##
{
  HP2022801_assay <- CreateChromatinAssay(HP2022801.counts, fragments = frags.HP2022801, min.features = 100)
  HP2022801_atac <- CreateSeuratObject(HP2022801_assay, assay = "ATAC", meta.data=md.HP2022801)
  
  SAMN15877725_assay <- CreateChromatinAssay(SAMN15877725.counts, fragments = frags.SAMN15877725, min.features = 100)
  SAMN15877725_atac <- CreateSeuratObject(SAMN15877725_assay, assay = "ATAC", meta.data=md.SAMN15877725)
  
  HP2024001_assay <- CreateChromatinAssay(HP2024001.counts, fragments = frags.HP2024001, min.features = 100)
  HP2024001_atac <- CreateSeuratObject(HP2024001_assay, assay = "ATAC", meta.data=md.HP2024001)
  
  HP2031401_assay <- CreateChromatinAssay(HP2031401.counts, fragments = frags.HP2031401, min.features = 100)
  HP2031401_atac <- CreateSeuratObject(HP2031401_assay, assay = "ATAC", meta.data=md.HP2031401)
  
  HP2105501_assay <- CreateChromatinAssay(HP2105501.counts, fragments = frags.HP2105501, min.features = 100)
  HP2105501_atac <- CreateSeuratObject(HP2105501_assay, assay = "ATAC", meta.data=md.HP2105501)
  
  HP2106201_assay <- CreateChromatinAssay(HP2106201.counts, fragments = frags.HP2106201, min.features = 100)
  HP2106201_atac <- CreateSeuratObject(HP2106201_assay, assay = "ATAC", meta.data=md.HP2106201)
  
  HP2107001_assay <- CreateChromatinAssay(HP2107001.counts, fragments = frags.HP2107001, min.features = 100)
  HP2107001_atac <- CreateSeuratObject(HP2107001_assay, assay = "ATAC", meta.data=md.HP2107001)
  
  HP2107901_assay <- CreateChromatinAssay(HP2107901.counts, fragments = frags.HP2107901, min.features = 100)
  HP2107901_atac <- CreateSeuratObject(HP2107901_assay, assay = "ATAC", meta.data=md.HP2107901)
  
  HP2108601_assay <- CreateChromatinAssay(HP2108601.counts, fragments = frags.HP2108601, min.features = 100)
  HP2108601_atac <- CreateSeuratObject(HP2108601_assay, assay = "ATAC", meta.data=md.HP2108601)
  
  HP2108901_assay <- CreateChromatinAssay(HP2108901.counts, fragments = frags.HP2108901, min.features = 100)
  HP2108901_atac <- CreateSeuratObject(HP2108901_assay, assay = "ATAC", meta.data=md.HP2108901)
  
  HP2110001_assay <- CreateChromatinAssay(HP2110001.counts, fragments = frags.HP2110001, min.features = 100)
  HP2110001_atac <- CreateSeuratObject(HP2110001_assay, assay = "ATAC", meta.data=md.HP2110001)
  
  HP2121601_assay <- CreateChromatinAssay(HP2121601.counts, fragments = frags.HP2121601, min.features = 100)
  HP2121601_atac <- CreateSeuratObject(HP2121601_assay, assay = "ATAC", meta.data=md.HP2121601)
  
  HP2123201_assay <- CreateChromatinAssay(HP2123201.counts, fragments = frags.HP2123201, min.features = 100)
  HP2123201_atac <- CreateSeuratObject(HP2123201_assay, assay = "ATAC", meta.data=md.HP2123201)
  
  HP2132801_assay <- CreateChromatinAssay(HP2132801.counts, fragments = frags.HP2132801, min.features = 100)
  HP2132801_atac <- CreateSeuratObject(HP2132801_assay, assay = "ATAC", meta.data=md.HP2132801)
  
  HP2202101_assay <- CreateChromatinAssay(HP2202101.counts, fragments = frags.HP2202101, min.features = 100)
  HP2202101_atac <- CreateSeuratObject(HP2202101_assay, assay = "ATAC", meta.data=md.HP2202101)
}

# Sample specific Metadata addition
{
  HP2022801_atac$sample <- "HP2022801"
  SAMN15877725_atac$sample <- "SAMN15877725"
  HP2024001_atac$sample <- "HP2024001"
  HP2031401_atac$sample <- "HP2031401"
  HP2105501_atac$sample <- "HP2105501"
  HP2106201_atac$sample <- "HP2106201"
  HP2107001_atac$sample <- "HP2107001"
  HP2107901_atac$sample <- "HP2107901"
  HP2108601_atac$sample <- "HP2108601"
  HP2108901_atac$sample <- "HP2108901"
  HP2110001_atac$sample <- "HP2110001"
  HP2121601_atac$sample <- "HP2121601"
  HP2123201_atac$sample <- "HP2123201"
  HP2132801_atac$sample <- "HP2132801"
  HP2202101_atac$sample <- "HP2202101"
  
  # Sex specific Metadata addition
  HP2022801_atac$sex <- "female"
  SAMN15877725_atac$sex <- "male"
  HP2024001_atac$sex <- "female"
  HP2031401_atac$sex <- "male"
  HP2105501_atac$sex <- "female"
  HP2106201_atac$sex <- "female"
  HP2107001_atac$sex <- "male"
  HP2107901_atac$sex <- "male"
  HP2108601_atac$sex <- "female"
  HP2108901_atac$sex <- "female"
  HP2110001_atac$sex <- "male"
  HP2121601_atac$sex <- "female"
  HP2123201_atac$sex <- "male"
  HP2132801_atac$sex <- "female"
  HP2202101_atac$sex <- "female"
  
  # Ancestry specific Metadata addition
  HP2022801_atac$ancestry <- "white"
  SAMN15877725_atac$ancestry <- "white"
  HP2024001_atac$ancestry <- "white"
  HP2031401_atac$ancestry <- "black"
  HP2105501_atac$ancestry <- "white"
  HP2106201_atac$ancestry <- "black"
  HP2107001_atac$ancestry <- "white"
  HP2107901_atac$ancestry <- "white"
  HP2108601_atac$ancestry <- "white"
  HP2108901_atac$ancestry <- "white"
  HP2110001_atac$ancestry <- "black"
  HP2121601_atac$ancestry <- "black"
  HP2123201_atac$ancestry <- "black"
  HP2132801_atac$ancestry <- "black"
  HP2202101_atac$ancestry <- "black"
            
  # Ancestry and sex specific Metadata addition
  HP2022801_atac$ancestry_sex <- "white_female"
  SAMN15877725_atac$ancestry_sex <- "white_male"
  HP2024001_atac$ancestry_sex <- "white_female"
  HP2031401_atac$ancestry_sex <- "black_male"
  HP2105501_atac$ancestry_sex <- "white_female"
  HP2106201_atac$ancestry_sex <- "black_female"
  HP2107001_atac$ancestry_sex <- "white_male"
  HP2107901_atac$ancestry_sex <- "white_male"
  HP2108601_atac$ancestry_sex <- "white_female"
  HP2108901_atac$ancestry_sex <- "white_female"
  HP2110001_atac$ancestry_sex <- "black_male"
  HP2121601_atac$ancestry_sex <- "black_female"
  HP2123201_atac$ancestry_sex <- "black_male"
  HP2132801_atac$ancestry_sex <- "black_female"
  HP2202101_atac$ancestry_sex <- "black_female"
                            
  # Ancestry, sex and assay specific Metadata addition
  HP2022801_atac$ancestry_sex_atac <- "white_female_atac"
  SAMN15877725_atac$ancestry_sex_atac <- "white_male_atac"
  HP2024001_atac$ancestry_sex_atac <- "white_female_atac"
  HP2031401_atac$ancestry_sex_atac <- "black_male_atac"
  HP2105501_atac$ancestry_sex_atac <- "white_female_atac"
  HP2106201_atac$ancestry_sex_atac <- "black_female_atac"
  HP2107001_atac$ancestry_sex_atac <- "white_male_atac"
  HP2107901_atac$ancestry_sex_atac <- "white_male_atac"
  HP2108601_atac$ancestry_sex_atac <- "white_female_atac"
  HP2108901_atac$ancestry_sex_atac <- "white_female_atac"
  HP2110001_atac$ancestry_sex_atac <- "black_male_atac"
  HP2121601_atac$ancestry_sex_atac <- "black_female_atac"
  HP2123201_atac$ancestry_sex_atac <- "black_male_atac"
  HP2132801_atac$ancestry_sex_atac <- "black_female_atac"
  HP2202101_atac$ancestry_sex_atac <- "black_female_atac"
}

# Add annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# Gen annotation
genome(annotations) <- "hg38"

# add gene information to the object
{
  Annotation(HP2022801_atac) <- annotations
  Annotation(SAMN15877725_atac) <- annotations
  Annotation(HP2024001_atac) <- annotations
  Annotation(HP2031401_atac) <- annotations
  Annotation(HP2105501_atac) <- annotations
  Annotation(HP2106201_atac) <- annotations
  Annotation(HP2107001_atac) <- annotations
  Annotation(HP2107901_atac) <- annotations
  Annotation(HP2108601_atac) <- annotations
  Annotation(HP2108901_atac) <- annotations
  Annotation(HP2110001_atac) <- annotations
  Annotation(HP2121601_atac) <- annotations
  Annotation(HP2123201_atac) <- annotations
  Annotation(HP2132801_atac) <- annotations
  Annotation(HP2202101_atac) <- annotations
  
  # compute nucleosome signal score per cell
  HP2022801_atac <- NucleosomeSignal(object = HP2022801_atac)
  SAMN15877725_atac <- NucleosomeSignal(object = SAMN15877725_atac)
  HP2024001_atac <- NucleosomeSignal(object = HP2024001_atac)
  HP2031401_atac <- NucleosomeSignal(object = HP2031401_atac)
  HP2105501_atac <- NucleosomeSignal(object = HP2105501_atac)
  HP2106201_atac <- NucleosomeSignal(object = HP2106201_atac)
  HP2107001_atac <- NucleosomeSignal(object = HP2107001_atac)
  HP2107901_atac <- NucleosomeSignal(object = HP2107901_atac)
  HP2108601_atac <- NucleosomeSignal(object = HP2108601_atac)
  HP2108901_atac <- NucleosomeSignal(object = HP2108901_atac)
  HP2110001_atac <- NucleosomeSignal(object = HP2110001_atac)
  HP2121601_atac <- NucleosomeSignal(object = HP2121601_atac)
  HP2123201_atac <- NucleosomeSignal(object = HP2123201_atac)
  HP2132801_atac <- NucleosomeSignal(object = HP2132801_atac)
  HP2202101_atac <- NucleosomeSignal(object = HP2202101_atac)
  
  # compute TSS enrichment score per cell
  HP2022801_atac <- TSSEnrichment(object = HP2022801_atac, fast = FALSE)
  SAMN15877725_atac <- TSSEnrichment(object = SAMN15877725_atac, fast = FALSE)
  HP2024001_atac <- TSSEnrichment(object = HP2024001_atac, fast = FALSE)
  HP2031401_atac <- TSSEnrichment(object = HP2031401_atac, fast = FALSE)
  HP2105501_atac <- TSSEnrichment(object = HP2105501_atac, fast = FALSE)
  HP2106201_atac <- TSSEnrichment(object = HP2106201_atac, fast = FALSE)
  HP2107001_atac <- TSSEnrichment(object = HP2107001_atac, fast = FALSE)
  HP2107901_atac <- TSSEnrichment(object = HP2107901_atac, fast = FALSE)
  HP2108601_atac <- TSSEnrichment(object = HP2108601_atac, fast = FALSE)
  HP2108901_atac <- TSSEnrichment(object = HP2108901_atac, fast = FALSE)
  HP2110001_atac <- TSSEnrichment(object = HP2110001_atac, fast = FALSE)
  HP2121601_atac <- TSSEnrichment(object = HP2121601_atac, fast = FALSE)
  HP2123201_atac <- TSSEnrichment(object = HP2123201_atac, fast = FALSE)
  HP2132801_atac <- TSSEnrichment(object = HP2132801_atac, fast = FALSE)
  HP2202101_atac <- TSSEnrichment(object = HP2202101_atac, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  HP2022801_atac$pct_reads_in_peaks <- HP2022801_atac$peak_region_fragments / HP2022801_atac$passed_filters * 100
  HP2022801_atac$blacklist_ratio <- HP2022801_atac$blacklist_region_fragments / HP2022801_atac$peak_region_fragments
  
  SAMN15877725_atac$pct_reads_in_peaks <- SAMN15877725_atac$peak_region_fragments / SAMN15877725_atac$passed_filters * 100
  SAMN15877725_atac$blacklist_ratio <- SAMN15877725_atac$blacklist_region_fragments / SAMN15877725_atac$peak_region_fragments
  
  HP2024001_atac$pct_reads_in_peaks <- HP2024001_atac$peak_region_fragments / HP2024001_atac$passed_filters * 100
  HP2024001_atac$blacklist_ratio <- HP2024001_atac$blacklist_region_fragments / HP2024001_atac$peak_region_fragments
  
  HP2031401_atac$pct_reads_in_peaks <- HP2031401_atac$peak_region_fragments / HP2031401_atac$passed_filters * 100
  HP2031401_atac$blacklist_ratio <- HP2031401_atac$blacklist_region_fragments / HP2031401_atac$peak_region_fragments
  
  HP2105501_atac$pct_reads_in_peaks <- HP2105501_atac$peak_region_fragments / HP2105501_atac$passed_filters * 100
  HP2105501_atac$blacklist_ratio <- HP2105501_atac$blacklist_region_fragments / HP2105501_atac$peak_region_fragments
  
  HP2106201_atac$pct_reads_in_peaks <- HP2106201_atac$peak_region_fragments / HP2106201_atac$passed_filters * 100
  HP2106201_atac$blacklist_ratio <- HP2106201_atac$blacklist_region_fragments / HP2106201_atac$peak_region_fragments
  
  HP2107001_atac$pct_reads_in_peaks <- HP2107001_atac$peak_region_fragments / HP2107001_atac$passed_filters * 100
  HP2107001_atac$blacklist_ratio <- HP2107001_atac$blacklist_region_fragments / HP2107001_atac$peak_region_fragments
  
  HP2107901_atac$pct_reads_in_peaks <- HP2107901_atac$peak_region_fragments / HP2107901_atac$passed_filters * 100
  HP2107901_atac$blacklist_ratio <- HP2107901_atac$blacklist_region_fragments / HP2107901_atac$peak_region_fragments
  
  HP2108601_atac$pct_reads_in_peaks <- HP2108601_atac$peak_region_fragments / HP2108601_atac$passed_filters * 100
  HP2108601_atac$blacklist_ratio <- HP2108601_atac$blacklist_region_fragments / HP2108601_atac$peak_region_fragments
  
  HP2108901_atac$pct_reads_in_peaks <- HP2108901_atac$peak_region_fragments / HP2108901_atac$passed_filters * 100
  HP2108901_atac$blacklist_ratio <- HP2108901_atac$blacklist_region_fragments / HP2108901_atac$peak_region_fragments
  
  HP2110001_atac$pct_reads_in_peaks <- HP2110001_atac$peak_region_fragments / HP2110001_atac$passed_filters * 100
  HP2110001_atac$blacklist_ratio <- HP2110001_atac$blacklist_region_fragments / HP2110001_atac$peak_region_fragments
  
  HP2121601_atac$pct_reads_in_peaks <- HP2121601_atac$peak_region_fragments / HP2121601_atac$passed_filters * 100
  HP2121601_atac$blacklist_ratio <- HP2121601_atac$blacklist_region_fragments / HP2121601_atac$peak_region_fragments
  
  HP2123201_atac$pct_reads_in_peaks <- HP2123201_atac$peak_region_fragments / HP2123201_atac$passed_filters * 100
  HP2123201_atac$blacklist_ratio <- HP2123201_atac$blacklist_region_fragments / HP2123201_atac$peak_region_fragments
  
  HP2132801_atac$pct_reads_in_peaks <- HP2132801_atac$peak_region_fragments / HP2132801_atac$passed_filters * 100
  HP2132801_atac$blacklist_ratio <- HP2132801_atac$blacklist_region_fragments / HP2132801_atac$peak_region_fragments
  
  HP2202101_atac$pct_reads_in_peaks <- HP2202101_atac$peak_region_fragments / HP2202101_atac$passed_filters * 100
  HP2202101_atac$blacklist_ratio <- HP2202101_atac$blacklist_region_fragments / HP2202101_atac$peak_region_fragments
}

# view QC
HP2022801_atac$high.tss <- ifelse(HP2022801_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2022801_atac, group.by = 'high.tss') + NoLegend()

SAMN15877725_atac$high.tss <- ifelse(SAMN15877725_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(SAMN15877725_atac, group.by = 'high.tss') + NoLegend()

HP2024001_atac$high.tss <- ifelse(HP2024001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2024001_atac, group.by = 'high.tss') + NoLegend()

HP2031401_atac$high.tss <- ifelse(HP2031401_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2031401_atac, group.by = 'high.tss') + NoLegend()

HP2105501_atac$high.tss <- ifelse(HP2105501_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2105501_atac, group.by = 'high.tss') + NoLegend()

HP2106201_atac$high.tss <- ifelse(HP2106201_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2106201_atac, group.by = 'high.tss') + NoLegend()

HP2107001_atac$high.tss <- ifelse(HP2107001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2107001_atac, group.by = 'high.tss') + NoLegend()

HP2107901_atac$high.tss <- ifelse(HP2107901_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2107901_atac, group.by = 'high.tss') + NoLegend()

HP2108601_atac$high.tss <- ifelse(HP2108601_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2108601_atac, group.by = 'high.tss') + NoLegend()

HP2108901_atac$high.tss <- ifelse(HP2108901_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2108901_atac, group.by = 'high.tss') + NoLegend()

HP2110001_atac$high.tss <- ifelse(HP2110001_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2110001_atac, group.by = 'high.tss') + NoLegend()

HP2121601_atac$high.tss <- ifelse(HP2121601_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2121601_atac, group.by = 'high.tss') + NoLegend()

HP2123201_atac$high.tss <- ifelse(HP2123201_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2123201_atac, group.by = 'high.tss') + NoLegend()

HP2132801_atac$high.tss <- ifelse(HP2132801_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2132801_atac, group.by = 'high.tss') + NoLegend()

HP2202101_atac$high.tss <- ifelse(HP2202101_atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(HP2202101_atac, group.by = 'high.tss') + NoLegend()

# View Nucleosome signal  
HP2022801_atac$nucleosome_group <- ifelse(HP2022801_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2022801_atac, group.by = 'nucleosome_group')

SAMN15877725_atac$nucleosome_group <- ifelse(SAMN15877725_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = SAMN15877725_atac, group.by = 'nucleosome_group')

HP2024001_atac$nucleosome_group <- ifelse(HP2024001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2024001_atac, group.by = 'nucleosome_group')

HP2031401_atac$nucleosome_group <- ifelse(HP2031401_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2031401_atac, group.by = 'nucleosome_group')

HP2105501_atac$nucleosome_group <- ifelse(HP2105501_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2105501_atac, group.by = 'nucleosome_group')

HP2106201_atac$nucleosome_group <- ifelse(HP2106201_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2106201_atac, group.by = 'nucleosome_group')

HP2107001_atac$nucleosome_group <- ifelse(HP2107001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2107001_atac, group.by = 'nucleosome_group')

HP2107901_atac$nucleosome_group <- ifelse(HP2107901_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2107901_atac, group.by = 'nucleosome_group')

HP2108601_atac$nucleosome_group <- ifelse(HP2108601_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2108601_atac, group.by = 'nucleosome_group')

HP2108901_atac$nucleosome_group <- ifelse(HP2108901_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2108901_atac, group.by = 'nucleosome_group')

HP2110001_atac$nucleosome_group <- ifelse(HP2110001_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2110001_atac, group.by = 'nucleosome_group')

HP2121601_atac$nucleosome_group <- ifelse(HP2121601_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2121601_atac, group.by = 'nucleosome_group')

HP2123201_atac$nucleosome_group <- ifelse(HP2123201_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2123201_atac, group.by = 'nucleosome_group')

HP2132801_atac$nucleosome_group <- ifelse(HP2132801_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2132801_atac, group.by = 'nucleosome_group')

HP2202101_atac$nucleosome_group <- ifelse(HP2202101_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HP2202101_atac, group.by = 'nucleosome_group')

# Count fragments
HP2022801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\1_220628_Fahd_snATAC1_HP-20228-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
SAMN15877725_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\2_220701_Fahd_snATAC2_SAMN15877725\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2024001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\3_220701_Fahd_snATAC3_HP-20240-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2031401_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\4_220630_Fahd_snATAC4_HP-20314-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2105501_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\5_220303_snATAC_F52_HP-21055-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2106201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\6_210401_snATAC_F62_HP-21062-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\7_210401_snATAC_F7a_HP-21070-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2107901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\9_210628_snATAC_F9a_HP-21079-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\10_210628_snATAC_F10a_HP-21086-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2108901_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\11_210714_snATAC_F11a_HP-21089-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2110001_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\12_210714_snATAC_F12a_HP-21100-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2121601_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\13_211208_snATAC_F13_HP-21216-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2123201_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\14_211208_snATAC_F14_HP-21232-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2132801_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\15_220303_snATAC_F15a_HP-21328-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)
HP2202101_atac_fraginfo <- CountFragments(r"(E:\1.SexbasedStudyrawdata\Cellranger_raw_data\snATACseq\16_220630_Fahd_snATAC16_HP-22021-01\fragments.tsv.gz)", cells = NULL, max_lines = NULL, verbose = TRUE)

rownames(HP2022801_atac_fraginfo) <- HP2022801_atac_fraginfo$CB
rownames(SAMN15877725_atac_fraginfo) <- SAMN15877725_atac_fraginfo$CB
rownames(HP2024001_atac_fraginfo) <- HP2024001_atac_fraginfo$CB
rownames(HP2031401_atac_fraginfo) <- HP2031401_atac_fraginfo$CB
rownames(HP2105501_atac_fraginfo) <- HP2105501_atac_fraginfo$CB
rownames(HP2106201_atac_fraginfo) <- HP2106201_atac_fraginfo$CB
rownames(HP2107001_atac_fraginfo) <- HP2107001_atac_fraginfo$CB
rownames(HP2107901_atac_fraginfo) <- HP2107901_atac_fraginfo$CB
rownames(HP2108601_atac_fraginfo) <- HP2108601_atac_fraginfo$CB
rownames(HP2108901_atac_fraginfo) <- HP2108901_atac_fraginfo$CB
rownames(HP2110001_atac_fraginfo) <- HP2110001_atac_fraginfo$CB
rownames(HP2121601_atac_fraginfo) <- HP2121601_atac_fraginfo$CB
rownames(HP2123201_atac_fraginfo) <- HP2123201_atac_fraginfo$CB
rownames(HP2132801_atac_fraginfo) <- HP2132801_atac_fraginfo$CB
rownames(HP2202101_atac_fraginfo) <- HP2202101_atac_fraginfo$CB

HP2022801_atac$fragments <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "frequency_count"]
HP2022801_atac$mononucleosomal <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "mononucleosomal"]
HP2022801_atac$nucleosome_free <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "nucleosome_free"]
HP2022801_atac$reads_count <- HP2022801_atac_fraginfo[colnames(HP2022801_atac), "reads_count"]

SAMN15877725_atac$fragments <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "frequency_count"]
SAMN15877725_atac$mononucleosomal <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "mononucleosomal"]
SAMN15877725_atac$nucleosome_free <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "nucleosome_free"]
SAMN15877725_atac$reads_count <- SAMN15877725_atac_fraginfo[colnames(SAMN15877725_atac), "reads_count"]

HP2024001_atac$fragments <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "frequency_count"]
HP2024001_atac$mononucleosomal <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "mononucleosomal"]
HP2024001_atac$nucleosome_free <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "nucleosome_free"]
HP2024001_atac$reads_count <- HP2024001_atac_fraginfo[colnames(HP2024001_atac), "reads_count"]

HP2031401_atac$fragments <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "frequency_count"]
HP2031401_atac$mononucleosomal <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "mononucleosomal"]
HP2031401_atac$nucleosome_free <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "nucleosome_free"]
HP2031401_atac$reads_count <- HP2031401_atac_fraginfo[colnames(HP2031401_atac), "reads_count"]

HP2105501_atac$fragments <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "frequency_count"]
HP2105501_atac$mononucleosomal <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "mononucleosomal"]
HP2105501_atac$nucleosome_free <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "nucleosome_free"]
HP2105501_atac$reads_count <- HP2105501_atac_fraginfo[colnames(HP2105501_atac), "reads_count"]

HP2106201_atac$fragments <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "frequency_count"]
HP2106201_atac$mononucleosomal <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "mononucleosomal"]
HP2106201_atac$nucleosome_free <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "nucleosome_free"]
HP2106201_atac$reads_count <- HP2106201_atac_fraginfo[colnames(HP2106201_atac), "reads_count"]

HP2107001_atac$fragments <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "frequency_count"]
HP2107001_atac$mononucleosomal <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "mononucleosomal"]
HP2107001_atac$nucleosome_free <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "nucleosome_free"]
HP2107001_atac$reads_count <- HP2107001_atac_fraginfo[colnames(HP2107001_atac), "reads_count"]

HP2107901_atac$fragments <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "frequency_count"]
HP2107901_atac$mononucleosomal <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "mononucleosomal"]
HP2107901_atac$nucleosome_free <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "nucleosome_free"]
HP2107901_atac$reads_count <- HP2107901_atac_fraginfo[colnames(HP2107901_atac), "reads_count"]

HP2108601_atac$fragments <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "frequency_count"]
HP2108601_atac$mononucleosomal <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "mononucleosomal"]
HP2108601_atac$nucleosome_free <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "nucleosome_free"]
HP2108601_atac$reads_count <- HP2108601_atac_fraginfo[colnames(HP2108601_atac), "reads_count"]

HP2108901_atac$fragments <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "frequency_count"]
HP2108901_atac$mononucleosomal <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "mononucleosomal"]
HP2108901_atac$nucleosome_free <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "nucleosome_free"]
HP2108901_atac$reads_count <- HP2108901_atac_fraginfo[colnames(HP2108901_atac), "reads_count"]

HP2110001_atac$fragments <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "frequency_count"]
HP2110001_atac$mononucleosomal <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "mononucleosomal"]
HP2110001_atac$nucleosome_free <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "nucleosome_free"]
HP2110001_atac$reads_count <- HP2110001_atac_fraginfo[colnames(HP2110001_atac), "reads_count"]

HP2121601_atac$fragments <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "frequency_count"]
HP2121601_atac$mononucleosomal <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "mononucleosomal"]
HP2121601_atac$nucleosome_free <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "nucleosome_free"]
HP2121601_atac$reads_count <- HP2121601_atac_fraginfo[colnames(HP2121601_atac), "reads_count"]

HP2123201_atac$fragments <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "frequency_count"]
HP2123201_atac$mononucleosomal <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "mononucleosomal"]
HP2123201_atac$nucleosome_free <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "nucleosome_free"]
HP2123201_atac$reads_count <- HP2123201_atac_fraginfo[colnames(HP2123201_atac), "reads_count"]

HP2132801_atac$fragments <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "frequency_count"]
HP2132801_atac$mononucleosomal <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "mononucleosomal"]
HP2132801_atac$nucleosome_free <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "nucleosome_free"]
HP2132801_atac$reads_count <- HP2132801_atac_fraginfo[colnames(HP2132801_atac), "reads_count"]

HP2202101_atac$fragments <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "frequency_count"]
HP2202101_atac$mononucleosomal <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "mononucleosomal"]
HP2202101_atac$nucleosome_free <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "nucleosome_free"]
HP2202101_atac$reads_count <- HP2202101_atac_fraginfo[colnames(HP2202101_atac), "reads_count"]


HP2022801_atac <- FRiP(
  object = HP2022801_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

SAMN15877725_atac <- FRiP(
  object = SAMN15877725_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2024001_atac <- FRiP(
  object = HP2024001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2031401_atac <- FRiP(
  object = HP2031401_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2105501_atac <- FRiP(
  object = HP2105501_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2106201_atac <- FRiP(
  object = HP2106201_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2107001_atac <- FRiP(
  object = HP2107001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2107901_atac <- FRiP(
  object = HP2107901_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2108601_atac <- FRiP(
  object = HP2108601_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2108901_atac <- FRiP(
  object = HP2108901_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2110001_atac <- FRiP(
  object = HP2110001_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2121601_atac <- FRiP(
  object = HP2121601_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2123201_atac <- FRiP(
  object = HP2123201_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2132801_atac <- FRiP(
  object = HP2132801_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

HP2202101_atac <- FRiP(
  object = HP2202101_atac,
  assay = 'ATAC',
  total.fragments = 'fragments')

# Add unique cell names otherwise integration will give errors
{
  HP2022801_atac <- RenameCells(object = HP2022801_atac, add.cell.id = "HP2022801")
  SAMN15877725_atac <- RenameCells(object = SAMN15877725_atac, add.cell.id = "SAMN15877725")
  HP2024001_atac <- RenameCells(object = HP2024001_atac, add.cell.id = "HP2024001")
  HP2031401_atac <- RenameCells(object = HP2031401_atac, add.cell.id = "HP2031401")
  HP2105501_atac <- RenameCells(object = HP2105501_atac, add.cell.id = "HP2105501")
  HP2106201_atac <- RenameCells(object = HP2106201_atac, add.cell.id = "HP2106201")
  HP2107001_atac <- RenameCells(object = HP2107001_atac, add.cell.id = "HP2107001")
  HP2107901_atac <- RenameCells(object = HP2107901_atac, add.cell.id = "HP2107901")
  HP2108601_atac <- RenameCells(object = HP2108601_atac, add.cell.id = "HP2108601")
  HP2108901_atac <- RenameCells(object = HP2108901_atac, add.cell.id = "HP2108901")
  HP2110001_atac <- RenameCells(object = HP2110001_atac, add.cell.id = "HP2110001")
  HP2121601_atac <- RenameCells(object = HP2121601_atac, add.cell.id = "HP2121601")
  HP2123201_atac <- RenameCells(object = HP2123201_atac, add.cell.id = "HP2123201")
  HP2132801_atac <- RenameCells(object = HP2132801_atac, add.cell.id = "HP2132801")
  HP2202101_atac <- RenameCells(object = HP2202101_atac, add.cell.id = "HP2202101")
}

# Visualize QC
VlnPlot(
  object = HP2022801_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = SAMN15877725_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2024001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2031401_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2105501_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2106201_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2107001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2107901_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2108601_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2108901_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2110001_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2121601_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2123201_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2132801_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

VlnPlot(
  object = HP2202101_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6
)

# Add Gene activity matrix
{
  gene.activities.HP2022801 <- GeneActivity(HP2022801_atac)
  gene.activities.SAMN15877725 <- GeneActivity(SAMN15877725_atac)
  gene.activities.HP2024001 <- GeneActivity(HP2024001_atac)
  gene.activities.HP2031401 <- GeneActivity(HP2031401_atac)
  gene.activities.HP2105501 <- GeneActivity(HP2105501_atac)
  gene.activities.HP2106201 <- GeneActivity(HP2106201_atac)
  gene.activities.HP2107001 <- GeneActivity(HP2107001_atac)
  gene.activities.HP2107901 <- GeneActivity(HP2107901_atac)
  gene.activities.HP2108601 <- GeneActivity(HP2108601_atac)
  gene.activities.HP2108901 <- GeneActivity(HP2108901_atac)
  gene.activities.HP2110001 <- GeneActivity(HP2110001_atac)
  gene.activities.HP2121601 <- GeneActivity(HP2121601_atac)
  gene.activities.HP2123201 <- GeneActivity(HP2123201_atac)
  gene.activities.HP2132801 <- GeneActivity(HP2132801_atac)
  gene.activities.HP2202101 <- GeneActivity(HP2202101_atac)
}

# add the gene activity matrix to the Seurat object as a new assay and normalize it
{
  HP2022801_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2022801)
  HP2022801_atac <- NormalizeData(
    object = HP2022801_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2022801_atac$nCount_RNA)
  )
  
  SAMN15877725_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.SAMN15877725)
  SAMN15877725_atac <- NormalizeData(
    object = SAMN15877725_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(SAMN15877725_atac$nCount_RNA)
  )
  
  HP2024001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2024001)
  HP2024001_atac <- NormalizeData(
    object = HP2024001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2024001_atac$nCount_RNA)
  )
  
  HP2031401_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2031401)
  HP2031401_atac <- NormalizeData(
    object = HP2031401_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2031401_atac$nCount_RNA)
  )
  
  HP2105501_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2105501)
  HP2105501_atac <- NormalizeData(
    object = HP2105501_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2105501_atac$nCount_RNA)
  )
  
  HP2106201_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2106201)
  HP2106201_atac <- NormalizeData(
    object = HP2106201_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2106201_atac$nCount_RNA)
  )
  
  HP2107001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2107001)
  HP2107001_atac <- NormalizeData(
    object = HP2107001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2107001_atac$nCount_RNA)
  )
  
  HP2107901_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2107901)
  HP2107901_atac <- NormalizeData(
    object = HP2107901_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2107901_atac$nCount_RNA)
  )
  
  HP2108601_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2108601)
  HP2108601_atac <- NormalizeData(
    object = HP2108601_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2108601_atac$nCount_RNA)
  )
  
  HP2108901_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2108901)
  HP2108901_atac <- NormalizeData(
    object = HP2108901_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2108901_atac$nCount_RNA)
  )
  
  HP2110001_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2110001)
  HP2110001_atac <- NormalizeData(
    object = HP2110001_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2110001_atac$nCount_RNA)
  )
  
  HP2121601_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2121601)
  HP2121601_atac <- NormalizeData(
    object = HP2121601_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2121601_atac$nCount_RNA)
  )
  
  HP2123201_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2123201)
  HP2123201_atac <- NormalizeData(
    object = HP2123201_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2123201_atac$nCount_RNA)
  )
  
  HP2132801_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2132801)
  HP2132801_atac <- NormalizeData(
    object = HP2132801_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2132801_atac$nCount_RNA)
  )
  
  HP2202101_atac[['RNA']] <- CreateAssayObject(counts = gene.activities.HP2202101)
  HP2202101_atac <- NormalizeData(
    object = HP2202101_atac,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(HP2202101_atac$nCount_RNA)
  )
}


# merge all datasets, adding a cell ID to make sure cell names are unique
combined_atac <- merge(
  x = HP2022801_atac,
  y = list(SAMN15877725_atac, HP2024001_atac, HP2031401_atac, HP2105501_atac,
           HP2106201_atac, HP2107001_atac, HP2107901_atac, HP2108601_atac, 
           HP2108901_atac, HP2110001_atac, HP2121601_atac, HP2123201_atac,
           HP2132801_atac, HP2202101_atac))

VlnPlot(
  object = combined_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 6, raster = FALSE
)

# QC Cleanup
combined_atac <- subset(
  x = combined_atac, 
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 30 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    FRiP > 0.2
)
combined_atac

#Save file
#qsave(combined_atac, file = r"(E:\2.SexbasedStudyCurrent\QS files\combined_atac.qs)")
#qsave(combined_atac_doublet, file = r"(E:\2.SexbasedStudyCurrent\QS files\combined_atac_doublet.qs)")
combined_atac <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac.qs)")
combined_atac_doublet <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac_doublet.qs)")

############################ STAGE ############################
############################   3   ############################
# NORMALIZATION AND BATCH CORRECTION
# Run TFDIF  
combined_atac <- FindTopFeatures(combined_atac, min.cutoff = 20)
combined_atac <- RunTFIDF(combined_atac)
combined_atac <- RunSVD(combined_atac)
combined_atac <- RunUMAP(combined_atac, dims = 2:30, reduction = 'lsi')
DimPlot(combined_atac, group.by = 'ancestry_sex', pt.size = 0.1)

# Batch correction using Harmony
hm.integrated <- RunHarmony(object = combined_atac, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated, group.by = 'ancestry_sex', pt.size = 0.1)

# We need to subset out all of those cell IDs in our doublet calculated data which are present in actual ATAC data
#create cell names as metadata colum
atac_cellnames <- colnames(combined_atac)

#subset the doublet file to have only those barcodes that have been identified as high quality cells
combined_atac_doublet <- subset(combined_atac_doublet, rownames(combined_atac_doublet) %in% atac_cellnames)

# Adding a column to define multiplet vs singlet
combined_atac_doublet$doublets <- ifelse(combined_atac_doublet$p.value >= 0.05, "singlet", "multiplet" )

# Check number of singlets and multiplets
sum(combined_atac_doublet$doublets == "singlet")
sum(combined_atac_doublet$doublets == "multiplet")

#Some more checking because AddMetadata was driving me INSANE
#cells_seuratobj <- Cells(combined_atac)
#cells_netadatobj <- rownames(combined_atac_doublet)
#setequal(cells_seuratobj, cells_netadatobj)

# subset out the columns you wish to add only, here we are only adding the doublets column
doubletdat <- subset(combined_atac_doublet, select = c("doublets"))

# Add doublet metadata to seurat object
hm.integrated <- AddMetaData(object = hm.integrated, metadata = doubletdat, col.name = 'doublets')
table(hm.integrated$doublets) #party time

# Check doublets
#Idents(hm.integrated) <- "doublets"
DimPlot(hm.integrated, group.by = 'doublets', pt.size = 0.1)

# Remove doublets
Idents(hm.integrated) <- "doublets" 
hm.integrated.dfree <- subset(hm.integrated, idents = c("singlet"))
hm.integrated.dfree

# Re-Run TFDIF  
hm.integrated.dfree <- FindTopFeatures(hm.integrated.dfree, min.cutoff = 20)
hm.integrated.dfree <- RunTFIDF(hm.integrated.dfree)
hm.integrated.dfree <- RunSVD(hm.integrated.dfree)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'lsi')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Batch correction using Harmony
hm.integrated.dfree <- RunHarmony(object = hm.integrated.dfree, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Normalize gene activities
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- NormalizeData(hm.integrated.dfree)
hm.integrated.dfree <- ScaleData(hm.integrated.dfree, features = rownames(hm.integrated.dfree))

############################ STAGE ############################
############################   4   ############################
# SUB_CLUSTERING QC AND RE-NORMALIZATION
# Clustering
DefaultAssay(hm.integrated.dfree) <- "ATAC"
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'harmony', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()

#Subset out poor quality cells
Idents(hm.integrated.dfree) <- "seurat_clusters"
hm.integrated.dfree <- subset(hm.integrated.dfree, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                                              #"11", 
                                                              "12", "13", #"14", "15", 
                                                              "16", #"17", 
                                                              "18", "19", "20", "21"))

# Re-Run TFDIF  
hm.integrated.dfree <- FindTopFeatures(hm.integrated.dfree, min.cutoff = 20)
hm.integrated.dfree <- RunTFIDF(hm.integrated.dfree)
hm.integrated.dfree <- RunSVD(hm.integrated.dfree)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'lsi')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Normalize gene activities
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- NormalizeData(hm.integrated.dfree)
hm.integrated.dfree <- ScaleData(hm.integrated.dfree, features = rownames(hm.integrated.dfree))

# Batch correction using Harmony
hm.integrated.dfree <- RunHarmony(object = hm.integrated.dfree, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated.dfree <- RunUMAP(hm.integrated.dfree, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated.dfree, group.by = 'ancestry_sex', pt.size = 0.1)

# Clusterning UMAP was created on basis of ATAC profile
DefaultAssay(hm.integrated.dfree) <- "ATAC"
hm.integrated.dfree <- FindNeighbors(object = hm.integrated.dfree, reduction = 'harmony', dims = 2:30)
hm.integrated.dfree <- FindClusters(object = hm.integrated.dfree, verbose = FALSE, algorithm = 3)
DimPlot(object = hm.integrated.dfree, label = TRUE) #+ NoLegend()

#Save file
#qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
#combined_atac_doublet <- readRDS(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\snATACseq\combined_atac_doublet.rds)")

# Open necessary scRNAseq and snATAC data
#hm.integrated.dfree <- readRDS(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
#pancreas.combined.h.s <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

############################ STAGE ############################
############################   5   ############################
# TRANSFERENCE OF ANCHORS and celltype annotation
# Identify anchors
# Read in RNA data
DefaultAssay(pancreas.combined.h.s) <- 'RNA'
system.time({
#  user      system  elapsed 
#  25456.72  2710.84 28241.15 ~7.8hrs
transfer.anchors <- FindTransferAnchors(reference = pancreas.combined.h.s, 
                                        query = hm.integrated.dfree, 
                                        features = pancreas.combined.h.s@assays$RNA@var.features,
                                        reference.assay = "RNA", 
                                        query.assay = "RNA", 
                                        reduction = "cca")
})

# Save
#qsave(transfer.anchors, file = r"(E:\2.SexbasedStudyCurrent\QS files\transfer.anchors.qs)")

#### START VERY EARLY FROM HERE ##
####                            ##
transfer.anchors <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\transfer.anchors.qs)")
pancreas.combined.h.s <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# Annotation of scATAC cells via label transfer  
# map query onto the reference dataset
DefaultAssay(pancreas.combined.h.s) <- "RNA"
# you dont need to run this because return.model = T when generating UMAP for scRNAseq
#DimPlot(pancreas.combined.h.s, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
#pancreas.combined.h.s <- RunUMAP(object = pancreas.combined.h.s, assay = "SCT", reduction = "harmony", dims = 1:30, return.model = TRUE) # return model = TRUE

# View clustering
DefaultAssay(hm.integrated.dfree) <- "ATAC"
DimPlot(object = hm.integrated.dfree, label = TRUE) + NoLegend()

# Query mapping
hm.integrated.dfree <- MapQuery(
  anchorset = transfer.anchors,
  reference = pancreas.combined.h.s,
  query = hm.integrated.dfree,
  refdata = pancreas.combined.h.s$celltype_qadir,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

# Clusterning UMAP was created on basis of ATAC profile
# Add predicted ID data to the new Signac seurat object
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = pancreas.combined.h.s$celltype_qadir,
                                     weight.reduction = hm.integrated.dfree[["harmony"]], 
                                     dims = 2:30)

hm.integrated.dfree <- AddMetaData(hm.integrated.dfree, metadata = celltype.predictions)

# View UMAP
DimPlot(hm.integrated.dfree, group.by = "predicted.id", reduction = "umap", label = TRUE) + ggtitle("Predicted annotation")
DimPlot(hm.integrated.dfree, group.by = "ATAC_snn_res.0.8", reduction = "umap", label = TRUE) + ggtitle("Res")# + nolegend()
DimPlot(pancreas.combined.h.s, group.by = "celltype_qadir", reduction = "umap", label = TRUE) + ggtitle("Celltype Classification")

# rename Idents ans save as celltype
Idents(hm.integrated.dfree) <- "ATAC_snn_res.0.8"
hm.integrated.dfree <- RenameIdents(hm.integrated.dfree,
                                    "0" = "alpha", 
                                    "1" = "beta",
                                    "2" = "alpha", 
                                    "3" = "alpha",
                                    "4" = "beta", 
                                    "5" = "activated_stellate", #(PDGFRA+)
                                    "6" = "ductal", 
                                    "7" = "delta",
                                    "8" = "acinar", 
                                    "9" = "alpha",
                                    "10" = "beta", 
                                    "11" = "gamma",
                                    "12" = "endothelial",
                                    "13" = "endothelial", 
                                    "14" = "quiescent_stellate", #(RGS5+)
                                    "15" = "macrophage",
                                    "16" = "lymphocyte"
)
table(Idents(hm.integrated.dfree))
DimPlot(hm.integrated.dfree, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(hm.integrated.dfree))
hm.integrated.dfree$celltype <- Idents(hm.integrated.dfree)
head(hm.integrated.dfree@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "alpha", "delta", "gamma",
               "ductal", "acinar", 
               "activated_stellate", "quiescent_stellate", "endothelial",
               "macrophage", "lymphocyte"
)
# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree@meta.data$celltype <- factor(x = hm.integrated.dfree@meta.data$celltype, levels = my_levels)
Idents(hm.integrated.dfree) <- "celltype"
table(Idents(hm.integrated.dfree))

# Observing cells
Idents(hm.integrated.dfree) <- "celltype"
DimPlot(hm.integrated.dfree, 
        #split.by = "ancestry_sex", 
        #group.by = "celltype", 
        label = FALSE, 
        ncol = 2,  
        cols = c("dodgerblue3", #"beta"
                 "lightseagreen", #"alpha"
                 "chartreuse3", #"delta"
                 "springgreen4", #"gamma"
                 "darkorange2", #"ductal"
                 "salmon3", #"acinar"
                 "orange", #"activated-stellate"
                 "salmon", #"quiescent-stellate"
                 "red", #"endothelial"
                 "magenta3", #"macrophages"
                 "orchid1" #"lymphocyte"
                 )
) + ggtitle("Celltype Classification")


# Plotting
DefaultAssay(hm.integrated.dfree) <- "RNA"
FeaturePlot(
  object = hm.integrated.dfree,
  features = c('SST', 'INS'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  reduction = "umap",
  ncol = 3
)

# Add metadata for DT
Idents(object = hm.integrated.dfree) <- "celltype"
hm.integrated.dfree$celltype_sex <- paste(Idents(hm.integrated.dfree), hm.integrated.dfree$sex, sep = "_")
hm.integrated.dfree$celltype_sex_ancestry <- paste(Idents(hm.integrated.dfree), hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, sep = "_")
table(hm.integrated.dfree@meta.data$celltype_sex)
table(hm.integrated.dfree@meta.data$celltype_sex_ancestry)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("beta_female_black", "beta_male_black", "beta_female_white", "beta_male_white", 
                "alpha_female_black", "alpha_male_black", "alpha_female_white", "alpha_male_white", 
                "delta_female_black", "delta_male_black", "delta_female_white", "delta_male_white",
                "gamma_female_black", "gamma_male_black", "gamma_female_white", "gamma_male_white",
                "ductal_female_black", "ductal_male_black", "ductal_female_white", "ductal_male_white",
                "acinar_female_black", "acinar_male_black", "acinar_female_white", "acinar_male_white",
                "activated_stellate_female_black", "activated_stellate_male_black", "activated_stellate_female_white", "activated_stellate_male_white", 
                "quiescent_stellate_female_black", "quiescent_stellate_male_black", "quiescent_stellate_female_white", "quiescent_stellate_male_white", 
                "endothelial_female_black", "endothelial_male_black", "endothelial_female_white", "endothelial_male_white", 
                "macrophage_female_black", "macrophage_male_black", "macrophage_female_white", "macrophage_male_white",
                "lymphocyte_female_black", "lymphocyte_male_black", "lymphocyte_female_white", "lymphocyte_male_white"
)

my_levels3 <- c("beta_female", "beta_male", 
                "alpha_female", "alpha_male", 
                "delta_female", "delta_male",
                "gamma_female", "gamma_male",
                "ductal_female", "ductal_male",
                "acinar_female", "acinar_male",
                "activated_stellate_female", "activated_stellate_male", 
                "quiescent_stellate_female", "quiescent_stellate_male", 
                "endothelial_female", "endothelial_male", 
                "macrophage_female", "macrophage_male",
                "lymphocyte_female", "lymphocyte_male"
)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree@meta.data$celltype_sex_ancestry <- factor(x = hm.integrated.dfree@meta.data$celltype_sex_ancestry, levels = my_levels2)
hm.integrated.dfree@meta.data$celltype_sex <- factor(x = hm.integrated.dfree@meta.data$celltype_sex, levels = my_levels3)
table(hm.integrated.dfree@meta.data$celltype_sex_ancestry)
table(hm.integrated.dfree@meta.data$celltype_sex)

#Peak calling
CoveragePlot(
  object = hm.integrated.dfree,
  group.by = "celltype",
  region = c("INS", 
             "GCG",
             "SST",
             "PPY"
             ),
  show.bulk = TRUE,
  #ymax = 80,
  ranges.title = "Ranges",
  window = 200,
  extend.upstream = 25000,
  extend.downstream = 10000,
) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) & 
  scale_fill_manual(values = c("dodgerblue3", #"beta"
                               "lightseagreen", #"alpha"
                               "chartreuse3", #"delta"
                               "springgreen4", #"gamma"
                               "darkorange2", #"ductal"
                               "salmon3", #"acinar"
                               "orange", #"activated-stellate"
                               "salmon", #"quiescent-stellate"
                               "red", #"endothelial"
                               "magenta3", #"macrophages"
                               "orchid1", #"lymphocyte"
                               "grey"
)) 

# Bask in the beauty of your optimised seurat object
hm.integrated.dfree[['ATAC']]
granges(hm.integrated.dfree)
table(hm.integrated.dfree@assays[["ATAC"]]@annotation@seqinfo@genome)

# Change back to peak data
DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype"

# Adding predicted RNA gene expression values
rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(pancreas.combined.h.s, assay = "RNA", slot = "data"),
  weight.reduction = hm.integrated.dfree[["lsi"]],
  dims = 2:30
)

# add predicted values as a new assay
hm.integrated.dfree[["predicted"]] <- rna

############################ STAGE ############################
############################   6   ############################
# MOTIF INFORMATION ADDITION
# Adding Motifs
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
hm.integrated.dfree <- AddMotifs(
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#Compute motif activity
# THIS STEP IS CRAZYYYYYYYYYYYY IT TAKES AGES OR BSOD's PC IF PARALLEL IS NOT TURNED OFF
register(SerialParam()) # VERY IMPORTANT TRANSITION FROM PARALLEL TO SERIAL COMPUTING
hm.integrated.dfree <- RunChromVAR(
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(hm.integrated.dfree) <- 'chromvar'

# SAVE YOUR FILE NOW
qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# look at the activity of Mef2c
FeaturePlot(
  object = hm.integrated.dfree,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

############################ STAGE ############################
############################   10  ############################
# LINKAGE ANALYSIS ##
#Load file
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
Idents(hm.integrated.dfree) <- "celltype"

# Signac issue #496
if (is.null(seqinfo(hm.integrated.dfree))) {
  message("OOPSssssSSSSssss SOMEONE FORGOT TO RUN ....genome = 'hg38'... WHEN THEY RAN creatChromatinAssay() HAHAHAHAHAHAHA.......")
  message("YOU CANT RUN CICERO HAHAHAHAHAHAHA.......WHAT YOU GONNA DO NOW?? REPEAT EVERYTHING AGAIN????")
  message("Dont worry, drink wome water. I got you boo....keep runing the code XD")
  } else {
  print(seqinfo(hm.integrated.dfree))
  print("Wow...magic")
}

hm.integrated.dfree@assays[["ATAC"]]@seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

if (is.null(seqinfo(hm.integrated.dfree))) {
  message("OOPSssssSSSSssss SOMEONE FORGOT TO RUN ....genome = 'hg38'... WHEN THEY RAN creatChromatinAssay() HAHAHAHAHAHAHA.......")
  message("YOU CANT RUN CICERO HAHAHAHAHAHAHA.......WHAT YOU GONNA DO NOW?? REPEAT EVERYTHING AGAIN????")
  message("Dont worry, drink wome water. I got you boo....keep runing the code XD")
} else {
  print(seqinfo(hm.integrated.dfree))
  print("Wow...magic")
}

# convert to CellDataSet format and make the cicero object
panc.cds <- as.cell_data_set(x = hm.integrated.dfree)
panc.cicero <- make_cicero_cds(panc.cds, reduced_coordinates = reducedDims(panc.cds)$UMAP)

# Cicero connections, bringing accessibility together
genome <- seqlengths(hm.integrated.dfree)

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(panc.cicero, genomic_coords = genome.df, sample_num = 100)
head(conns)

# Find co-acessibility networks
ccans <- generate_ccans(conns)

# Add links to seurat object
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(hm.integrated.dfree) <- links

# View
CoveragePlot(hm.integrated.dfree, region = c("INS", "GCG", "SST", "PPY"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex")

#Save connections takes ages to run this
qsave(conns, file = r"(E:\2.SexbasedStudyCurrent\QS files\conns.qs)")
qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

############################ STAGE ############################
############################   11  ############################
# TF footprinting
# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
hm.integrated.dfree <- AddMotifs(hm.integrated.dfree, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)


# gather the footprinting information for sets of motifs
grep("ARX", hm.integrated.dfree@assays[["ATAC"]]@motifs@motif.names, value = TRUE) 
hm.integrated.dfree <- Footprint(
  object = hm.integrated.dfree,
  motif.name = c("MAFA", "PDX1"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(hm.integrated.dfree, features = c("MAFA"))


# gather the footprinting information for sets of motifs
hm.integrated.dfree <- Footprint(
  object = hm.integrated.dfree,
  motif.name = c("PDX1", "ARX", "MAFA"),
  genome = BSgenome.Hsapiens.UCSC.hg19
)

############################ STAGE ############################
############################   12  ############################
#Trajectory analysis
DefaultAssay(hm.integrated.dfree) <- "ATAC"
pancreas.cds <- as.cell_data_set(hm.integrated.dfree)
pancreas.cds <- cluster_cells(cds = pancreas.cds, reduction_method = "UMAP")
pancreas.cds <- learn_graph(pancreas.cds, use_partition = TRUE)

# order cells
pancreas.cds <- order_cells(pancreas.cds, reduction_method = "UMAP", root_cells = hsc)

# plot trajectories colored by pseudotime
plot_cells(
  cds = pancreas.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

hm.integrated.dfree <- AddMetaData(
  object = hm.integrated.dfree,
  metadata = pancreas.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "monocle"
)

FeaturePlot(hm.integrated.dfree, "monocle", pt.size = 0.1) & scale_color_viridis_c()


############################ STAGE ############################
############################   10  ############################
# DA ANALYSIS
# ALL DATA
#Save file
#qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
table(hm.integrated.dfree@meta.data[["celltype"]])

DefaultAssay(hm.integrated.dfree) <- "ATAC"
Idents(hm.integrated.dfree) <- "celltype"

# DE testing
#plan()
#plan("multisession", workers = 20)
options(future.globals.maxSize = 20 * 1024^3)


beta_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "beta", 
  ident.2 = c("alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(beta_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\beta_peaks.csv)")


alpha_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "alpha", 
  ident.2 = c("beta", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(alpha_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\alpha_peaks.csv)")

delta_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "delta", 
  ident.2 = c("beta", "alpha", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(delta_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\delta_peaks.csv)")

gamma_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "gamma", 
  ident.2 = c("beta", "alpha", "delta", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(gamma_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\gamma_peaks.csv)")

acinar_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "acinar", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(acinar_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\acinar_peaks.csv)")

ductal_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "ductal", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(ductal_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\ductal_peaks.csv)")

activatedstellate_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "activated_stellate", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(activatedstellate_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\activatedstellate_peaks.csv)")

quiescentstellate_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "quiescent_stellate", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(quiescentstellate_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\quiescentstellate_peaks.csv)")

endothelial_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "endothelial", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(endothelial_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\endothelial_peaks.csv)")

macrophage_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "macrophage", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(macrophage_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\macrophage_peaks.csv)")

lymphocyte_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034,
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "lymphocyte", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage"),
  min.pct = 0.05,
  test.use = 'wilcox',
  #latent.vars = 'peak_region_fragments'
)
write.csv(lymphocyte_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex\lymphocyte_peaks.csv)")

############################ STAGE ############################
############################   11  ############################
# DA ANALYSIS
# SEX/ANCESTRY SPECIFIC
# RUN ANALYSIS
system.time({
  ###Step 1: Make Pseudobulk Matrices
  #Read in final Seurat object
  adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
  Idents(adata) <- adata@meta.data$celltype
  samples <- unique(adata@meta.data$sample)
  
  #Pull out list of all cell types
  unique_cell_types <- unique(adata$celltype)
  DefaultAssay(adata) <- 'ATAC'
  
  #Get counts data
  da.counts <- GetAssayData(adata,slot='counts')
  
  dim(da.counts)
  head(da.counts)
  adata_matrices <- adata
  
  ##Pull out barcodes
  sample_bcs <- list()
  for (sample in samples){
    sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$sample == sample,])
  }
  
  lengths(sample_bcs)
  head(sample_bcs[[1]])
  
  #Looping through cell types by making ^ into a function
  get_per_sample_da_SUMS <- function(cell.type, mtx.fp){
    print(paste(cell.type))
    
    #pull out rows of da.counts where BC Ident matches cell.type
    bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
    counts <- da.counts[,colnames(da.counts) %in% bcs]
    print(dim(counts))
    
    #initialize the matrix of sample da
    counts.df <- as.data.frame(rep(0,length(row.names(da.counts))))
    row.names(counts.df) <- row.names(da.counts)
    colnames(counts.df) <- c('temp')
    
    #go through samples and calculate sum of da values
    for (sample in samples){
      sample_cols <- colnames(counts) %in% sample_bcs[[sample]]
      counts.cut <- counts[,sample_cols]
      
      #if only one bc, this becomes a vector which is an issue
      if (typeof(counts.cut) == 'double'){
        mean.counts <- counts.cut
        #if there are NO bcs, this will return NA (just return 0 for everything)
      } else if(length(colnames(counts.cut)) == 0){
        mean.counts <- rep(0,length(row.names(counts)))
      } else {
        mean.counts <- rowSums(counts.cut)
      }
      counts.df <- cbind(counts.df,as.data.frame(mean.counts))
    }
    fin.counts.df <- counts.df[,-c(1)]
    colnames(fin.counts.df) <- samples
    head(fin.counts.df)
    
    #export df
    mtx.fp <- sprintf('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/%s_sample_da_total_counts.txt',cell.type) # change to save dir
    write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
  }
  
  #Run function to make matrices
  unique_cell_types <- unique(adata$celltype)
  for (cell.type in unique_cell_types){
    fp = sprintf('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/%s_pseudobulk.txt',cell.type) # change to save dir as above
    get_per_sample_da_SUMS(cell.type, fp)
  }
  
  
  ###Step 3: DESeq
  #Pseudobulk matrices directory
  dir = 'E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/'
  files = list.files(dir, pattern =".txt")
  cells = gsub("_sample_da_total_counts.txt","", files)
  
  #Create outdir for results
  outdir <- 'E:/2.SexbasedStudyCurrent/test_env/DEtesting/cells/' #changes based on analysis
  
  #Create a metadata table
  meta <- adata@meta.data[,c('sample', 'sex', 'ancestry')]
  colnames(meta) <- c('sample', 'sex', 'ancestry')
  rownames(meta) <- NULL
  meta <- meta[!duplicated(meta),]
  meta$sex_ancestry <- paste(meta$sex, meta$ancestry, sep = '_')
  
  # list of pseudobulk files
  files <- list.files('E:/2.SexbasedStudyCurrent/test_env/pseudobulk/cells/', pattern='da')
  
  ##Create matrices for results
  sumres <- matrix(nrow=length(cells), ncol = 3)
  rownames(sumres) <- cells
  
  # testing for sex_ancestry_diabetes
  for (FILE in files) {
    cell <- gsub('_sample_da_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$sample)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$sample)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('male' %in% meta2$sex && 'female' %in% meta2$sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      accessibilities_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          accessibilities_to_keep <- c(accessibilities_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% accessibilities_to_keep),] 
      
      if ((length(which(meta2$sex == 'male')) > 1) | (length(which(meta2$sex == 'female')) > 1)) { #fix
        my_design <- as.formula ('~sex_ancestry') # alldata
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      }
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('male_white', 'female_white', 'male_white', 'male_black')
      tests2 <- c('male_black', 'female_black', 'female_white', 'female_black')
      
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_ancestry', tests1[[x]],tests2[[x]]) # This should not change when you test subsetted data
        numoft1 <- length(which(meta2$sex_ancestry==t1))
        numoft2 <- length(which(meta2$sex_ancestry==t2))
        
        if (numoft1 < 3 | numoft2 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t1, "INSUFFICIENT N", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft1 > 2 & numoft2 > 2) {
          #sprintf("%s and %s are present in the dataset", t1, t2)
          #sprintf("Find data here: %s", outdir)
          res <- results(dds, contrast=c(test), cooksCutoff=FALSE, independentFiltering=FALSE) #cooksCutoff = FALSE see here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier
          res <- as.data.frame(res)
          res <- res[order(res$pvalue),]
          outfile <- paste0(cell, '.deseq.WaldTest.', tests1[[x]],'.vs.',tests2[[x]],'.tsv')
          write.table(res,paste0(outdir, outfile) , sep='\t', quote=F)
          #print(paste(t1, 'and', t2, 'are present in dataset metadata', sep = " "))
          print(sprintf('%s and %s are present in the dataset metadata', t1, t2)) #just because I wanted to understand using sprintf
          print(paste("Data copied here:", outdir, sep = " "))
          print(paste('####'))
          print(paste('####'))
        }
        
      }
    }
  }
  
  # testing for sex_diabetes
  for (FILE in files) {
    cell <- gsub('_sample_da_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$sample)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$sample)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('male' %in% meta2$sex && 'female' %in% meta2$sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      genes_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
      
      if ((length(which(meta2$sex == 'male')) > 1) | (length(which(meta2$sex == 'female')) > 1)) {
        my_design <- as.formula ('~sex') # design for sex_diabetes
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } 
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('male')
      tests2 <- c('female')
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex', tests1[[x]],tests2[[x]]) # For Sex_diabetes
        numoft1 <- length(which(meta2$sex==t1))
        numoft2 <- length(which(meta2$sex==t2))
        
        
        if (numoft1 > 2 & numoft2 > 2) {
          #sprintf("%s and %s are present in the dataset", t1, t2)
          #sprintf("Find data here: %s", outdir)
          res <- results(dds, contrast=c(test), cooksCutoff=FALSE, independentFiltering=FALSE) #cooksCutoff = FALSE see here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier
          res <- as.data.frame(res)
          res <- res[order(res$pvalue),]
          outfile <- paste0(cell, '.deseq.WaldTest.', tests1[[x]],'.vs.',tests2[[x]],'.tsv')
          write.table(res,paste0(outdir, outfile) , sep='\t', quote=F)
          #print(paste(t1, 'and', t2, 'are present in dataset metadata', sep = " "))
          print(sprintf('%s cells and %s cells are present in the dataset metadata', t1, t2)) #just because I wanted to understand using sprintf
          print(paste("Data copied here:", outdir, sep = " "))
          print(paste('####'))
          print(paste('####'))
        }
        
      }
    }
  }
}) # Sys float time


############################ STAGE ############################
############################   12  ############################
# GO analysis
# All data
system.time({
  dgelist <- list.files(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks\allsex)", # Tulane
                        all.files = FALSE, 
                        full.names = FALSE, 
                        pattern = "*.csv")
  
  # Point towards WD using a function
  for (sample in dgelist){
    wd <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/snATAC/DE_accessible_sites/DA_peaks/allsex/%s', dgelist) # Tulane
  }
  
  # Run iterative function to perform GO on all data
  for (x in wd) {
    sample_name <- str_split_fixed(x, "/", n=13)[13] # needs to be number of / in wd +1 (for tulane = 12 + 1)
    #datfile <- read.table(file.path(x), sep = '\t', row.names = 1) 
    datfile <- read.csv(file.path(x), sep = ",", row.names = 1)
    
    # List of significanta accessibilities
    datfile$p_val_adj[datfile$p_val_adj == 0] <- 2e-302
    datfile <- dplyr::filter(datfile, p_val_adj < 1e-5) 
    datfile <- rownames(datfile[datfile$avg_log2FC > 1, ])
    
    # Closest features NEEDS SIGNA OBJECT
    datfile <- ClosestFeature(hm.integrated.dfree, regions = datfile)
    datfile <- dplyr::filter(datfile, distance < 100000) 
    datfile <- distinct(datfile, gene_name, .keep_all = TRUE)

    # All genes
    datfile <- as.character(datfile$gene_name)
    
    # Run GO enrichment analysis genes up
    GO.up <- enrichGO(gene = datfile, 
                      #universe = all_genes, 
                      keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                      OrgDb = org.Hs.eg.db, 
                      ont = c("ALL"), 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, #if not set default is at 0.05
                      readable = TRUE)
    
    go_data_up <- data.frame(GO.up)
    
    if (nrow(go_data_up) > 0) {
      go_data_up <- dplyr::filter(go_data_up, qvalue < 0.1) }
    
    # Save outputs
    adjusted_name <- gsub('.{4}$', '', sample_name)
    adjusted_name <- gsub('deseq.WaldTest.', '', adjusted_name)
    if (nrow(go_data_up) > 0) {
      write.csv(go_data_up, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/snATAC/ORA/allsex/%s.csv", adjusted_name), row.names = FALSE)} #Tulane
    print(sprintf('%s analysis run', adjusted_name))
  }
})















 
