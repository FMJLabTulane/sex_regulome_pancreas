# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 03/23/2022
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
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("dittoSeq")
BiocManager::install("escape")
BiocManager::install("ComplexHeatmap")
BiocManager::install(c("DropletUtils", "Nebulosa"))
BiocManager::install("hdf5r", force = TRUE)
BiocManager::install("DESeq2")
BiocManager::install('multtest')
BiocManager::install("MAST")

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
    library(scCustomize)
    library(Nebulosa)
    library(DropletUtils)
    library(ggvenn)
    library(ggVennDiagram)
    library(devtools)
    library(R.utils)
    library(qs)
    library(parallel)
    library(readr)
    library(DESeq2)
    library(beeswarm)
    library(limma)
    library(edgeR)
    library(GenomicFeatures)
    library(data.table)
    library(RColorBrewer)
    library(pheatmap)
    library(multtest)
    library(metap)
    library(MAST)
  }
)

# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())


# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")

############################ STAGE ############################
############################   1   ############################
# OBJECT SETUP AND NORMALIZATION ###
# STEP 1: Load 
system.time({
#user    system elapsed 
#8328.56 254.45 8669.53 
  {
    samples <- c("HP2022801", "SAMN15877725", "HP2107001", "HP2107901",
                 "HP2024001", "HP2105501", "HP2108601", "HP2108901",
                 "HP2031401", "HP2110001", "HP2123201",
                 "HP2106201", "HP2121601", "HP2132801", "HP2202101")
  }
  
  {
    for (sample in samples){
      wd <- sprintf('E:/1.SexbasedStudyrawdata/Cellranger_raw_data/scRNAseq/%s', samples)
    }
    }
  {
    # Automated SoupX created
    for (x in wd){
    sample_name <- str_split_fixed(x, "/", n=5)[5] #Adjust this to output your sample name
    rna_counts <- Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    data <- CreateSeuratObject(counts=rna_counts)
    data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
    data <- subset(x = data, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 20)
  
    #Running sctransform takes into account sequencing depth at each cell
    #data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE, vst.flavor = "v2")
    #data <- SCTransform(data, verbose = FALSE)
  
    #Log normalization alternative to sctransform
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
  
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, dims = 1:30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
    data <- FindClusters(data, algorithm=4, resolution = 1, verbose=FALSE)
  
    #Read in RNA assay counts from our filtered seurat object
    DefaultAssay(data) <- 'RNA'
    toc <- GetAssayData(object = data, slot = "counts") #with nFeature >500 filter
    tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5')) #raw count matrix
  
    #Pull out the required metadata from the clustered filtered adata object
    #We need the UMAP coordinates (RD1 and RD2) and the cluster assignments at minimum
    metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                       as.data.frame(Idents(data)),
                       as.data.frame(Idents(data))))
    colnames(metadata) <- c("RD1","RD2","Cluster","Annotation")
  
    #Create the SoupChannel Object
    sc <- SoupChannel(tod,toc)
  
    #Add in the metadata (dimensionality reduction UMAP coords and cluster assignments)
    sc <- setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
    sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc, roundToInt = TRUE) # very important step, use roundToInt = TRUE
  
    #Create post-SoupX Seurat Object
    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    #data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000)  #Can be changed to sctransform
    #data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    #all.genes <- rownames(data2)
    #data2 <- ScaleData(data2, features = all.genes)
    data2 <- SCTransform(data2, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE, vst.flavor = "v2")
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    qsave(data2, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/RDS files/current/1_soupx_rounttoint/%s.qs",sample_name))
    }
  }
})
  

############################ STAGE ############################
############################   2   ############################
# Load data
system.time({
#user     system elapsed 
#30790.28 562.56 31476.69 
{
  HP2022801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2022801.qs)")
  SAMN15877725 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\SAMN15877725.qs)")
  HP2024001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2024001.qs)")
  HP2031401 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2031401.qs)")
  HP2105501 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2105501.qs)")
  HP2106201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2106201.qs)")
  HP2107001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2107001.qs)")
  HP2107901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2107901.qs)")
  HP2108601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2108601.qs)")
  HP2108901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2108901.qs)")
  HP2110001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2110001.qs)")
  HP2121601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2121601.qs)")
  HP2123201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2123201.qs)")
  HP2132801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2132801.qs)")
  HP2202101 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\1_soupx_rounttoint\HP2202101.qs)")
}


# Doublet removal
# Optimization
{
  sweep.res <- paramSweep_v3(HP2022801, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2022801) * 0.05)  # expect 4% doublets
  HP2022801 <- doubletFinder_v3(HP2022801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(SAMN15877725, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(SAMN15877725) * 0.05)  # expect 4% doublets
  SAMN15877725 <- doubletFinder_v3(SAMN15877725, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2024001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2024001) * 0.05)  # expect 4% doublets
  HP2024001 <- doubletFinder_v3(HP2024001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2031401, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2031401) * 0.05)  # expect 4% doublets
  HP2031401 <- doubletFinder_v3(HP2031401, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2105501, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2105501) * 0.05)  # expect 4% doublets
  HP2105501 <- doubletFinder_v3(HP2105501, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2106201, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2106201) * 0.05)  # expect 4% doublets
  HP2106201 <- doubletFinder_v3(HP2106201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2107001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2107001) * 0.05)  # expect 4% doublets
  HP2107001 <- doubletFinder_v3(HP2107001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2107901, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2107901) * 0.05)  # expect 4% doublets
  HP2107901 <- doubletFinder_v3(HP2107901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2108601, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2108601) * 0.05)  # expect 4% doublets
  HP2108601 <- doubletFinder_v3(HP2108601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2108901, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2108901) * 0.05)  # expect 4% doublets
  HP2108901 <- doubletFinder_v3(HP2108901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2110001, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2110001) * 0.05)  # expect 4% doublets
  HP2110001 <- doubletFinder_v3(HP2110001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2121601, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2121601) * 0.05)  # expect 4% doublets
  HP2121601 <- doubletFinder_v3(HP2121601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2123201, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2123201) * 0.05)  # expect 4% doublets
  HP2123201 <- doubletFinder_v3(HP2123201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2132801, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2132801) * 0.05)  # expect 4% doublets
  HP2132801 <- doubletFinder_v3(HP2132801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  sweep.res <- paramSweep_v3(HP2202101, sct = TRUE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats)
  barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
  nExp <- round(ncol(HP2202101) * 0.05)  # expect 4% doublets
  HP2202101 <- doubletFinder_v3(HP2202101, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  }

# Setup one metadata column
HP2022801$doublets <- HP2022801$DF.classifications_0.25_0.09_213
SAMN15877725$doublets <- SAMN15877725$DF.classifications_0.25_0.09_201
HP2107001$doublets <- HP2107001$DF.classifications_0.25_0.09_212
HP2107901$doublets <- HP2107901$DF.classifications_0.25_0.09_165
HP2024001$doublets <- HP2024001$DF.classifications_0.25_0.09_151
HP2105501$doublets <- HP2105501$DF.classifications_0.25_0.09_156
HP2108601$doublets <- HP2108601$DF.classifications_0.25_0.09_273
HP2108901$doublets <- HP2108901$DF.classifications_0.25_0.09_214
HP2031401$doublets <- HP2031401$DF.classifications_0.25_0.09_233
HP2110001$doublets <- HP2110001$DF.classifications_0.25_0.09_290
HP2123201$doublets <- HP2123201$DF.classifications_0.25_0.09_78
HP2106201$doublets <- HP2106201$DF.classifications_0.25_0.09_325
HP2121601$doublets <- HP2121601$DF.classifications_0.25_0.09_179
HP2132801$doublets <- HP2132801$DF.classifications_0.25_0.09_120
HP2202101$doublets <- HP2202101$DF.classifications_0.25_0.09_200

# Step 4: Add cell IDs ####
# Add cell IDs
{
  HP2022801 <- RenameCells(HP2022801, add.cell.id = "HP2022801")
  SAMN15877725 <- RenameCells(SAMN15877725, add.cell.id = "SAMN15877725")
  HP2024001 <- RenameCells(HP2024001, add.cell.id = "HP2024001")
  HP2031401 <- RenameCells(HP2031401, add.cell.id = "HP2031401")
  HP2105501 <- RenameCells(HP2105501, add.cell.id = "HP2105501")
  HP2106201 <- RenameCells(HP2106201, add.cell.id = "HP2106201")
  HP2107001 <- RenameCells(HP2107001, add.cell.id = "HP2107001")
  HP2107901 <- RenameCells(HP2107901, add.cell.id = "HP2107901")
  HP2108601 <- RenameCells(HP2108601, add.cell.id = "HP2108601")
  HP2108901 <- RenameCells(HP2108901, add.cell.id = "HP2108901")
  HP2110001 <- RenameCells(HP2110001, add.cell.id = "HP2110001")
  HP2121601 <- RenameCells(HP2121601, add.cell.id = "HP2121601")
  HP2123201 <- RenameCells(HP2123201, add.cell.id = "HP2123201")
  HP2132801 <- RenameCells(HP2132801, add.cell.id = "HP2132801")
  HP2202101 <- RenameCells(HP2202101, add.cell.id = "HP2202101")
}

# Save point
{
  qsave(HP2022801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2022801.qs)")
  qsave(SAMN15877725, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\SAMN15877725.qs)")
  qsave(HP2024001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2024001.qs)")
  qsave(HP2031401, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2031401.qs)")
  qsave(HP2105501, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2105501.qs)")
  qsave(HP2106201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2106201.qs)")
  qsave(HP2107001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107001.qs)")
  qsave(HP2107901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107901.qs)")
  qsave(HP2108601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108601.qs)")
  qsave(HP2108901, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108901.qs)")
  qsave(HP2110001, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2110001.qs)")
  qsave(HP2121601, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2121601.qs)")
  qsave(HP2123201, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2123201.qs)")
  qsave(HP2132801, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2132801.qs)")
  qsave(HP2202101, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2202101.qs)")
}
})

############################ STAGE ############################
############################   3   ############################

# Analysis of Ruth's data
# Hpap dataset
system.time({
  #user    system elapsed 
  #9294.29 614.52 9908.50 (~2.75hrs)
hpap <- readRDS(file=r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\Raw_hpap_file\hpap_islet_scRNAseq.rds)")
table(hpap@meta.data[["Library"]])

# Subsetting data
Idents(hpap) <- "Library"
{
  HPAP_019 <- subset(hpap, idents = c("HPAP-019"))
  HPAP_020 <- subset(hpap, idents = c("HPAP-020"))
  HPAP_021 <- subset(hpap, idents = c("HPAP-021"))
  HPAP_022 <- subset(hpap, idents = c("HPAP-022"))
  HPAP_023 <- subset(hpap, idents = c("HPAP-023"))
  HPAP_024 <- subset(hpap, idents = c("HPAP-024"))
  HPAP_026 <- subset(hpap, idents = c("HPAP-026"))
  HPAP_028 <- subset(hpap, idents = c("HPAP-028"))
  HPAP_029 <- subset(hpap, idents = c("HPAP-029"))
  HPAP_032 <- subset(hpap, idents = c("HPAP-032"))
  HPAP_034 <- subset(hpap, idents = c("HPAP-034"))
  HPAP_035 <- subset(hpap, idents = c("HPAP-035"))
  HPAP_036 <- subset(hpap, idents = c("HPAP-036"))
  HPAP_037 <- subset(hpap, idents = c("HPAP-037"))
  HPAP_038 <- subset(hpap, idents = c("HPAP-038"))
  HPAP_039 <- subset(hpap, idents = c("HPAP-039"))
  HPAP_040 <- subset(hpap, idents = c("HPAP-040"))
  HPAP_042 <- subset(hpap, idents = c("HPAP-042"))
  HPAP_043 <- subset(hpap, idents = c("HPAP-043"))
  HPAP_044 <- subset(hpap, idents = c("HPAP-044"))
  HPAP_045 <- subset(hpap, idents = c("HPAP-045"))
  HPAP_047 <- subset(hpap, idents = c("HPAP-047"))
  HPAP_049 <- subset(hpap, idents = c("HPAP-049"))
  HPAP_050 <- subset(hpap, idents = c("HPAP-050"))
  HPAP_051 <- subset(hpap, idents = c("HPAP-051"))
  HPAP_052 <- subset(hpap, idents = c("HPAP-052"))
  HPAP_053 <- subset(hpap, idents = c("HPAP-053"))
  HPAP_054 <- subset(hpap, idents = c("HPAP-054"))
  HPAP_055 <- subset(hpap, idents = c("HPAP-055"))
  HPAP_056 <- subset(hpap, idents = c("HPAP-056"))
  HPAP_057 <- subset(hpap, idents = c("HPAP-057"))
  HPAP_058 <- subset(hpap, idents = c("HPAP-058"))
  HPAP_059 <- subset(hpap, idents = c("HPAP-059"))
  HPAP_061 <- subset(hpap, idents = c("HPAP-061"))
  HPAP_063 <- subset(hpap, idents = c("HPAP-063"))
  HPAP_064 <- subset(hpap, idents = c("HPAP-064"))
  HPAP_065 <- subset(hpap, idents = c("HPAP-065"))
  HPAP_070 <- subset(hpap, idents = c("HPAP-070"))
  HPAP_071 <- subset(hpap, idents = c("HPAP-071"))
  HPAP_072 <- subset(hpap, idents = c("HPAP-072"))
  HPAP_074 <- subset(hpap, idents = c("HPAP-074"))
  HPAP_075 <- subset(hpap, idents = c("HPAP-075"))
  HPAP_077 <- subset(hpap, idents = c("HPAP-077"))
  HPAP_079 <- subset(hpap, idents = c("HPAP-079"))
  HPAP_080 <- subset(hpap, idents = c("HPAP-080"))
  HPAP_081 <- subset(hpap, idents = c("HPAP-081"))
  HPAP_082 <- subset(hpap, idents = c("HPAP-082"))
  HPAP_083 <- subset(hpap, idents = c("HPAP-083"))
  HPAP_084 <- subset(hpap, idents = c("HPAP-084"))
  HPAP_085 <- subset(hpap, idents = c("HPAP-085"))
  HPAP_087 <- subset(hpap, idents = c("HPAP-087"))
  HPAP_088 <- subset(hpap, idents = c("HPAP-088"))
  HPAP_090 <- subset(hpap, idents = c("HPAP-090"))
  HPAP_091 <- subset(hpap, idents = c("HPAP-091"))
  HPAP_092 <- subset(hpap, idents = c("HPAP-092"))
  HPAP_099 <- subset(hpap, idents = c("HPAP-099"))
  HPAP_100 <- subset(hpap, idents = c("HPAP-100"))
  HPAP_101 <- subset(hpap, idents = c("HPAP-101"))
  HPAP_103 <- subset(hpap, idents = c("HPAP-103"))
  HPAP_104 <- subset(hpap, idents = c("HPAP-104"))
  HPAP_105 <- subset(hpap, idents = c("HPAP-105"))
  HPAP_106 <- subset(hpap, idents = c("HPAP-106"))
  HPAP_107 <- subset(hpap, idents = c("HPAP-107"))
  HPAP_108 <- subset(hpap, idents = c("HPAP-108"))
  HPAP_109 <- subset(hpap, idents = c("HPAP-109"))
}

pancreas.list <- list("HPAP-019" = HPAP_019, "HPAP-020" = HPAP_020, "HPAP-021" = HPAP_021, "HPAP-022" = HPAP_022, "HPAP-023" = HPAP_023,
                      "HPAP-024" = HPAP_024, "HPAP-026" = HPAP_026, "HPAP-028" = HPAP_028, "HPAP-029" = HPAP_029, "HPAP-032" = HPAP_032,
                      "HPAP-034" = HPAP_034, "HPAP-035" = HPAP_035, "HPAP-036" = HPAP_036, "HPAP-037" = HPAP_037, "HPAP-038" = HPAP_038,
                      "HPAP-039" = HPAP_039, "HPAP-040" = HPAP_040, "HPAP-042" = HPAP_042, "HPAP-043" = HPAP_043, "HPAP-044" = HPAP_044,
                      "HPAP-045" = HPAP_045, "HPAP-047" = HPAP_047, "HPAP-049" = HPAP_049, "HPAP-050" = HPAP_050, "HPAP-051" = HPAP_051,
                      "HPAP-052" = HPAP_052, "HPAP-053" = HPAP_053, "HPAP-054" = HPAP_054, "HPAP-055" = HPAP_055, "HPAP-056" = HPAP_056,
                      "HPAP-057" = HPAP_057, "HPAP-058" = HPAP_058, "HPAP-059" = HPAP_059, "HPAP-061" = HPAP_061, "HPAP-063" = HPAP_063,
                      "HPAP-064" = HPAP_064, "HPAP-065" = HPAP_065, "HPAP-070" = HPAP_070, "HPAP-071" = HPAP_071, "HPAP-072" = HPAP_072,
                      "HPAP-074" = HPAP_074, "HPAP-075" = HPAP_075, "HPAP-077" = HPAP_077, "HPAP-079" = HPAP_079, "HPAP-080" = HPAP_080,
                      "HPAP-081" = HPAP_081, "HPAP-082" = HPAP_082, "HPAP-083" = HPAP_083, "HPAP-084" = HPAP_084, "HPAP-085" = HPAP_085,
                      "HPAP-087" = HPAP_087, "HPAP-088" = HPAP_088, "HPAP-090" = HPAP_090, "HPAP-091" = HPAP_091, "HPAP-092" = HPAP_092,
                      "HPAP-099" = HPAP_099, "HPAP-100" = HPAP_100, "HPAP-101" = HPAP_101, "HPAP-103" = HPAP_103, "HPAP-104" = HPAP_104,
                      "HPAP-105" = HPAP_105, "HPAP-106" = HPAP_106, "HPAP-107" = HPAP_107, "HPAP-108" = HPAP_108, "HPAP-109" = HPAP_109
)

# Pulling data for analysis
# Subset out AAB+, T1D and metadata error individuals
pancreas_subset <- pancreas.list[c("HPAP-022", "HPAP-026", "HPAP-035", "HPAP-036", "HPAP-037", "HPAP-040", "HPAP-051", "HPAP-052", 
                                   "HPAP-053", "HPAP-054", "HPAP-056", "HPAP-057", "HPAP-058", "HPAP-059", "HPAP-061", "HPAP-063", 
                                   "HPAP-065", "HPAP-070", "HPAP-074", "HPAP-075", "HPAP-077", "HPAP-079", "HPAP-080", "HPAP-081",
                                   "HPAP-082", "HPAP-083", "HPAP-085", "HPAP-088", "HPAP-091", "HPAP-099", "HPAP-100", "HPAP-101", 
                                   "HPAP-103", "HPAP-105", "HPAP-106", "HPAP-108", "HPAP-109")]

options(future.globals.maxSize = 1000 * 1024^2)
pancreas_subset <- lapply(X = pancreas_subset, 
                            FUN = function(x){
                              x[['percent.mt']] <- PercentageFeatureSet(x, 
                                                                        pattern = '^MT-')
                              x <- SCTransform(x,
                                               vars.to.regress = "percent.mt", 
                                               verbose = TRUE, 
                                               return.only.var.genes = FALSE, 
                                               vst.flavor = "v2")
                            })
options(future.globals.maxSize = 1000 * 1024^2)

# Read in SoupX corrected/doublet classified data
{
  HP2022801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2022801.qs)")
  SAMN15877725 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\SAMN15877725.qs)")
  HP2024001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2024001.qs)")
  HP2031401 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2031401.qs)")
  HP2105501 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2105501.qs)")
  HP2106201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2106201.qs)")
  HP2107001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107001.qs)")
  HP2107901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2107901.qs)")
  HP2108601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108601.qs)")
  HP2108901 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2108901.qs)")
  HP2110001 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2110001.qs)")
  HP2121601 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2121601.qs)")
  HP2123201 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2123201.qs)")
  HP2132801 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2132801.qs)")
  HP2202101 <- qread(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\2_doublet_defined\HP2202101.qs)")
}


# Adjust and add metadata
# Sex specific Metadata addition
{
  HP2022801$Sex <- "F"
  SAMN15877725$Sex <- "M"
  HP2024001$Sex <- "F"
  HP2031401$Sex <- "M"
  HP2105501$Sex <- "F"
  HP2106201$Sex <- "F"
  HP2107001$Sex <- "M"
  HP2107901$Sex <- "M"
  HP2108601$Sex <- "F"
  HP2108901$Sex <- "F"
  HP2110001$Sex <- "M"
  HP2121601$Sex <- "F"
  HP2123201$Sex <- "M"
  HP2132801$Sex <- "F"
  HP2202101$Sex <- "F"


  HP2022801@meta.data[["Tissue Source"]] <- "Tulane"
  SAMN15877725@meta.data[["Tissue Source"]] <- "Tulane"
  HP2024001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2031401@meta.data[["Tissue Source"]] <- "Tulane"
  HP2105501@meta.data[["Tissue Source"]] <- "Tulane"
  HP2106201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2110001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2121601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2123201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2132801@meta.data[["Tissue Source"]] <- "Tulane"
  HP2202101@meta.data[["Tissue Source"]] <- "Tulane"


HP2022801@meta.data[["Library"]] <- deparse(substitute(HP2022801))
SAMN15877725@meta.data[["Library"]] <- deparse(substitute(SAMN15877725))
HP2024001@meta.data[["Library"]] <- deparse(substitute(HP2024001))
HP2031401@meta.data[["Library"]] <- deparse(substitute(HP2031401))
HP2105501@meta.data[["Library"]] <- deparse(substitute(HP2105501))
HP2106201@meta.data[["Library"]] <- deparse(substitute(HP2106201))
HP2107001@meta.data[["Library"]] <- deparse(substitute(HP2107001))
HP2107901@meta.data[["Library"]] <- deparse(substitute(HP2107901))
HP2108601@meta.data[["Library"]] <- deparse(substitute(HP2108601))
HP2108901@meta.data[["Library"]] <- deparse(substitute(HP2108901))
HP2110001@meta.data[["Library"]] <- deparse(substitute(HP2110001))
HP2121601@meta.data[["Library"]] <- deparse(substitute(HP2121601))
HP2123201@meta.data[["Library"]] <- deparse(substitute(HP2123201))
HP2132801@meta.data[["Library"]] <- deparse(substitute(HP2132801))
HP2202101@meta.data[["Library"]] <- deparse(substitute(HP2202101))


  HP2022801@meta.data[["sample"]] <- NULL
  SAMN15877725@meta.data[["sample"]] <- NULL
  HP2024001@meta.data[["sample"]] <- NULL
  HP2031401@meta.data[["sample"]] <- NULL
  HP2105501@meta.data[["sample"]] <- NULL
  HP2106201@meta.data[["sample"]] <- NULL
  HP2107001@meta.data[["sample"]] <- NULL
  HP2107901@meta.data[["sample"]] <- NULL
  HP2108601@meta.data[["sample"]] <- NULL
  HP2108901@meta.data[["sample"]] <- NULL
  HP2110001@meta.data[["sample"]] <- NULL
  HP2121601@meta.data[["sample"]] <- NULL
  HP2123201@meta.data[["sample"]] <- NULL
  HP2132801@meta.data[["sample"]] <- NULL
  HP2202101@meta.data[["sample"]] <- NULL


  HP2022801@meta.data[["Chemistry"]] <- "10Xv3"
  SAMN15877725@meta.data[["Chemistry"]] <- "10Xv3"
  HP2024001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2031401@meta.data[["Chemistry"]] <- "10Xv3"
  HP2105501@meta.data[["Chemistry"]] <- "10Xv3"
  HP2106201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2110001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2121601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2123201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2132801@meta.data[["Chemistry"]] <- "10Xv3"
  HP2202101@meta.data[["Chemistry"]] <- "10Xv3"


  HP2022801@meta.data[["Diabetes Status"]] <- "ND"
  SAMN15877725@meta.data[["Diabetes Status"]] <- "ND"
  HP2024001@meta.data[["Diabetes Status"]] <- "ND"
  HP2031401@meta.data[["Diabetes Status"]] <- "ND"
  HP2105501@meta.data[["Diabetes Status"]] <- "ND"
  HP2106201@meta.data[["Diabetes Status"]] <- "ND"
  HP2107001@meta.data[["Diabetes Status"]] <- "ND"
  HP2107901@meta.data[["Diabetes Status"]] <- "ND"
  HP2108601@meta.data[["Diabetes Status"]] <- "ND"
  HP2108901@meta.data[["Diabetes Status"]] <- "ND"
  HP2110001@meta.data[["Diabetes Status"]] <- "ND"
  HP2121601@meta.data[["Diabetes Status"]] <- "ND"
  HP2123201@meta.data[["Diabetes Status"]] <- "ND"
  HP2132801@meta.data[["Diabetes Status"]] <- "ND"
  HP2202101@meta.data[["Diabetes Status"]] <- "ND"
}

# Make list
pancreas_qadir <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                       "HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901, 
                       "HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201,
                       "HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
)

# Subset out all single cells
pancreas_qadir <- lapply(X = pancreas_qadir, 
                         FUN = function(x){
                           
                           # Subset out all singlets
                           Idents(x) <- "doublets"
                           x <- subset(x = x, idents = c("Singlet"))
                         })

# Join datasets
pancreas_combined <- c(pancreas_subset, pancreas_qadir)

# merge data sets
sampleset <- names(pancreas_combined)
# sampleset <- c("HPAP-022", "HPAP-026", "HPAP-035", "HPAP-036", "HPAP-037", "HPAP-040", "HPAP-051", "HPAP-052", 
#                "HPAP-053", "HPAP-054", "HPAP-056", "HPAP-057", "HPAP-058", "HPAP-059", "HPAP-061", "HPAP-063", 
#                "HPAP-065", "HPAP-070", "HPAP-074", "HPAP-075", "HPAP-077", "HPAP-079", "HPAP-080", "HPAP-081",
#                "HPAP-082", "HPAP-083", "HPAP-085", "HPAP-088", "HPAP-091", "HPAP-099", "HPAP-100", "HPAP-101", 
#                "HPAP-103", "HPAP-105", "HPAP-106", "HPAP-108", "HPAP-109", 
#                "HP2022801", "SAMN15877725", "HP2107001", "HP2107901", "HP2024001", "HP2105501", "HP2108601", "HP2108901", 
#                "HP2031401", "HP2110001", "HP2123201", "HP2106201", "HP2121601", "HP2132801", "HP2202101")

pancreas_rna <- merge(pancreas_combined[[sampleset[[1]]]], 
                      y=pancreas_combined[sampleset[2:length(sampleset)]], 
                      project='pancreas', 
                      merge.data = TRUE)

# Inspect data
Idents(pancreas_rna) <- "Tissue Source"
VlnPlot(object = pancreas_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Subset
pancreas_rna <- subset(x = pancreas_rna, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 15)

# Inspect data
pancreas_rna
table(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Sex"]])
table(pancreas_rna@meta.data[["Library"]])
table(pancreas_rna@meta.data[["Tissue Source"]])
table(pancreas_rna@meta.data[["Diabetes Status"]])
table(pancreas_rna@meta.data[["ancestry"]])
table(pancreas_rna@meta.data[["doublets"]])
table(pancreas_rna@meta.data[["Cell Type"]])

# Metadata classification, paranoia edition
Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
pancreas_rna$Sex <- NULL
pancreas_rna$Sex <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("F", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
         "M", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
         "M", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
         "F", #Female	16	23	Caucasian	5.2	T1D control	nPOD
         "F", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
         "M", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
         "F", #Female	45.49	43	African American	6.7	T2DM	UPENN
         "M", #Male	38.72	27	African American	5.2	T2D control	UPENN
         "F", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
         "F", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
         "M", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
         "F", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
         "F", #Female	29.26	34	African American	9.4	T2DM	nPOD
         "M", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
         "F", #Female	38.27	59	African American	5.9	T2DM	nPOD
         "F", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
         "M", #Male	37.47	40	African American	9.5	T2DM	nPOD
         "M", #Male	17.09	55	African American	7	T2DM	UPENN
         "F", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
         "M", #Male	27.52	35	Caucasian	6	T2D control	UPENN
         "M", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
         "F", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
         "M", #Male	35.71	22	African American	5.4	T2D control	nPOD
         "F", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
         "M", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
         "M", #Male	35.62	45	African American	5	T2DM	UPENN
         "F", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
         "M", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
         "F", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
         "F", #Female	24.7	28	Hispanic	5	T1D control	UPENN
         "M", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
         "F", #Female	38.01	55	Hispanic	5	T2D control	nPOD
         "F", #Female	36.44	48	Caucasian	6	T2D control	UPENN
         "F", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
         "M", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
         "M", #Male	33	42	African American	10.8	T2DM	nPOD
         "F", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
         "F", #F 25.6 36 W 5.5 ND Tulane
         "M", #M 27.4 31 White N/A ND Tulane
         "M", #M 22.7 39 white 5.1 ND Tulane
         "M", #M 22.6 59 white 5.2 ND Tulane
         "F", #F 23.5 42 white 5.4 ND Tulane
         "F", #F 31.6 36 white 5.6 ND Tulane
         "F", #F 25.6 31 white 5.1 ND Tulane
         "F", #F 30.8 42 white 5.9 ND Tulane
         "M", #M 28.0 37 black 5.4 ND Tulane
         "M", #M 34.3 66 black 5.5 ND Tulane
         "M", #M 29.8 52 black 5.3 ND Tulane
         "F", #F 28.0 41 black 5.3 ND Tulane
         "F", #F 31.8 54 black 5.8 ND Tulane
         "F", #F 35.9 40 black 5.5 ND Tulane
         "F")) #F 23.2 52 black 5.5 ND Tulane
table(pancreas_rna@meta.data[["Sex"]])

Idents(pancreas_rna) <- "Library"
table(pancreas_rna$Library)
pancreas_rna$bmi <- NULL
pancreas_rna$bmi <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("34.7", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "20.8", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "26.91", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "16", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "21.9", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "23.98", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "45.49", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "38.72", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "24.2", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "30.85", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "32.89", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "30.49", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "29.26", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "37.96", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "38.27", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "38.41", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "37.47", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "17.09", #Male	17.09	55	African American	7	T2DM	UPENN
           "36.88", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "27.52", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "32.78", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "28.38", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "35.71", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "28.91", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "23.96", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "35.62", #Male	35.62	45	African American	5	T2DM	UPENN
           "39.78", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "32.81", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "35.58", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "24.7", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "28.83", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "38.01", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "36.44", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "28.1", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "28.12", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "33", #Male	33	42	African American	10.8	T2DM	nPOD
           "29.49", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "25.6", #F 25.6 36 W 5.5 ND Tulane
           "27.4", #M 27.4 31 White N/A ND Tulane
           "22.7", #M 22.7 39 white 5.1 ND Tulane
           "22.6", #M 22.6 59 white 5.2 ND Tulane
           "23.5", #F 23.5 42 white 5.4 ND Tulane
           "31.6", #F 31.6 36 white 5.6 ND Tulane
           "25.6", #F 25.6 31 white 5.1 ND Tulane
           "30.8", #F 30.8 42 white 5.9 ND Tulane
           "28.0", #M 28.0 37 black 5.4 ND Tulane
           "34.3", #M 34.3 66 black 5.5 ND Tulane
           "29.8", #M 29.8 52 black 5.3 ND Tulane
           "28.0", #F 28.0 41 black 5.3 ND Tulane
           "31.8", #F 31.8 54 black 5.8 ND Tulane
           "35.9", #F 35.9 40 black 5.5 ND Tulane
           "23.2")) #F 23.2 52 black 5.5 ND Tulane))
unique(pancreas_rna$bmi)

Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
pancreas_rna$age <- NULL
pancreas_rna$age <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("39",
         "24",
         "35",
         "23",
         "35",
         "35",
         "43",
         "27",
         "58",
         "40",
         "33",
         "50",
         "34",
         "35",
         "59",
         "45",
         "40",
         "55",
         "40",
         "35",
         "47",
         "52",
         "22",
         "45",
         "25",
         "45",
         "48",
         "37",
         "50",
         "28",
         "41",
         "55",
         "48",
         "51",
         "55",
         "42",
         "59",
         "36",
         "31",
         "39",
         "59",
         "42",
         "36",
         "31",
         "42",
         "37",
         "66",
         "52",
         "41",
         "54",
         "40",
         "52"))
unique(pancreas_rna$age)

Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
pancreas_rna$ancestry <- NULL
pancreas_rna$ancestry <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("white", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
                "white", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
                "white", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
                "white", #Female	16	23	Caucasian	5.2	T1D control	nPOD
                "white", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
                "white", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
                "black", #Female	45.49	43	African American	6.7	T2DM	UPENN
                "black", #Male	38.72	27	African American	5.2	T2D control	UPENN
                "white", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
                "white", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
                "white", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
                "white", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
                "black", #Female	29.26	34	African American	9.4	T2DM	nPOD
                "white", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
                "black", #Female	38.27	59	African American	5.9	T2DM	nPOD
                "white", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
                "black", #Male	37.47	40	African American	9.5	T2DM	nPOD
                "black", #Male	17.09	55	African American	7	T2DM	UPENN
                "white", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
                "white", #Male	27.52	35	Caucasian	6	T2D control	UPENN
                "white", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
                "hispanic", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
                "black", #Male	35.71	22	African American	5.4	T2D control	nPOD
                "white", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
                "white", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
                "black", #Male	35.62	45	African American	5	T2DM	UPENN
                "white", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
                "white", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
                "hispanic", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
                "hispanic", #Female	24.7	28	Hispanic	5	T1D control	UPENN
                "white", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
                "hispanic", #Female	38.01	55	Hispanic	5	T2D control	nPOD
                "white", #Female	36.44	48	Caucasian	6	T2D control	UPENN
                "hispanic", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
                "white", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
                "black", #Male	33	42	African American	10.8	T2DM	nPOD
                "hispanic", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
                "white", #F 25.6 36 White 5.5 ND Tulane
                "white", #M 27.4 31 White N/A ND Tulane
                "white", #M 22.7 39 white 5.1 ND Tulane
                "white", #M 22.6 59 white 5.2 ND Tulane
                "white", #F 23.5 42 white 5.4 ND Tulane
                "white", #F 31.6 36 white 5.6 ND Tulane
                "white", #F 25.6 31 white 5.1 ND Tulane
                "white", #F 30.8 42 white 5.9 ND Tulane
                "black", #M 28.0 37 black 5.4 ND Tulane
                "black", #M 34.3 66 black 5.5 ND Tulane
                "black", #M 29.8 52 black 5.3 ND Tulane
                "black", #F 28.0 41 black 5.3 ND Tulane
                "black", #F 31.8 54 black 5.8 ND Tulane
                "black", #F 35.9 40 black 5.5 ND Tulane
                "black")) #F 23.2 52 black 5.5 ND Tulane
unique(pancreas_rna$ancestry) 

Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
pancreas_rna$'Tissue Source' <- NULL
pancreas_rna$tissue_source <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("UPENN", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
         "nPOD", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
         "UPENN", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
         "nPOD", #Female	16	23	Caucasian	5.2	T1D control	nPOD
         "UPENN", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
         "UPENN", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
         "UPENN", #Female	45.49	43	African American	6.7	T2DM	UPENN
         "UPENN", #Male	38.72	27	African American	5.2	T2D control	UPENN
         "UPENN", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
         "UPENN", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
         "UPENN", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
         "UPENN", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
         "nPOD", #Female	29.26	34	African American	9.4	T2DM	nPOD
         "UPENN", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
         "nPOD", #Female	38.27	59	African American	5.9	T2DM	nPOD
         "nPOD", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
         "nPOD", #Male	37.47	40	African American	9.5	T2DM	nPOD
         "UPENN", #Male	17.09	55	African American	7	T2DM	UPENN
         "UPENN", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
         "UPENN", #Male	27.52	35	Caucasian	6	T2D control	UPENN
         "UPENN", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
         "nPOD", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
         "nPOD", #Male	35.71	22	African American	5.4	T2D control	nPOD
         "nPOD", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
         "nPOD", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
         "UPENN", #Male	35.62	45	African American	5	T2DM	UPENN
         "UPENN", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
         "nPOD", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
         "nPOD", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
         "UPENN", #Female	24.7	28	Hispanic	5	T1D control	UPENN
         "nPOD", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
         "nPOD", #Female	38.01	55	Hispanic	5	T2D control	nPOD
         "UPENN", #Female	36.44	48	Caucasian	6	T2D control	UPENN
         "nPOD", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
         "UPENN", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
         "nPOD", #Male	33	42	African American	10.8	T2DM	nPOD
         "nPOD", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
         "Tulane", #F 25.6 36 W 5.5 ND Tulane
         "Tulane", #M 27.4 31 White N/A ND Tulane
         "Tulane", #M 22.7 39 white 5.1 ND Tulane
         "Tulane", #M 22.6 59 white 5.2 ND Tulane
         "Tulane", #F 23.5 42 white 5.4 ND Tulane
         "Tulane", #F 31.6 36 white 5.6 ND Tulane
         "Tulane", #F 25.6 31 white 5.1 ND Tulane
         "Tulane", #F 30.8 42 white 5.9 ND Tulane
         "Tulane", #M 28.0 37 black 5.4 ND Tulane
         "Tulane", #M 34.3 66 black 5.5 ND Tulane
         "Tulane", #M 29.8 52 black 5.3 ND Tulane
         "Tulane", #F 28.0 41 black 5.3 ND Tulane
         "Tulane", #F 31.8 54 black 5.8 ND Tulane
         "Tulane", #F 35.9 40 black 5.5 ND Tulane
         "Tulane")) #F 23.2 52 black 5.5 ND Tulane
unique(pancreas_rna$tissue_source)

Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
pancreas_rna$diabetes_status <- plyr::mapvalues(
  x= pancreas_rna$Library,
  from = c("HPAP-022", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
           "HPAP-026", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
           "HPAP-035", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
           "HPAP-036", #Female	16	23	Caucasian	5.2	T1D control	nPOD
           "HPAP-037", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
           "HPAP-040", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
           "HPAP-051", #Female	45.49	43	African American	6.7	T2DM	UPENN
           "HPAP-052", #Male	38.72	27	African American	5.2	T2D control	UPENN
           "HPAP-053", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
           "HPAP-054", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
           "HPAP-056", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
           "HPAP-057", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
           "HPAP-058", #Female	29.26	34	African American	9.4	T2DM	nPOD
           "HPAP-059", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
           "HPAP-061", #Female	38.27	59	African American	5.9	T2DM	nPOD
           "HPAP-063", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
           "HPAP-065", #Male	37.47	40	African American	9.5	T2DM	nPOD
           "HPAP-070", #Male	17.09	55	African American	7	T2DM	UPENN
           "HPAP-074", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
           "HPAP-075", #Male	27.52	35	Caucasian	6	T2D control	UPENN
           "HPAP-077", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
           "HPAP-079", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
           "HPAP-080", #Male	35.71	22	African American	5.4	T2D control	nPOD
           "HPAP-081", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
           "HPAP-082", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
           "HPAP-083", #Male	35.62	45	African American	5	T2DM	UPENN
           "HPAP-085", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
           "HPAP-088", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
           "HPAP-091", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
           "HPAP-099", #Female	24.7	28	Hispanic	5	T1D control	UPENN
           "HPAP-100", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
           "HPAP-101", #Female	38.01	55	Hispanic	5	T2D control	nPOD
           "HPAP-103", #Female	36.44	48	Caucasian	6	T2D control	UPENN
           "HPAP-105", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
           "HPAP-106", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
           "HPAP-108", #Male	33	42	African American	10.8	T2DM	nPOD
           "HPAP-109", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
           "HP2022801", #F 25.6 36 W 5.5 ND Tulane
           "SAMN15877725", #M 27.4 31 White N/A ND Tulane
           "HP2107001", #M 22.7 39 white 5.1 ND Tulane
           "HP2107901", #M 22.6 59 white 5.2 ND Tulane
           "HP2024001", #F 23.5 42 white 5.4 ND Tulane
           "HP2105501", #F 31.6 36 white 5.6 ND Tulane
           "HP2108601", #F 25.6 31 white 5.1 ND Tulane
           "HP2108901", #F 30.8 42 white 5.9 ND Tulane
           "HP2031401", #M 28.0 37 black 5.4 ND Tulane
           "HP2110001", #M 34.3 66 black 5.5 ND Tulane
           "HP2123201", #M 29.8 52 black 5.3 ND Tulane
           "HP2106201", #F 28.0 41 black 5.3 ND Tulane
           "HP2121601", #F 31.8 54 black 5.8 ND Tulane
           "HP2132801", #F 35.9 40 black 5.5 ND Tulane
           "HP2202101"), #F 23.2 52 black 5.5 ND Tulane
  to = c("ND", #Female	34.7	39	Caucasian	4.7	T1D control	UPENN
         "ND", #Male	20.8	24	Caucasian	4.9	T1D control	nPOD
         "ND", #Male	26.91	35	Caucasian	5.2	T1D control	UPENN
         "ND", #Female	16	23	Caucasian	5.2	T1D control	nPOD
         "ND", #Female	21.9	35	Caucasian	5.3	T1D control	UPENN
         "ND", #Male	23.98	35	Caucasian	5.4	T1D control	UPENN
         "T2D", #Female	45.49	43	African American	6.7	T2DM	UPENN
         "ND", #Male	38.72	27	African American	5.2	T2D control	UPENN
         "ND", #Female	24.2	58	Caucasian	5.6	T2D control	UPENN
         "ND", #Female	30.85	40	Caucasian	4.8	T2D control	UPENN
         "ND", #Male	32.89	33	Caucasian	5.6	T1D control	UPENN
         "T2D", #Female	30.49	50	Caucasian	5.8	T2DM	UPENN
         "T2D", #Female	29.26	34	African American	9.4	T2DM	nPOD
         "ND", #Male	37.96	35	Caucasian	5.1	T2D control	UPENN
         "T2D", #Female	38.27	59	African American	5.9	T2DM	nPOD
         "ND", #Female	38.41	45	Caucasian	6.3	T2D control	nPOD
         "T2D", #Male	37.47	40	African American	9.5	T2DM	nPOD
         "T2D", #Male	17.09	55	African American	7	T2DM	UPENN
         "ND", #Female	36.88	40	Caucasian	6.3	T2D control	UPENN
         "ND", #Male	27.52	35	Caucasian	6	T2D control	UPENN
         "ND", #Male	32.78	47	Caucasian	5.7	T2D control	UPENN
         "T2D", #Female	28.38	52	Hispanic	6.8	T2DM	nPOD
         "ND", #Male	35.71	22	African American	5.4	T2D control	nPOD
         "T2D", #Female	28.91	45	Caucasian	7.9	T2DM	nPOD
         "ND", #Male	23.96	25	Caucasian	5.6	T1D control	nPOD
         "T2D", #Male	35.62	45	African American	5	T2DM	UPENN
         "T2D", #Female	39.78	48	Caucasian	7.3	T2DM	UPENN
         "T2D", #Male	32.81	37	Caucasian	10.3	T2DM	nPOD
         "T2D", #Female	35.58	50	Hispanic	6.9	T2DM	nPOD
         "ND", #Female	24.7	28	Hispanic	5	T1D control	UPENN
         "T2D", #Male	28.83	41	Caucasian	10.7	T2DM	nPOD
         "ND", #Female	38.01	55	Hispanic	5	T2D control	nPOD
         "ND", #Female	36.44	48	Caucasian	6	T2D control	UPENN
         "ND", #Female	28.1	51	Hispanic	5.2	T2D control	nPOD
         "T2D", #Male	28.12	55	Caucasian	5.8	T2DM	UPENN
         "T2D", #Male	33	42	African American	10.8	T2DM	nPOD
         "T2D", #Female	29.49	59	Hispanic	7.5	T2DM	nPOD
         "ND", #F 25.6 36 W 5.5 ND Tulane
         "ND", #M 27.4 31 White N/A ND Tulane
         "ND", #M 22.7 39 white 5.1 ND Tulane
         "ND", #M 22.6 59 white 5.2 ND Tulane
         "ND", #F 23.5 42 white 5.4 ND Tulane
         "ND", #F 31.6 36 white 5.6 ND Tulane
         "ND", #F 25.6 31 white 5.1 ND Tulane
         "ND", #F 30.8 42 white 5.9 ND Tulane
         "ND", #M 28.0 37 black 5.4 ND Tulane
         "ND", #M 34.3 66 black 5.5 ND Tulane
         "ND", #M 29.8 52 black 5.3 ND Tulane
         "ND", #F 28.0 41 black 5.3 ND Tulane
         "ND", #F 31.8 54 black 5.8 ND Tulane
         "ND", #F 35.9 40 black 5.5 ND Tulane
         "ND")) #F 23.2 52 black 5.5 ND Tulane
unique(pancreas_rna$diabetes_status)

# Remove Unecessary metadata
pancreas_rna@meta.data[["sex"]] <- NULL
pancreas_rna$`Diabetes Status` <- NULL
pancreas_rna$`Tissue Source` <- NULL
pancreas_rna$celltype_sex_ancestry_disease <- NULL
pancreas_rna$celltype_sex_disease <- NULL
pancreas_rna$celltype_sex_ancestry_disease <- NULL
pancreas_rna$celltype_sex_ancestry <- NULL
pancreas_rna$celltype_sex <- NULL
pancreas_rna$ancestry_sex <- NULL

# #Split data on basis of disease status and calculate integration features on object list
# DefaultAssay(pancreas_rna) <- "SCT"
# integrationfeatures <- SelectIntegrationFeatures(pancreas_combined, nfeatures = 3000, verbose = TRUE)
# 
# # replacing variable features with integrationfeatures
# VariableFeatures(pancreas_rna, assay = "SCT") <- integrationfeatures
# VariableFeatures(pancreas_rna, assay = "RNA") <- integrationfeatures


#Perform basic threshold filtering and log normalization
DefaultAssay(pancreas_rna) <- "RNA"
pancreas_rna <- NormalizeData(pancreas_rna, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas_rna <- FindVariableFeatures(pancreas_rna, selection.method = "vst", nfeatures = 2000)
pancreas_rna <- ScaleData(pancreas_rna, verbose = FALSE) %>% 
  RunPCA(pc.genes = pancreas_rna@assays$RNA@var.features, npcs = 20, verbose = FALSE)

# # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
# Idents(pancreas_rna) <- "ancestry"
# pancreas_rna$ancestry_sex <- paste(Idents(pancreas_rna), pancreas_rna$Sex, sep = "_")
# table(pancreas_rna@meta.data[["ancestry_sex"]])

#Run Harmony batch correction with library and tissue source covariates
Idents(pancreas_rna) <- "Library"
unique(pancreas_rna$Library)
unique(pancreas_rna$tissue_source)
pancreas_rna <- RunHarmony(pancreas_rna, 
                           assay.use = "RNA",
                           reduction = "pca",
                           dims.use = 1:20,
                           group.by.vars = c("Library", "Chemistry", "tissue_source"),
                           kmeans_init_nstart=20, kmeans_init_iter_max=100,
                           plot_convergence = TRUE)

unique(pancreas_rna@meta.data[["Library"]])
unique(pancreas_rna@meta.data[["Tissue Source"]])
unique(pancreas_rna@meta.data[["ancestry_sex"]])
unique(pancreas_rna@meta.data[["Cell Type"]])
unique(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Chemistry"]])

# Run UMAP
pancreas_rna <- RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=TRUE)

# Clustering
pancreas_rna <- pancreas_rna %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(6), method = 'igraph') #25 res

# Save file
qsave(pancreas_rna, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\pancreas_rna.qs)")
qsave(pancreas_rna, file = r"(E:\2.SexbasedStudyCurrent\QS files\pancreas_rna.qs)")
})


############################ STAGE ############################
############################   4   ############################

system.time({
# Load data
pancreas_rna <- qread(file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\pancreas_rna.qs)")

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
#clustree(pancreas_rna, prefix = "RNA_snn_res.")

# View clustering
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=TRUE)
FeaturePlot(object = pancreas_rna,
            features = c("TPSAB1"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 10000,
            slot = 'counts',
            order = TRUE,
            raster=TRUE)

# Subclustering
Idents(pancreas_rna) <- "RNA_snn_res.6"
subset_clust <- subset(pancreas_rna, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", 
                                                "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", 
                                                "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", 
                                                #"51", 
                                                "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", 
                                                #"62", #"63", 
                                                "64", "65", "66", "67", #"68", 
                                                "69", "70", "71", #"72", 
                                                "73", "74", "75", "76", #"77", 
                                                "78", "79", "80", "81", "82", "83", #"84", 
                                                #"85", 
                                                "86"))

# Checking cluster loss
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=TRUE)
DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=TRUE)

# # As cells were subsetted re-run batch correction and cluster assignmnet
# subset_clust <- RunHarmony(subset_clust, 
#                            assay.use = "RNA",
#                            reduction = "pca",
#                            dims.use = 1:20,
#                            group.by.vars = c("Library", "Chemistry", "tissue_source"),
#                            kmeans_init_nstart=20, kmeans_init_iter_max=100,
#                            plot_convergence = TRUE)
# 
# 
# # UMAP 
# subset_clust <- RunUMAP(subset_clust, reduction = "harmony", dims = 1:20, return.model = TRUE)
# 
# #Neighbours + Clustering
# subset_clust <- subset_clust %>% 
#   FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
#   FindClusters(algorithm=4,resolution = c(6), method = 'igraph')
# 
# DefaultAssay(subset_clust) <- "RNA"
# DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
# FeaturePlot(object = subset_clust,
#             features = c("HSPB1", "DNAJB6", "HSPH1", "GADD45B"
#             ),
#             pt.size = 0.01,
#             cols = c("darkgrey", "red"),
#             min.cutoff = 100,
#             max.cutoff = 500,
#             slot = 'counts',
#             order = TRUE,
#             raster=TRUE)

DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = subset_clust,
            features = c("PPY"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

# Cluster assignment
table(subset_clust@meta.data$RNA_snn_res.6)
Idents(subset_clust) <- "RNA_snn_res.6"
subset_clust <- RenameIdents(subset_clust, 
                             "1" = "alpha",
                             "2" = "ductal",
                             "3" = "alpha",
                             "4" = "beta+alpha",
                             "5" = "beta",
                             "6" = "alpha",
                             "7" = "beta+alpha",
                             "8" = "alpha",
                             "9" = "acinar",
                             "10" = "alpha",
                             "11" = "beta",
                             "12" = "acinar",
                             "13" = "alpha",
                             "14" = "quiescent_stellate",
                             "15" = "alpha",
                             "16" = "ductal",
                             "17" = "alpha",
                             "18" = "ductal",
                             "19" = "beta",
                             "20" = "beta",
                             "21" = "activated_stellate",
                             "22" = "acinar",
                             "23" = "beta",
                             "24" = "acinar",
                             "25" = "activated_stellate",
                             "26" = "acinar",
                             "27" = "beta",
                             "28" = "beta",
                             "29" = "alpha",
                             "30" = "alpha",
                             "31" = "alpha",
                             "32" = "endothelial",
                             "33" = "alpha",
                             "34" = "delta",
                             "35" = "endothelial",
                             "36" = "beta",
                             "37" = "alpha",
                             "38" = "alpha",
                             "39" = "beta+delta",
                             "40" = "beta",
                             "41" = "acinar",
                             "42" = "alpha",
                             "43" = "beta",
                             "44" = "alpha",
                             "45" = "acinar",
                             "46" = "gamma",
                             "47" = "beta+alpha",
                             "48" = "alpha",
                             "49" = "beta",
                             "50" = "activated_stellate",
                             #"51" = "",
                             "52" = "alpha",
                             "53" = "beta",
                             "54" = "acinar",
                             "55" = "acinar",
                             "56" = "alpha",
                             "57" = "lymphocytes",
                             "58" = "macrophages",
                             "59" = "acinar",
                             "60" = "beta",
                             "61" = "ductal",
                             #"62" = "",
                             #"63" = "",
                             "64" = "activated_stellate",
                             "65" = "ductal",
                             "66" = "endothelial",
                             "67" = "delta",
                             #"68" = "",
                             "69" = "alpha",
                             "70" = "acinar",
                             "71" = "ductal",
                             #"72" = "",
                             "73" = "cycling_endo",
                             "74" = "acinar",
                             "75" = "endothelial",
                             "76" = "beta+alpha",
                             #"77" = "",
                             "78" = "acinar",
                             "79" = "acinar",
                             "80" = "endothelial",
                             "81" = "macrophages",
                             "82" = "acinar",
                             "83" = "schwann",
                             #"84" = "",
                             #"85" = "",
                             "86" = "acinar"
                             )

# Check renaming
table(subset_clust@active.ident)
unique(subset_clust@active.ident)


############################## #
########## SUBSET1 ########### #
############################## #

# Subset epsilon cells
DefaultAssay(object = subset_clust) <- "RNA"
Idents(subset_clust, WhichCells(object = subset_clust, expression = GHRL > 50, slot = 'counts')) <- 'epsilon'
subset_clust$celltype_qadir <- Idents(subset_clust)
Idents(subset_clust) <- "celltype_qadir"
DimPlot(subset_clust, reduction = "umap", label = TRUE)

############################## #
########## SUBSET2 ########### #
############################## #

# Subsetting Lymphocytes from Mast cells
# First subset all cells, they are called "Mast" in the primary object
lymphocytes <- subset(subset_clust, idents = "lymphocytes")
table(lymphocytes@active.ident)
lymphocytes <- lymphocytes %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0.5), method = 'igraph')

# Run UMAP
lymphocytes <- RunUMAP(lymphocytes, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(lymphocytes, reduction = 'umap', label = FALSE, pt.size = 1, raster=FALSE)
# Cluster assignment
table(lymphocytes@meta.data$RNA_snn_res.0.5)
Idents(lymphocytes) <- "RNA_snn_res.0.5"
DimPlot(lymphocytes, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE, pt.size = 1, raster=FALSE)
FeaturePlot(object = lymphocytes,
            features = c("TPSAB1"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

table(lymphocytes@active.ident)
lymphocytes <- RenameIdents(lymphocytes, 
                     "1" = "mast",
                     "2" = "lymphocyte",
                     "3" = "mast",
                     "4" = "lymphocyte"
)

table(lymphocytes@active.ident)

# Generate a new column called celltype_qadir in the metadata copying all Ident info there, this is active idents so check
table(subset_clust@active.ident)
subset_clust$celltype_qadir <- as.character(Idents(subset_clust)) #as.character imp
table(subset_clust$celltype_qadir)

# Change the information of cells containing sub-cluster information
subset_clust$celltype_qadir[Cells(lymphocytes)] <- paste(Idents(lymphocytes))
table(subset_clust$celltype_qadir)
DimPlot(subset_clust, 
        #split.by = "Tissue Source", 
        group.by = "celltype_qadir", 
        label = FALSE, 
        ncol = 1,  
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange",       #ductal
                 "salmon3",          #acinar
                 "orangered",        #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "beta+alpha", "alpha", "cycling_endo", "epsilon", "gamma", "delta", "beta+delta",
               "ductal", "acinar", 
               "activated_stellate", "quiescent_stellate", "endothelial",
               "macrophages", "lymphocyte", "mast",
               "schwann")
table(subset_clust@meta.data$celltype_qadir)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
subset_clust@meta.data$celltype_qadir <- factor(x = subset_clust@meta.data$celltype_qadir, levels = my_levels)
table(subset_clust@meta.data$celltype_qadir)
unique(subset_clust@meta.data$celltype_qadir)

# Set celltype_qadir as default
Idents(subset_clust) <- "celltype_qadir"

# Change the information of cells containing sub-cluster information
DimPlot(subset_clust, 
        #split.by = "Diabetes Status", 
        group.by = "celltype_qadir", 
        label = TRUE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.5,
        cols = c("dodgerblue3",      #beta
                 "turquoise2",       #beta+alpha
                 "lightseagreen",    #alpha
                 "darkseagreen2",    #cycling_endo
                 "khaki2",           #epsilon 
                 "springgreen4",     #gamma
                 "chartreuse3",      #delta
                 "burlywood3",       #beta+delta
                 "darkorange",       #ductal
                 "salmon3",          #acinar
                 "orangered",        #activated_setallate
                 "salmon",           #quiescent_stellate
                 "red",              #endothelial
                 "magenta3",         #macrophages
                 "orchid1",          #lymphocytes
                 "red4",             #mast
                 "grey30"            #schwann
        )
)


# Create a metadata slot for celltype_sex, celltype_sex_ancestry and celltype_sex_ancestry_disease
processed_rna <- subset_clust
Idents(processed_rna) <- "celltype_qadir"
processed_rna$celltype_sex <- paste(Idents(processed_rna), processed_rna$Sex, sep = "_")
Idents(processed_rna) <- "celltype_sex"
processed_rna$celltype_sex_ancestry <- paste(Idents(processed_rna), processed_rna$ancestry, sep = "_")
Idents(processed_rna) <- "celltype_sex_ancestry"
processed_rna$celltype_sex_ancestry_disease <- paste(Idents(processed_rna), processed_rna$'diabetes_status', sep = "_")
Idents(processed_rna) <- "celltype_sex"
processed_rna$celltype_sex_disease <- paste(Idents(processed_rna), processed_rna$'diabetes_status', sep = "_")
table(processed_rna$celltype_qadir)
table(processed_rna$celltype_sex)
table(processed_rna$celltype_sex_ancestry)
table(processed_rna$celltype_sex_ancestry_disease)
table(processed_rna$celltype_sex_disease)

# Save data
#qsave(processed_rna, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\current\3_seuratobj\processed_rna.qs)")
#qsave(processed_rna, file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
})



############################ STAGE ############################
############################   5   ############################

# after performing QC on the dataset we discover that some sex strata are incorrect.
# Correcting metadata, for Deseq2 analysis
# some individuals had incorrect metadata
# adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\old_data\processed_rna.qs)")
# 
# # Subset data
# Idents(adata) <- "Library"
# adata <- subset(adata, idents = c("HPAP-022", "HPAP-026", "HPAP-035", "HPAP-036", "HPAP-037", "HPAP-040",   
#                                   #"HPAP-044", 
#                                   "HPAP-051", "HPAP-052", "HPAP-053", "HPAP-054", "HPAP-056",   
#                                   "HPAP-057", "HPAP-058", "HPAP-059", "HPAP-061", "HPAP-063", "HPAP-065",   
#                                   "HPAP-070", "HPAP-074", "HPAP-075", "HPAP-077", "HPAP-079", "HPAP-080",    
#                                   "HPAP-081", "HPAP-082", "HPAP-083", "HPAP-085", "HPAP-088", #"HPAP-090",   
#                                   "HPAP-091", "HPAP-099", "HPAP-100", "HPAP-101", "HPAP-103", "HPAP-105",   
#                                   "HPAP-106", "HPAP-108", "HPAP-109", "HP2022801", "SAMN15877725", "HP2107001",  
#                                   "HP2107901", "HP2024001", "HP2105501", "HP2108601", "HP2108901", "HP2031401",    
#                                   "HP2110001", "HP2123201", "HP2106201", "HP2121601", "HP2132801", "HP2202101"))
# 
# # As cells were subsetted re-run batch correction and cluster assignmnet
# subset_clust <- RunHarmony(adata, 
#                            assay.use = "RNA",
#                            reduction = "pca",
#                            dims.use = 1:20,
#                            group.by.vars = c('Library','Tissue Source', 'Chemistry'),
#                            kmeans_init_nstart=20, kmeans_init_iter_max=100,
#                            plot_convergence = TRUE)
# 
# # UMAP 
# subset_clust <- RunUMAP(subset_clust, reduction = "harmony", dims = 1:20, return.model = TRUE)
# 
# # Celltype annotations should remain consistent
# DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
# DimPlot(subset_clust, reduction = 'umap', group.by = 'celltype_qadir', label = TRUE, pt.size = 0.01, raster=FALSE)
# 
# # Change the information of cells containing sub-cluster information
# DimPlot(subset_clust, 
#         #split.by = "Diabetes Status", 
#         group.by = "celltype_qadir", 
#         label = TRUE, 
#         ncol = 1, 
#         raster = FALSE,
#         pt.size = 0.5,
#         cols = c("dodgerblue3",      #beta
#                  "turquoise2",       #beta+alpha
#                  "lightseagreen",    #alpha
#                  "darkseagreen2",    #cycling_endo
#                  "khaki2",           #epsilon 
#                  "springgreen4",     #gamma
#                  "chartreuse3",      #delta
#                  "burlywood3",       #beta+delta
#                  "darkorange",       #ductal
#                  "salmon3",          #acinar
#                  "orangered",        #activated_setallate
#                  "salmon",           #quiescent_stellate
#                  "red",              #endothelial
#                  "magenta3",         #macrophages
#                  "orchid1",          #lymphocytes
#                  "red4",             #mast
#                  "grey30"            #schwann
#         )
# )
# 
# 
# # Create a metadata slot for celltype_sex, celltype_sex_ancestry and celltype_sex_ancestry_disease
# processed_rna <- subset_clust
# Idents(processed_rna) <- "celltype_qadir"
# processed_rna$celltype_sex <- paste(Idents(processed_rna), processed_rna$Sex, sep = "_")
# Idents(processed_rna) <- "celltype_sex"
# processed_rna$celltype_sex_ancestry <- paste(Idents(processed_rna), processed_rna$ancestry, sep = "_")
# Idents(processed_rna) <- "celltype_sex_ancestry"
# processed_rna$celltype_sex_ancestry_disease <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
# Idents(processed_rna) <- "celltype_sex"
# processed_rna$celltype_sex_disease <- paste(Idents(processed_rna), processed_rna$'Diabetes Status', sep = "_")
# table(processed_rna$celltype_qadir)
# table(processed_rna$celltype_sex)
# table(processed_rna$celltype_sex_ancestry)
# table(processed_rna$celltype_sex_ancestry_disease)
# table(processed_rna$celltype_sex_disease)
# 
# # Save file
# qsave(processed_rna, r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

############################ STAGE ############################
############################   6   ############################
# ALL DATA
# RUN ANALYSIS
system.time({
###Step 1: Make Pseudobulk Matrices
#Read in final Seurat object
adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(adata) <- "Tissue Source"
Idents(adata) <- adata@meta.data$celltype_qadir
samples <- unique(adata@meta.data$Library)

#Pull out list of all cell types
unique_cell_types <- unique(adata$celltype_qadir)
DefaultAssay(adata) <- 'RNA'

#Get counts data
gex.counts <- GetAssayData(adata,slot='counts')

dim(gex.counts)
head(gex.counts)
adata_matrices <- adata

##Pull out barcodes
sample_bcs <- list()
for (sample in samples){
  sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$Library == sample,])
}

lengths(sample_bcs)
head(sample_bcs[[1]])

#Looping through cell types by making ^ into a function
get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
  print(paste(cell.type))
  
  #pull out rows of gex.counts where BC Ident matches cell.type
  bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
  counts <- gex.counts[,colnames(gex.counts) %in% bcs]
  print(dim(counts))
  
  #initialize the matrix of sample gex
  counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
  row.names(counts.df) <- row.names(gex.counts)
  colnames(counts.df) <- c('temp')
  
  #go through samples and calculate sum of gex values
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
  mtx.fp <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/alldata/%s_sample_gex_total_counts.txt',cell.type) # change to save dir
  write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
}

#Run function to make matrices
unique_cell_types <- unique(adata$celltype_qadir)
for (cell.type in unique_cell_types){
  fp = sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/alldata/%s_pseudobulk.txt',cell.type) # change to save dir as above
  get_per_sample_gex_SUMS(cell.type, fp)
}

###Step 2: Make TPM Matrices
#Pull out gene exon info and calculate effective length
gene_annotations_gtf_fp <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.annotation.gtf'
suppressMessages(txdb <- makeTxDbFromGFF(gene_annotations_gtf_fp,format="gtf"))
exons.list.per.gene <- exonsBy(txdb,by="gene") #Collect the exons per gene_id
#Reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

#checking <- gene.info
gene.info <- rtracklayer::import('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.basic.annotation.gtf')
gene.info <- as.data.frame(gene.info)
gene.info <- gene.info %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
head(gene.info)

# Calculate entire gene boundaries
gene.info.start <- gene.info %>% group_by(gene_id) %>% slice_min(order_by = start)
gene.info.start <- gene.info.start[!duplicated(gene.info.start$gene_id),]
gene.info.start <- gene.info.start %>% select(gene_id, gene_name, gene_type, seqnames, end, strand, source, level)

gene.info.end <- gene.info %>% group_by(gene_id) %>% slice_max(order_by = end)
gene.info.end <- gene.info.end[!duplicated(gene.info.end$gene_id),]
gene.info.end <- gene.info.end %>% select(gene_id, start)

gene.info.comp <- merge(gene.info.end, gene.info.start, by = "gene_id")
gene.info.comp <- gene.info.comp %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
gene.info.comp$check <- ifelse(gene.info.comp$end > gene.info.comp$start, 'TRUE',
                               ifelse(gene.info.comp$end < gene.info.comp$start, 'FALSE'))

unique(gene.info.comp$check)
gene.info <- gene.info.comp

#Add the effective lengths to the original gene.info dataframe
temp_df <- gene.info
rownames(temp_df) <- gene.info$gene_id
temp_df2 <- as.data.frame(exonic.gene.sizes)
temp_df2$gene_id <- rownames(temp_df2)

new_df <- merge(temp_df,temp_df2, by='gene_id', all=TRUE)
#Remove duplicate rows from gene info df
fin.gene.info <- new_df[!duplicated(new_df$gene_name),]

#Read in psedobulk matrices from above 
dir = 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/alldata/'

files = list.files(dir, pattern =".txt")
cells = gsub("_sample_gex_total_counts.txt","", files)

make_tpm = function(raw_counts, gene_sizes){
  rpk <- raw_counts / gene_sizes
  tpm <- rpk
  for (i in 1:ncol(rpk)){
    tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
    
  }
  return(tpm)
}

#Output dir for TPM matrices
outdir = "C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/TPM/alldata/"

for (FILE in files){
  cell <- cells[which(files == FILE)]
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  raw_counts <- subset(raw_counts ,rownames(raw_counts) %in% fin.gene.info$gene_name)
  gene_sizes <- fin.gene.info$exonic.gene.sizes[match(rownames(raw_counts), fin.gene.info$gene_name )]
  
  tpm_mat <- make_tpm(raw_counts, gene_sizes)
  write.table(tpm_mat, paste0(outdir,  cell, "_TPM_per_sample.txt"), sep="\t", quote=F)
}


###Step 3: DESeq
#Create a metadata table
meta <- adata@meta.data[,c('Library', 'Sex', 'tissue_source', 'Chemistry', 'ancestry', 'diabetes_status')]
colnames(meta) <- c('Library', 'Sex', 'Tissue_Source', 'Chemistry', 'ancestry', 'Diabetes_Status')
rownames(meta) <- NULL
meta <- meta[!duplicated(meta),]
meta$sex_ancestry_diabetes <- paste0(meta$Sex, '_', meta$ancestry, '_', meta$Diabetes_Status)
meta$sex_diabetes <- paste0(meta$Sex, '_', meta$Diabetes_Status)

#Pseudobulk matrices directory
dir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/alldata/'

#Create outdir for results
outdir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/alldata/' #changes based on analysis
#dir.create(outdir) #works like mkdir

# list of pseudobulk files
files <- list.files('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/alldata/', pattern='gex')

##Create matrices for results
sumres <- matrix(nrow=length(cells), ncol = 3)
rownames(sumres) <- cells

# testing for sex_ancestry_diabetes
for (FILE in files) {
  cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  sample_names <- unique(adata@meta.data$Library)
  sample_names <- gsub('-','.', sample_names)
  raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
  raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
  meta$Library2 <- gsub('-', '.', meta$Library)
  meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
  
  if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
    print(cell)
    print('Data for 2 sex present, however not all data may be present will check this at a later step')
    
    genes_to_keep <- c()
    for (i in 1:nrow(raw_counts)) {
      if (sum(raw_counts[i, ] >= 5) >= 2) {
        genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
      }
    }
    counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
    
    if (length(which(meta2$Chemistry == '10Xv2')) > 1 &
        length(which(meta2$Tissue_Source == 'UPENN')) > 1 & #Alldata 
        length(which(meta2$Tissue_Source == 'nPOD')) > 1 &
        length(which(meta2$Tissue_Source == 'Tulane')) > 1) {
      my_design <- as.formula ('~Chemistry + Tissue_Source + sex_ancestry_diabetes') # alldata
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    } else if (length(which(meta2$Tissue_Source == 'UPENN')) > 1 & #Alldata 
               length(which(meta2$Tissue_Source == 'nPOD')) > 1 &
               length(which(meta2$Tissue_Source == 'Tulane')) > 1) {
      my_design <- as.formula ('~Tissue_Source + sex_ancestry_diabetes') # alldata
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    } else {
      my_design <- as.formula ('~sex_ancestry_diabetes') # alldata
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500)}
    
    # No need for conditional formatting here
    # Specifying test combinations
    tests1 <- c('M_white_ND', 'M_white_ND', 'M_black_ND', 'F_white_ND', 'F_white_ND', 'F_black_ND', 'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'M_white_T2D', 
                'M_white_T2D', 'M_black_T2D', 'F_white_T2D', 'F_white_T2D', 'F_black_T2D', 'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 
                'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D')
    
    tests2 <- c('M_hispanic_ND', 'M_black_ND', 'M_hispanic_ND', 'F_hispanic_ND', 'F_black_ND', 'F_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND', 'M_hispanic_T2D', 
                'M_black_T2D', 'M_hispanic_T2D', 'F_hispanic_T2D', 'F_black_T2D', 'F_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D', 
                'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND')
    
    print('Preparing to run DESeq2')
    
    for (x in 1:length(tests1)){
      # No need for conditional formatting here
      t1 <- tests1[[x]]
      t2 <- tests2[[x]]
      test <- c('sex_ancestry_diabetes', tests1[[x]],tests2[[x]]) # This should not change when you test subsetted data
      numoft1 <- length(which(meta2$sex_ancestry_diabetes==t1))
      numoft2 <- length(which(meta2$sex_ancestry_diabetes==t2))
      
      if (numoft1 < 3) {
        message(paste("!!WARNING!!"))
        message(paste(t1, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste('####'))
        message(paste('####'))
      } else if (numoft2 < 3) {
        message(paste("!!WARNING!!"))
        message(paste(t2, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste("####"))
        message(paste("####"))
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
  cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
  raw_counts <- read.table(paste0(dir, FILE), row.names=1)
  sample_names <- unique(adata@meta.data$Library)
  sample_names <- gsub('-','.', sample_names)
  raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
  raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
  meta$Library2 <- gsub('-', '.', meta$Library)
  meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
  
  if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
    print(cell)
    print('Data for 2 sex present, however not all data may be present will check this at a later step')
    
    genes_to_keep <- c()
    for (i in 1:nrow(raw_counts)) {
      if (sum(raw_counts[i, ] >= 5) >= 2) {
        genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
      }
    }
    counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
   
     if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && #Alldata 
         length(which(meta2$sex_diabetes == 'M_T2D')) > 1 && 
         length(which(meta2$sex_diabetes == 'F_ND')) > 1 && 
         length(which(meta2$sex_diabetes == 'F_T2D')) > 1) {
      my_design <- as.formula ('~Chemistry + Tissue_Source + sex_diabetes') # design for sex_diabetes
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    } else {
      print(sprintf('%s does not have sufficient diabetes samples to test, bypassing to test ND only', cell))
      meta2 <- subset(meta2, Diabetes_Status == 'ND') # it is possible some T2D are present so eliminate them from your dataset since you are restricted to sex
      counts <- counts[,meta2$Library2]
      my_design <- as.formula ('~Tissue_Source + sex_diabetes') # design for sex_diabetes
      dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
    }
    
    # No need for conditional formatting here
    # Specifying test combinations
    tests1 <- c('M_ND', 'M_T2D', 'F_T2D', 'M_T2D')
    tests2 <- c('F_ND', 'M_ND', 'F_ND', 'F_T2D')
    print('Preparing to run DESeq2')
    
    for (x in 1:length(tests1)){
      # No need for conditional formatting here
      t1 <- tests1[[x]]
      t2 <- tests2[[x]]
      test <- c('sex_diabetes', tests1[[x]],tests2[[x]]) # For Sex_diabetes
      numoft1 <- length(which(meta2$sex_diabetes==t1))
      numoft2 <- length(which(meta2$sex_diabetes==t2))
      
      
      if (numoft1 < 3) {
        message(paste("!!WARNING CHECK METADATA!!"))
        message(paste(t1, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste('####'))
        message(paste('####'))
      } else if (numoft2 < 3) {
        message(paste("!!WARNING CHECK METADATA!!"))
        message(paste(t2, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
        message(paste("####"))
        message(paste("####"))
      } else if (numoft1 > 2 & numoft2 > 2) {
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
############################   7   ############################
# Gene ontology analysis Rapid Gene ontology Auto Loader (Rapid GOAL)
# Create a list of all files in directory
system.time({
  dgelist <- list.files(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata)", 
                        all.files = FALSE, 
                        full.names = FALSE, 
                        pattern = "*.tsv")
  
  # Point towards WD using a function
  for (sample in dgelist){
    wd <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/alldata/%s', dgelist)
    }
  
  # Run iterative function to perform GO on all data
  for (x in wd) {
    sample_name <- str_split_fixed(x, "/", n=13)[13] # needs to be number of / in wd +1 (for alldata = 12 + 1)
    datfile <- read.table(file.path(x), sep = '\t', row.names = 1) 
    #datfile <- read.csv(file.path(x), row.names = 1)
    
    # Gene list of genes going UP
    sig_df_up <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange > 0.000000000014)
    sig_genes_up <- rownames(sig_df_up)
    
    # Gene list of genes going DOWN
    sig_df_down <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange < -0.000000000014) 
    sig_genes_down <- rownames(sig_df_down)
    
    # All genes
    all_genes <- rownames(datfile)
    
    # Run GO enrichment analysis genes up
    GO.up <- enrichGO(gene = sig_genes_up, 
                      universe = all_genes, 
                      keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                      OrgDb = org.Hs.eg.db, 
                      ont = c("ALL"), 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, #if not set default is at 0.05
                      readable = TRUE)
    
    # Run GO enrichment analysis genes down
    GO.down <- enrichGO(gene = sig_genes_down, 
                        universe = all_genes, 
                        keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                        OrgDb = org.Hs.eg.db, 
                        ont = c("ALL"), 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 1, 
                        qvalueCutoff = 1, #if not set default is at 0.05
                        readable = TRUE)
    
    go_data_up <- data.frame(GO.up)
    go_data_down <- data.frame(GO.down)
    
    if (nrow(go_data_up) > 0) {
    go_data_up <- dplyr::filter(go_data_up, qvalue < 0.1) }
    if (nrow(go_data_down) > 0) {
    go_data_down <- dplyr::filter(go_data_down, qvalue < 0.1)}
    
    # Save outputs
    adjusted_name <- gsub('.{4}$', '', sample_name)
    adjusted_name <- gsub('deseq.WaldTest.', '', adjusted_name)
    if (nrow(go_data_up) > 0) {
    write.csv(go_data_up, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/alldata/UP/%s.csv", adjusted_name), row.names = FALSE)}
    if (nrow(go_data_down) > 0) {
    write.csv(go_data_down, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/alldata/DOWN/%s.csv", adjusted_name), row.names = FALSE)}
    }
}) # Sys float time

############################ STAGE ############################
############################   8   ############################
# TULANE TEST
###Step 1: Make Pseudobulk Matrices
#Read in final Seurat object
# RUN ANALYSIS
system.time({
#user    system elapsed 
#3098.69 60.33  3371.21 
adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(adata) <- "Tissue Source"
tulane <- subset(adata, idents = c("Tulane")) # tulane
adata <- tulane #tulane
Idents(adata) <- adata@meta.data$celltype_qadir
samples <- unique(adata@meta.data$Library)

  #Pull out list of all cell types
  unique_cell_types <- unique(adata$celltype_qadir)
  
  DefaultAssay(adata) <- 'RNA'
  #Get counts data
  gex.counts <- GetAssayData(adata,slot='counts')
  
  dim(gex.counts)
  head(gex.counts)
  adata_matrices <- adata
  
  ##Pull out barcodes
  sample_bcs <- list()
  for (sample in samples){
    sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$Library == sample,])
  }
  
  lengths(sample_bcs)
  head(sample_bcs[[1]])
  
  #Looping through cell types by making ^ into a function
  get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
    print(paste(cell.type))
    
    #pull out rows of gex.counts where BC Ident matches cell.type
    bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
    counts <- gex.counts[,colnames(gex.counts) %in% bcs]
    print(dim(counts))
    
    #initialize the matrix of sample gex
    counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
    row.names(counts.df) <- row.names(gex.counts)
    colnames(counts.df) <- c('temp')
    
    #go through samples and calculate sum of gex values
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
    mtx.fp <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/tulane/%s_sample_gex_total_counts.txt',cell.type) # Tulane
    write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
  }
  
  #Run function to make matrices
  unique_cell_types <- unique(adata$celltype_qadir)
  for (cell.type in unique_cell_types){
    fp = sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/tulane/%s_pseudobulk.txt',cell.type) # Tulane
    get_per_sample_gex_SUMS(cell.type, fp)
  }
  
  ###Step 2: Make TPM Matrices
  #Pull out gene exon info and calculate effective length
  gene_annotations_gtf_fp <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.annotation.gtf'
  suppressMessages(txdb <- makeTxDbFromGFF(gene_annotations_gtf_fp,format="gtf"))
  exons.list.per.gene <- exonsBy(txdb,by="gene") #Collect the exons per gene_id
  #Reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
  exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
  
  #checking <- gene.info
  gene.info <- rtracklayer::import('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.basic.annotation.gtf')
  gene.info <- as.data.frame(gene.info)
  gene.info <- gene.info %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
  head(gene.info)
  
  # Calculate entire gene boundaries
  gene.info.start <- gene.info %>% group_by(gene_id) %>% slice_min(order_by = start)
  gene.info.start <- gene.info.start[!duplicated(gene.info.start$gene_id),]
  gene.info.start <- gene.info.start %>% select(gene_id, gene_name, gene_type, seqnames, end, strand, source, level)
  
  gene.info.end <- gene.info %>% group_by(gene_id) %>% slice_max(order_by = end)
  gene.info.end <- gene.info.end[!duplicated(gene.info.end$gene_id),]
  gene.info.end <- gene.info.end %>% select(gene_id, start)
  
  gene.info.comp <- merge(gene.info.end, gene.info.start, by = "gene_id")
  gene.info.comp <- gene.info.comp %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
  gene.info.comp$check <- ifelse(gene.info.comp$end > gene.info.comp$start, 'TRUE',
                                 ifelse(gene.info.comp$end < gene.info.comp$start, 'FALSE'))
  
  unique(gene.info.comp$check)
  gene.info <- gene.info.comp
  
  #Add the effective lengths to the original gene.info dataframe
  temp_df <- gene.info
  rownames(temp_df) <- gene.info$gene_id
  temp_df2 <- as.data.frame(exonic.gene.sizes)
  temp_df2$gene_id <- rownames(temp_df2)
  
  new_df <- merge(temp_df,temp_df2, by='gene_id', all=TRUE)
  #Remove duplicate rows from gene info df
  fin.gene.info <- new_df[!duplicated(new_df$gene_name),]
  
  #Read in psedobulk matrices from above 
  dir = 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/tulane/' # Tulane
  
  files = list.files(dir, pattern =".txt")
  cells = gsub("_sample_gex_total_counts.txt","", files)
  
  make_tpm = function(raw_counts, gene_sizes){
    rpk <- raw_counts / gene_sizes
    tpm <- rpk
    for (i in 1:ncol(rpk)){
      tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
      
    }
    return(tpm)
  }
  
  #Output dir for TPM matrices
  outdir = "C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/TPM/tulane/" # Tulane
  
  for (FILE in files){
    cell <- cells[which(files == FILE)]
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    raw_counts <- subset(raw_counts ,rownames(raw_counts) %in% fin.gene.info$gene_name)
    gene_sizes <- fin.gene.info$exonic.gene.sizes[match(rownames(raw_counts), fin.gene.info$gene_name )]
    
    tpm_mat <- make_tpm(raw_counts, gene_sizes)
    write.table(tpm_mat, paste0(outdir,  cell, "_TPM_per_sample.txt"), sep="\t", quote=F)
  }

###Step 3: DESeq
  #Create a metadata table
  meta <- adata@meta.data[,c('Library', 'Sex', 'Tissue Source', 'Chemistry', 'ancestry', 'Diabetes Status')]
  colnames(meta) <- c('Library', 'Sex', 'Tissue_Source', 'Chemistry', 'ancestry', 'Diabetes_Status')
  rownames(meta) <- NULL
  meta <- meta[!duplicated(meta),]
  meta$sex_ancestry_diabetes <- paste0(meta$Sex, '_', meta$ancestry, '_', meta$Diabetes_Status)
  meta$sex_diabetes <- paste0(meta$Sex, '_', meta$Diabetes_Status)
  
  #Pseudobulk matrices directory
  dir <- "C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/tulane/" # Tulane
  
  #Create outdir for results
  outdir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/tulane/' #Tulane
  #dir.create(outdir) #works like mkdir
  
  # list of pseudobulk files
  files <- list.files('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/tulane/', pattern='gex') #Tulane
  
  ##Create matrices for results
  sumres <- matrix(nrow=length(cells), ncol = 3)
  rownames(sumres) <- cells
  
  # testing for sex_ancestry_diabetes
  for (FILE in files) {
    cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$Library)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$Library)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      genes_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
      
      if (length(unique(meta2$Chemistry)) > 1) {
        my_design <- as.formula ('~sex_ancestry_diabetes') # Tulane
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else {
        my_design <- as.formula ('~sex_ancestry_diabetes') # Tulane
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      }
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('M_white_ND', 'M_white_ND', 'M_black_ND', 'F_white_ND', 'F_white_ND', 'F_black_ND', 'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'M_white_T2D', 
                  'M_white_T2D', 'M_black_T2D', 'F_white_T2D', 'F_white_T2D', 'F_black_T2D', 'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 
                  'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D')
      
      tests2 <- c('M_hispanic_ND', 'M_black_ND', 'M_hispanic_ND', 'F_hispanic_ND', 'F_black_ND', 'F_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND', 'M_hispanic_T2D', 
                  'M_black_T2D', 'M_hispanic_T2D', 'F_hispanic_T2D', 'F_black_T2D', 'F_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D', 
                  'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND')
      
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_ancestry_diabetes', tests1[[x]],tests2[[x]]) # This should not change when you test subsetted data
        numoft1 <- length(which(meta2$sex_ancestry_diabetes==t1))
        numoft2 <- length(which(meta2$sex_ancestry_diabetes==t2))
        
        if (numoft1 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t1, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft2 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t2, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste("####"))
          message(paste("####"))
        } else if (numoft1 > 2 & numoft2 > 2) {
          #sprintf("%s and %s are present in the dataset", t1, t2)
          #sprintf("Find data here: %s", outdir)
          res <- results(dds, contrast=c(test), cooksCutoff=FALSE) #cooksCutoff = FALSE see here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier
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
    cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$Library)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$Library)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      genes_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
      
      if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && #Tulane 
         length(which(meta2$sex_diabetes == 'F_ND')) > 1) {
        my_design <- as.formula ('~sex_diabetes') # Tulane
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else {
        print(sprintf('%s does not have sufficient diabetes samples to test, bypassing to test ND only', cell))
        meta2 <- subset(meta2, Diabetes_Status == 'ND') # it is possible some T2D are present so eliminate them from your dataset since you are restricted to sex
        counts <- counts[,meta2$Library2]
        my_design <- as.formula ('~sex_diabetes') # Tulane
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      }
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('M_ND', 'M_T2D', 'F_T2D', 'M_T2D')
      tests2 <- c('F_ND', 'M_ND', 'F_ND', 'F_T2D')
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_diabetes', tests1[[x]],tests2[[x]]) # For Sex_diabetes
        numoft1 <- length(which(meta2$sex_diabetes==t1))
        numoft2 <- length(which(meta2$sex_diabetes==t2))

        
        if (numoft1 < 3) {
          message(paste("!!WARNING CHECK METADATA!!"))
          message(paste(t1, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft2 < 3) {
          message(paste("!!WARNING CHECK METADATA!!"))
          message(paste(t2, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste("####"))
          message(paste("####"))
        } else if (numoft1 > 2 & numoft2 > 2) {
          #sprintf("%s and %s are present in the dataset", t1, t2)
          #sprintf("Find data here: %s", outdir)
          res <- results(dds, contrast=c(test), cooksCutoff=FALSE) #cooksCutoff = FALSE see here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier
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
############################   9   ############################
# Gene ontology analysis Rapid Gene ontology Auto Loader (Rapid GOAL)
# Create a list of all files in directory
system.time({
dgelist <- list.files(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\tulane)", # Tulane
                      all.files = FALSE, 
                      full.names = FALSE, 
                      pattern = "*.tsv")

# Point towards WD using a function
for (sample in dgelist){
  wd <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/tulane/%s', dgelist) # Tulane
}

# Run iterative function to perform GO on all data
for (x in wd) {
  sample_name <- str_split_fixed(x, "/", n=13)[13] # needs to be number of / in wd +1 (for tulane = 12 + 1)
  datfile <- read.table(file.path(x), sep = '\t', row.names = 1) 
  #datfile <- read.csv(file.path(x), row.names = 1)
  
  # Gene list of genes going UP
  sig_df_up <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange > 0.000000000014)
  sig_genes_up <- rownames(sig_df_up)
  
  # Gene list of genes going DOWN
  sig_df_down <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange < -0.000000000014) 
  sig_genes_down <- rownames(sig_df_down)
  
  # All genes
  all_genes <- rownames(datfile)
  
  # Run GO enrichment analysis genes up
  GO.up <- enrichGO(gene = sig_genes_up, 
                    universe = all_genes, 
                    keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                    OrgDb = org.Hs.eg.db, 
                    ont = c("ALL"), 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1, #if not set default is at 0.05
                    readable = TRUE)
  
  # Run GO enrichment analysis genes down
  GO.down <- enrichGO(gene = sig_genes_down, 
                      universe = all_genes, 
                      keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                      OrgDb = org.Hs.eg.db, 
                      ont = c("ALL"), 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, #if not set default is at 0.05
                      readable = TRUE)
  
  go_data_up <- data.frame(GO.up)
  go_data_down <- data.frame(GO.down)
  
  if (nrow(go_data_up) > 0) {
    go_data_up <- dplyr::filter(go_data_up, qvalue < 0.1) }
  if (nrow(go_data_down) > 0) {
    go_data_down <- dplyr::filter(go_data_down, qvalue < 0.1)}
  
  # Save outputs
  adjusted_name <- gsub('.{4}$', '', sample_name)
  adjusted_name <- gsub('deseq.WaldTest.', '', adjusted_name)
  if (nrow(go_data_up) > 0) {
  write.csv(go_data_up, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/tulane/UP/%s.csv", adjusted_name), row.names = FALSE)} #Tulane
  if (nrow(go_data_down) > 0) {
  write.csv(go_data_down, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/tulane/DOWN/%s.csv", adjusted_name), row.names = FALSE)}
  print(sprintf('%s analysis run', adjusted_name))
}
}) #Sys float time

############################ STAGE ############################
############################   10  ############################
# RUN ANALYSIS
system.time({
# ALL DATA
###Step 1: Make Pseudobulk Matrices
#Read in final Seurat object
adata <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
Idents(adata) <- "Tissue Source"
hpap <- subset(adata, idents = c("nPod", "UPenn")) # tulane
adata <- hpap #hpap
Idents(adata) <- adata@meta.data$celltype_qadir
samples <- unique(adata@meta.data$Library)

  #Pull out list of all cell types
  unique_cell_types <- unique(adata$celltype_qadir)
  DefaultAssay(adata) <- 'RNA'
  
  #Get counts data
  gex.counts <- GetAssayData(adata,slot='counts')
  
  dim(gex.counts)
  head(gex.counts)
  adata_matrices <- adata
  
  ##Pull out barcodes
  sample_bcs <- list()
  for (sample in samples){
    sample_bcs[[sample]] <- row.names(adata[[]][adata[[]]$Library == sample,])
  }
  
  lengths(sample_bcs)
  head(sample_bcs[[1]])
  
  #Looping through cell types by making ^ into a function
  get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
    print(paste(cell.type))
    
    #pull out rows of gex.counts where BC Ident matches cell.type
    bcs <- names(Idents(adata_matrices)[Idents(adata_matrices) == cell.type])
    counts <- gex.counts[,colnames(gex.counts) %in% bcs]
    print(dim(counts))
    
    #initialize the matrix of sample gex
    counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
    row.names(counts.df) <- row.names(gex.counts)
    colnames(counts.df) <- c('temp')
    
    #go through samples and calculate sum of gex values
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
    mtx.fp <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/hpap/%s_sample_gex_total_counts.txt',cell.type) # change to save dir
    write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
  }
  
  #Run function to make matrices
  unique_cell_types <- unique(adata$celltype_qadir)
  for (cell.type in unique_cell_types){
    fp = sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/hpap/%s_pseudobulk.txt',cell.type) # change to save dir as above
    get_per_sample_gex_SUMS(cell.type, fp)
  }
  
  ###Step 2: Make TPM Matrices
  #Pull out gene exon info and calculate effective length
  gene_annotations_gtf_fp <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.annotation.gtf'
  suppressMessages(txdb <- makeTxDbFromGFF(gene_annotations_gtf_fp,format="gtf"))
  exons.list.per.gene <- exonsBy(txdb,by="gene") #Collect the exons per gene_id
  #Reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
  exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
  
  #checking <- gene.info
  gene.info <- rtracklayer::import('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/publicdata/gencode_v38/gencode.v38.basic.annotation.gtf')
  gene.info <- as.data.frame(gene.info)
  gene.info <- gene.info %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
  head(gene.info)
  
  # Calculate entire gene boundaries
  gene.info.start <- gene.info %>% group_by(gene_id) %>% slice_min(order_by = start)
  gene.info.start <- gene.info.start[!duplicated(gene.info.start$gene_id),]
  gene.info.start <- gene.info.start %>% select(gene_id, gene_name, gene_type, seqnames, end, strand, source, level)
  
  gene.info.end <- gene.info %>% group_by(gene_id) %>% slice_max(order_by = end)
  gene.info.end <- gene.info.end[!duplicated(gene.info.end$gene_id),]
  gene.info.end <- gene.info.end %>% select(gene_id, start)
  
  gene.info.comp <- merge(gene.info.end, gene.info.start, by = "gene_id")
  gene.info.comp <- gene.info.comp %>% select(gene_id, gene_name, gene_type, seqnames, start, end, strand, source, level)
  gene.info.comp$check <- ifelse(gene.info.comp$end > gene.info.comp$start, 'TRUE',
                                 ifelse(gene.info.comp$end < gene.info.comp$start, 'FALSE'))
  
  unique(gene.info.comp$check)
  gene.info <- gene.info.comp
  
  #Add the effective lengths to the original gene.info dataframe
  temp_df <- gene.info
  rownames(temp_df) <- gene.info$gene_id
  temp_df2 <- as.data.frame(exonic.gene.sizes)
  temp_df2$gene_id <- rownames(temp_df2)
  
  new_df <- merge(temp_df,temp_df2, by='gene_id', all=TRUE)
  #Remove duplicate rows from gene info df
  fin.gene.info <- new_df[!duplicated(new_df$gene_name),]
  
  #Read in psedobulk matrices from above 
  dir = 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/hpap/'
  
  files = list.files(dir, pattern =".txt")
  cells = gsub("_sample_gex_total_counts.txt","", files)
  
  make_tpm = function(raw_counts, gene_sizes){
    rpk <- raw_counts / gene_sizes
    tpm <- rpk
    for (i in 1:ncol(rpk)){
      tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
      
    }
    return(tpm)
  }
  
  #Output dir for TPM matrices
  outdir = "C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/TPM/hpap/"
  
  for (FILE in files){
    cell <- cells[which(files == FILE)]
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    raw_counts <- subset(raw_counts ,rownames(raw_counts) %in% fin.gene.info$gene_name)
    gene_sizes <- fin.gene.info$exonic.gene.sizes[match(rownames(raw_counts), fin.gene.info$gene_name )]
    
    tpm_mat <- make_tpm(raw_counts, gene_sizes)
    write.table(tpm_mat, paste0(outdir,  cell, "_TPM_per_sample.txt"), sep="\t", quote=F)
  }


###Step 3: DESeq
  #Create a metadata table
  meta <- adata@meta.data[,c('Library', 'Sex', 'Tissue Source', 'Chemistry', 'ancestry', 'Diabetes Status')]
  colnames(meta) <- c('Library', 'Sex', 'Tissue_Source', 'Chemistry', 'ancestry', 'Diabetes_Status')
  rownames(meta) <- NULL
  meta <- meta[!duplicated(meta),]
  meta$sex_ancestry_diabetes <- paste0(meta$Sex, '_', meta$ancestry, '_', meta$Diabetes_Status)
  meta$sex_diabetes <- paste0(meta$Sex, '_', meta$Diabetes_Status)
  
  #Pseudobulk matrices directory
  dir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/hpap/'
  
  #Create outdir for results
  outdir <- 'C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/hpap/' #changes based on analysis
  #dir.create(outdir) #works like mkdir
  
  # list of pseudobulk files
  files <- list.files('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/pseudobulk_counts/hpap/', pattern='gex')
  
  ##Create matrices for results
  sumres <- matrix(nrow=length(cells), ncol = 3)
  rownames(sumres) <- cells
  
  # testing for sex_ancestry_diabetes
  for (FILE in files) {
    cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$Library)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$Library)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      genes_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
      
      if (length(unique(meta2$Chemistry)) > 1) {
        my_design <- as.formula ('~Chemistry + Tissue_Source + sex_ancestry_diabetes') # alldata
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else if (length(which(meta2$sex_ancestry_diabetes == 'M_white_ND')) > 1 || #Alldata 
                 length(which(meta2$sex_ancestry_diabetes == 'F_white_ND')) > 1) {
        my_design <- as.formula ('~Tissue_Source + sex_ancestry_diabetes') # alldata
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else {
        message(paste("!!WARNING!!"))
        message(paste(cell, "cell containing samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))}
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('M_white_ND', 'M_white_ND', 'M_black_ND', 'F_white_ND', 'F_white_ND', 'F_black_ND', 'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'M_white_T2D', 
                  'M_white_T2D', 'M_black_T2D', 'F_white_T2D', 'F_white_T2D', 'F_black_T2D', 'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 
                  'M_white_T2D', 'M_black_T2D', 'M_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D')
      
      tests2 <- c('M_hispanic_ND', 'M_black_ND', 'M_hispanic_ND', 'F_hispanic_ND', 'F_black_ND', 'F_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND', 'M_hispanic_T2D', 
                  'M_black_T2D', 'M_hispanic_T2D', 'F_hispanic_T2D', 'F_black_T2D', 'F_hispanic_T2D', 'F_white_T2D', 'F_black_T2D', 'F_hispanic_T2D', 
                  'M_white_ND', 'M_black_ND', 'M_hispanic_ND', 'F_white_ND', 'F_black_ND', 'F_hispanic_ND')
      
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_ancestry_diabetes', tests1[[x]],tests2[[x]]) # This should not change when you test subsetted data
        numoft1 <- length(which(meta2$sex_ancestry_diabetes==t1))
        numoft2 <- length(which(meta2$sex_ancestry_diabetes==t2))
        
        if (numoft1 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t1, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft2 < 3) {
          message(paste("!!WARNING!!"))
          message(paste(t2, "is <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste("####"))
          message(paste("####"))
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
    cell <- gsub('_sample_gex_total_counts.txt', '', FILE)
    raw_counts <- read.table(paste0(dir, FILE), row.names=1)
    sample_names <- unique(adata@meta.data$Library)
    sample_names <- gsub('-','.', sample_names)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% sample_names)]
    meta$Library2 <- gsub('-', '.', meta$Library)
    meta2 <- meta[which(meta$Library2 %in% colnames(raw_counts)),]
    
    if ('M' %in% meta2$Sex && 'F' %in% meta2$Sex){
      print(cell)
      print('Data for 2 sex present, however not all data may be present will check this at a later step')
      
      genes_to_keep <- c()
      for (i in 1:nrow(raw_counts)) {
        if (sum(raw_counts[i, ] >= 5) >= 2) {
          genes_to_keep <- c(genes_to_keep, rownames(raw_counts[i, ]))
        }
      }
      counts <- raw_counts[which(rownames(raw_counts) %in% genes_to_keep),] 
      
      if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && #Alldata 
          length(which(meta2$sex_diabetes == 'M_T2D')) > 1 && 
          length(which(meta2$sex_diabetes == 'F_ND')) > 1 && 
          length(which(meta2$sex_diabetes == 'F_T2D')) > 1) {
        my_design <- as.formula ('~Chemistry + Tissue_Source + sex_diabetes') # design for sex_diabetes
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else if (length(which(meta2$sex_diabetes == 'M_ND')) > 1 && 
                 length(which(meta2$sex_diabetes == 'F_ND')) > 1) {
        print(sprintf('%s does not have sufficient diabetes samples to test, bypassing to test ND only', cell))
        meta2 <- subset(meta2, Diabetes_Status == 'ND') # it is possible some T2D are present so eliminate them from your dataset since you are restricted to sex
        counts <- counts[,meta2$Library2]
        my_design <- as.formula ('~Tissue_Source + sex_diabetes') # design for sex_diabetes
        dds <- DESeqDataSetFromMatrix(counts, colData = meta2, design = my_design) #colData is where design columns are found
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds, maxit = 500) # https://support.bioconductor.org/p/65091/
      } else {print('samples not sufficient to test')}
      
      # No need for conditional formatting here
      # Specifying test combinations
      tests1 <- c('M_ND', 'M_T2D', 'F_T2D', 'M_T2D')
      tests2 <- c('F_ND', 'M_ND', 'F_ND', 'F_T2D')
      print('Preparing to run DESeq2')
      
      for (x in 1:length(tests1)){
        # No need for conditional formatting here
        t1 <- tests1[[x]]
        t2 <- tests2[[x]]
        test <- c('sex_diabetes', tests1[[x]],tests2[[x]]) # For Sex_diabetes
        numoft1 <- length(which(meta2$sex_diabetes==t1))
        numoft2 <- length(which(meta2$sex_diabetes==t2))
       
        if (numoft1 < 3) {
          message(paste("!!WARNING CHECK METADATA!!"))
          message(paste(t1, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste('####'))
          message(paste('####'))
        } else if (numoft2 < 3) {
          message(paste("!!WARNING CHECK METADATA!!"))
          message(paste(t2, "samples are <3 in the dataset, statistical threshold not met, analysis bypassed continuing with next iteration", sep= " "))
          message(paste("####"))
          message(paste("####"))
        } else if (numoft1 > 2 & numoft2 > 2) {
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
}) # System time

############################ STAGE ############################
############################   11  ############################
# Gene ontology analysis Rapid Gene ontology Auto Loader (Rapid GOAL)
# Create a list of all files in directory
system.time({
  dgelist <- list.files(r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\hpap)", 
                        all.files = FALSE, 
                        full.names = FALSE, 
                        pattern = "*.tsv")
  
  # Point towards WD using a function
  for (sample in dgelist){
    wd <- sprintf('C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/DETesting/DE_testing/hpap/%s', dgelist)
  }
  
  # Run iterative function to perform GO on all data
  for (x in wd) {
    sample_name <- str_split_fixed(x, "/", n=13)[13] # needs to be number of / in wd +1 (for alldata = 12 + 1)
    datfile <- read.table(file.path(x), sep = '\t', row.names = 1) 
    #datfile <- read.csv(file.path(x), row.names = 1)
    
    # Gene list of genes going UP
    sig_df_up <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange > 0.000000000014)
    sig_genes_up <- rownames(sig_df_up)
    
    # Gene list of genes going DOWN
    sig_df_down <- dplyr::filter(datfile, padj < 0.1 & log2FoldChange < -0.000000000014) 
    sig_genes_down <- rownames(sig_df_down)
    
    # All genes
    all_genes <- rownames(datfile)
    
    # Run GO enrichment analysis genes up
    GO.up <- enrichGO(gene = sig_genes_up, 
                      universe = all_genes, 
                      keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                      OrgDb = org.Hs.eg.db, 
                      ont = c("ALL"), 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, #if not set default is at 0.05
                      readable = TRUE)
    
    # Run GO enrichment analysis genes down
    GO.down <- enrichGO(gene = sig_genes_down, 
                        universe = all_genes, 
                        keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                        OrgDb = org.Hs.eg.db, 
                        ont = c("ALL"), 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 1, 
                        qvalueCutoff = 1, #if not set default is at 0.05
                        readable = TRUE)
    
    go_data_up <- data.frame(GO.up)
    go_data_down <- data.frame(GO.down)
    
    if (nrow(go_data_up) > 0) {
      go_data_up <- dplyr::filter(go_data_up, qvalue < 0.1) }
    if (nrow(go_data_down) > 0) {
      go_data_down <- dplyr::filter(go_data_down, qvalue < 0.1)}
    
    # Save outputs
    adjusted_name <- gsub('.{4}$', '', sample_name)
    adjusted_name <- gsub('deseq.WaldTest.', '', adjusted_name)
    if (nrow(go_data_up) > 0) {
    write.csv(go_data_up, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/hpap/UP/%s.csv", adjusted_name), row.names = FALSE)}
    if (nrow(go_data_down) > 0) {
    write.csv(go_data_down, file = sprintf("C:/Users/QadirMirzaMuhammadFa/Box/Lab 2301/1. R_Coding Scripts/Sex Biology Study/Data Output/scRNA/ORA/hpap/DOWN/%s.csv", adjusted_name), row.names = FALSE)}
  }
}) # Sys float time

############################ STAGE ############################
############################   12  ############################
# Load dataset
processed_rna <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
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

# Find markers to make lists of genes different across cells
Idents(processed_rna) <- "Diabetes Status"
nd.pancreas <- subset(processed_rna, idents = c("ND"))
#var.genes <- nd.pancreas@assays[["RNA"]]@var.features

# DE testing to determine celltype specificity
DefaultAssay(nd.pancreas) <- "SCT"
Idents(nd.pancreas) <- "ancestry_sex"
hispanic_M <- subset(nd.pancreas, idents = c("hispanic_M"))
hispanic_F <- subset(nd.pancreas, idents = c("hispanic_F"))
white_M <- subset(nd.pancreas, idents = c("white_M"))
white_F <- subset(nd.pancreas, idents = c("white_F"))
black_M <- subset(nd.pancreas, idents = c("black_M"))
black_F <- subset(nd.pancreas, idents = c("black_F"))

test_seurat <- black_F
Idents(test_seurat) <- "celltype_qadir"


plan(strategy = "multicore", workers = 80)
beta.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "beta", 
                                      latent.vars = "Library", 
                                      group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001, test.use = "MAST", 
                                      only.pos = TRUE)
beta.conserved.markers$p_val_adj[beta.conserved.markers$p_val_adj == 0] <- 2e-302
beta.conserved.markers <- dplyr::filter(beta.conserved.markers, p_val_adj < 1e-10)

alpha.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "alpha", 
                                       latent.vars = "Library", 
                                       group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                       test.use = "MAST", only.pos = TRUE)
alpha.conserved.markers$p_val_adj[alpha.conserved.markers$p_val_adj == 0] <- 2e-302
alpha.conserved.markers <- dplyr::filter(alpha.conserved.markers, p_val_adj < 1e-10) 

delta.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "delta", 
                                       latent.vars = "Library", 
                                       group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                       test.use = "MAST", only.pos = TRUE)
delta.conserved.markers$p_val_adj[delta.conserved.markers$p_val_adj == 0] <- 2e-302
delta.conserved.markers <- dplyr::filter(delta.conserved.markers, p_val_adj < 1e-10) 

gamma.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "gamma", 
                                       latent.vars = "Library", 
                                       group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                       test.use = "MAST", only.pos = TRUE)
gamma.conserved.markers$p_val_adj[gamma.conserved.markers$p_val_adj == 0] <- 2e-302
gamma.conserved.markers <- dplyr::filter(gamma.conserved.markers, p_val_adj < 1e-10) 

epsilon.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "epsilon", 
                                         latent.vars = "Library", 
                                         group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                         test.use = "MAST", only.pos = TRUE)
epsilon.conserved.markers$p_val_adj[epsilon.conserved.markers$p_val_adj == 0] <- 2e-302
epsilon.conserved.markers <- dplyr::filter(epsilon.conserved.markers, p_val_adj < 1e-10) 

betaalpha.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "beta+alpha", 
                                           latent.vars = "Library", 
                                           group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                           test.use = "MAST", only.pos = TRUE)
betaalpha.conserved.markers$p_val_adj[betaalpha.conserved.markers$p_val_adj == 0] <- 2e-302
betaalpha.conserved.markers <- dplyr::filter(betaalpha.conserved.markers, p_val_adj < 1e-10) 

betadelta.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "beta+delta", 
                                           latent.vars = "Library", 
                                           group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                           test.use = "MAST", only.pos = TRUE)
betadelta.conserved.markers$p_val_adj[betadelta.conserved.markers$p_val_adj == 0] <- 2e-302
betadelta.conserved.markers <- dplyr::filter(betadelta.conserved.markers, p_val_adj < 1e-10) 

cycling_endo.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "cycling_endo", 
                                              latent.vars = "Library", 
                                              group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                              test.use = "MAST", only.pos = TRUE)
cycling_endo.conserved.markers$p_val_adj[cycling_endo.conserved.markers$p_val_adj == 0] <- 2e-302
cycling_endo.conserved.markers <- dplyr::filter(cycling_endo.conserved.markers, p_val_adj < 1e-10) 

acinar.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "acinar", 
                                        latent.vars = "Library", 
                                        group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                        test.use = "MAST", only.pos = TRUE)
acinar.conserved.markers$p_val_adj[acinar.conserved.markers$p_val_adj == 0] <- 2e-302
acinar.conserved.markers <- dplyr::filter(acinar.conserved.markers, p_val_adj < 1e-10) 

ductal.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "ductal", 
                                        latent.vars = "Library", 
                                        group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                        test.use = "MAST", only.pos = TRUE)
ductal.conserved.markers$p_val_adj[ductal.conserved.markers$p_val_adj == 0] <- 2e-302
ductal.conserved.markers <- dplyr::filter(ductal.conserved.markers, p_val_adj < 1e-10) 

activated_stellate.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.40, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "activated_stellate", 
                                                    latent.vars = "Library", 
                                                    group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                                    test.use = "MAST", only.pos = TRUE)
activated_stellate.conserved.markers$p_val_adj[activated_stellate.conserved.markers$p_val_adj == 0] <- 2e-302
activated_stellate.conserved.markers <- dplyr::filter(activated_stellate.conserved.markers, p_val_adj < 1e-10) 

quiescent_stellate.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "quiescent_stellate", 
                                                    latent.vars = "Library", 
                                                    group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                                    test.use = "MAST", only.pos = TRUE)
quiescent_stellate.conserved.markers$p_val_adj[quiescent_stellate.conserved.markers$p_val_adj == 0] <- 2e-302
quiescent_stellate.conserved.markers <- dplyr::filter(quiescent_stellate.conserved.markers, p_val_adj < 1e-10) 

endothelial.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "endothelial", 
                                             latent.vars = "Library", 
                                             group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                             test.use = "MAST", only.pos = TRUE)
endothelial.conserved.markers$p_val_adj[endothelial.conserved.markers$p_val_adj == 0] <- 2e-302
endothelial.conserved.markers <- dplyr::filter(endothelial.conserved.markers, p_val_adj < 1e-10) 

lymphocyte.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "lymphocyte", 
                                            latent.vars = "Library", 
                                            group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                            test.use = "MAST", only.pos = TRUE)
lymphocyte.conserved.markers$p_val_adj[lymphocyte.conserved.markers$p_val_adj == 0] <- 2e-302
lymphocyte.conserved.markers <- dplyr::filter(lymphocyte.conserved.markers, p_val_adj < 1e-10) 

macrophages.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "macrophages", 
                                             latent.vars = "Library", 
                                             group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                             test.use = "MAST", only.pos = TRUE)
macrophages.conserved.markers$p_val_adj[macrophages.conserved.markers$p_val_adj == 0] <- 2e-302
macrophages.conserved.markers <- dplyr::filter(macrophages.conserved.markers, p_val_adj < 1e-10) 

mast.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "mast", 
                                      latent.vars = "Library", 
                                      group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                      test.use = "MAST", only.pos = TRUE)
mast.conserved.markers$p_val_adj[mast.conserved.markers$p_val_adj == 0] <- 2e-302
mast.conserved.markers <- dplyr::filter(mast.conserved.markers, p_val_adj < 1e-10) 

schwann.conserved.markers <- FindMarkers(test_seurat, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "data", ident.1 = "schwann", 
                                         latent.vars = "Library", 
                                         group.by = "celltype_qadir", min.cells.group = 1, pseudocount.use = 0.001,
                                         test.use = "MAST", only.pos = TRUE)
schwann.conserved.markers$p_val_adj[schwann.conserved.markers$p_val_adj == 0] <- 2e-302
schwann.conserved.markers <- dplyr::filter(schwann.conserved.markers, p_val_adj < 1e-10) 


#Save files:
write.csv(beta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\beta.csv)")
write.csv(alpha.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\alpha.csv)")
write.csv(delta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\delta.csv)")
write.csv(gamma.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\gamma.csv)")
write.csv(epsilon.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\epsilon.csv)")
write.csv(betaalpha.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\betaalpha.csv)")
write.csv(betadelta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\betadelta.csv)")
write.csv(cycling_endo.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\cycling_endo.csv)")
write.csv(acinar.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\acinar.csv)")
write.csv(ductal.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\ductal.csv)")
write.csv(activated_stellate.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\activated_stellate.csv)")
write.csv(quiescent_stellate.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\quiescent_stellate.csv)")
write.csv(endothelial.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\endothelial.csv)")
write.csv(lymphocyte.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\lymphocyte.csv)")
write.csv(macrophages.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\macrophages.csv)")
write.csv(mast.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\mast.csv)")
write.csv(schwann.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\black_F\schwann.csv)")

############################ STAGE ############################
############################   13  ############################
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

# Add pseudocount
combined_processed_rna[["RNA"]]@counts<-as.matrix(combined_processed_rna[["RNA"]]@counts)+1
{
plan(strategy = "multicore", workers = 80)
beta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "beta", 
                                      group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                      only.pos = TRUE)
beta.conserved.markers <- dplyr::filter(beta.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 1)
beta.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "beta", 
                                      latent.vars ="Library",
                                      group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                      only.pos = TRUE)
beta.conserved.markers.mast <- dplyr::filter(beta.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


alpha.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "alpha", 
                                       group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
alpha.conserved.markers <- dplyr::filter(alpha.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 1) 
alpha.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "alpha", 
                                           latent.vars ="Library",
                                           group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                           only.pos = TRUE)
alpha.conserved.markers.mast <- dplyr::filter(alpha.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


delta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "delta", 
                                       group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
delta.conserved.markers <- dplyr::filter(delta.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
delta.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                            latent.vars ="Library",
                                            group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                            only.pos = TRUE)
delta.conserved.markers.mast <- dplyr::filter(delta.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)



gamma.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "gamma", 
                                       group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                       only.pos = TRUE)
gamma.conserved.markers <- dplyr::filter(gamma.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
gamma.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                            latent.vars ="Library",
                                            group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                            only.pos = TRUE)
gamma.conserved.markers.mast <- dplyr::filter(gamma.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)



epsilon.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "epsilon", 
                                         group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                         only.pos = TRUE)
epsilon.conserved.markers <- dplyr::filter(epsilon.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
epsilon.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                            latent.vars ="Library",
                                            group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                            only.pos = TRUE)
epsilon.conserved.markers.mast <- dplyr::filter(epsilon.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


betaalpha.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "beta+alpha", 
                                           group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                           only.pos = TRUE)
betaalpha.conserved.markers <- dplyr::filter(betaalpha.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
betaalpha.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                              latent.vars ="Library",
                                              group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                              only.pos = TRUE)
betaalpha.conserved.markers.mast <- dplyr::filter(betaalpha.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


betadelta.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "beta+delta", 
                                           group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                           only.pos = TRUE)
betadelta.conserved.markers <- dplyr::filter(betadelta.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
betadelta.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                latent.vars ="Library",
                                                group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                only.pos = TRUE)
betadelta.conserved.markers.mast <- dplyr::filter(betadelta.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


cycling_endo.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "cycling-endo", 
                                              group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                              only.pos = TRUE)
cycling_endo.conserved.markers <- dplyr::filter(cycling_endo.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
cycling_endo.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                latent.vars ="Library",
                                                group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                only.pos = TRUE)
cycling_endo.conserved.markers.mast <- dplyr::filter(cycling_endo.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


acinar.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "acinar", 
                                        group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                        only.pos = TRUE)
acinar.conserved.markers <- dplyr::filter(acinar.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
acinar.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                   latent.vars ="Library",
                                                   group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                   only.pos = TRUE)
acinar.conserved.markers.mast <- dplyr::filter(acinar.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


ductal.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "ductal", 
                                        group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                        only.pos = TRUE)
ductal.conserved.markers <- dplyr::filter(ductal.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
ductal.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                             latent.vars ="Library",
                                             group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                             only.pos = TRUE)
ductal.conserved.markers.mast <- dplyr::filter(ductal.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


activated_stellate.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "activated-stellate", 
                                                    group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                                    only.pos = TRUE)
activated_stellate.conserved.markers <- dplyr::filter(activated_stellate.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
activated_stellate.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                             latent.vars ="Library",
                                             group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                             only.pos = TRUE)
activated_stellate.conserved.markers.mast <- dplyr::filter(activated_stellate.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


quiescent_stellate.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "quiescent-stellate", 
                                                    group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                                    only.pos = TRUE)
quiescent_stellate.conserved.markers <- dplyr::filter(quiescent_stellate.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
quiescent_stellate.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                         latent.vars ="Library",
                                                         group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                         only.pos = TRUE)
quiescent_stellate.conserved.markers.mast <- dplyr::filter(quiescent_stellate.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


endothelial.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "endothelial", 
                                             group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                             only.pos = TRUE)
endothelial.conserved.markers <- dplyr::filter(endothelial.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
endothelial.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                         latent.vars ="Library",
                                                         group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                         only.pos = TRUE)
endothelial.conserved.markers.mast <- dplyr::filter(endothelial.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


lymphocyte.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "lymphocyte", 
                                            group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                            only.pos = TRUE)
lymphocyte.conserved.markers <- dplyr::filter(lymphocyte.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
lymphocyte.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                  latent.vars ="Library",
                                                  group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                  only.pos = TRUE)
lymphocyte.conserved.markers.mast <- dplyr::filter(lymphocyte.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


macrophages.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "macrophages", 
                                             group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                             only.pos = TRUE)
macrophages.conserved.markers <- dplyr::filter(macrophages.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
macrophages.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                 latent.vars ="Library",
                                                 group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                 only.pos = TRUE)
macrophages.conserved.markers.mast <- dplyr::filter(macrophages.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


mast.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "mast", 
                                      group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                      only.pos = TRUE)
mast.conserved.markers <- dplyr::filter(mast.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
mast.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                                  latent.vars ="Library",
                                                  group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                                  only.pos = TRUE)
mast.conserved.markers.mast <- dplyr::filter(mast.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)


schwann.conserved.markers <- FindMarkers(combined_processed_rna, min.pct = 0.25, logfc.threshold = 1, assay = "RNA", slot = "counts", ident.1 = "schwann", 
                                         group.by = "celltype", min.cells.group = 1, test.use = "DESeq2", 
                                         only.pos = TRUE)
schwann.conserved.markers <- dplyr::filter(schwann.conserved.markers, p_val_adj < 1e-2 & avg_log2FC >= 0.001) 
schwann.conserved.markers.mast <- FindMarkers(nd.pancreas, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data", ident.1 = "delta", 
                                           latent.vars ="Library",
                                           group.by = "celltype_qadir", min.cells.group = 1, test.use = "MAST", 
                                           only.pos = TRUE)
schwann.conserved.markers.mast <- dplyr::filter(schwann.conserved.markers.mast, p_val_adj < 1e-2 & avg_log2FC >= 0.25)
}

#Save files:
write.csv(beta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\beta.csv)")
write.csv(alpha.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\alpha.csv)")
write.csv(delta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\delta.csv)")
write.csv(gamma.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\gamma.csv)")
write.csv(epsilon.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\epsilon.csv)")
write.csv(betaalpha.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\betaalpha.csv)")
write.csv(betadelta.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\betadelta.csv)")
write.csv(cycling_endo.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\cycling_endo.csv)")
write.csv(acinar.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\acinar.csv)")
write.csv(ductal.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\ductal.csv)")
write.csv(activated_stellate.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\activated_stellate.csv)")
write.csv(quiescent_stellate.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\quiescent_stellate.csv)")
write.csv(endothelial.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\endothelial.csv)")
write.csv(lymphocyte.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\lymphocyte.csv)")
write.csv(macrophages.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\macrophages.csv)")
write.csv(mast.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\mast.csv)")
write.csv(schwann.conserved.markers, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\schwann.csv)")

write.csv(beta.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\beta.csv)")
write.csv(alpha.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\alpha.csv)")
write.csv(delta.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\delta.csv)")
write.csv(gamma.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\gamma.csv)")
write.csv(epsilon.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\epsilon.csv)")
write.csv(betaalpha.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\betaalpha.csv)")
write.csv(betadelta.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\betadelta.csv)")
write.csv(cycling_endo.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\cycling_endo.csv)")
write.csv(acinar.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\acinar.csv)")
write.csv(ductal.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\ductal.csv)")
write.csv(activated_stellate.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\activated_stellate.csv)")
write.csv(quiescent_stellate.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\quiescent_stellate.csv)")
write.csv(endothelial.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\endothelial.csv)")
write.csv(lymphocyte.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\lymphocyte.csv)")
write.csv(macrophages.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\macrophages.csv)")
write.csv(mast.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\mast.csv)")
write.csv(schwann.conserved.markers.mast, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\mast\schwann.csv)")











