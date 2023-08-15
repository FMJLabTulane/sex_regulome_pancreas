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

packageVersion("monocle3")
packageVersion("monocle3")

############################ STAGE ############################
############################   12  ############################
# Analysis using macs2 counts

############################ STAGE ############################
############################   13  ############################
# Load package
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
DefaultAssay(hm.integrated.dfree) <- "macs2"

# Switch frag paths
frags <- Fragments(hm.integrated.dfree)
Fragments(hm.integrated.dfree) <- NULL  # remove fragment information from assay

# We need to find the correct order of all paths
# First store old paths so we can switch back later
linux_paths <- c(frags[[1]]@path, frags[[2]]@path, frags[[3]]@path, frags[[4]]@path, frags[[5]]@path, 
                 frags[[6]]@path, frags[[7]]@path, frags[[8]]@path, frags[[9]]@path, frags[[10]]@path,
                 frags[[11]]@path, frags[[12]]@path, frags[[13]]@path, frags[[14]]@path, frags[[15]]@path)
linux_paths

windows_paths <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\windows_paths.qs)")
windows_paths

# # Store paths for linux
# linux_paths <- c("//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//1_220628_Fahd_snATAC1_HP-20228-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//2_220701_Fahd_snATAC2_SAMN15877725//fragments.tsv.gz", 
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//3_220701_Fahd_snATAC3_HP-20240-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//4_220630_Fahd_snATAC4_HP-20314-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//5_220303_snATAC_F52_HP-21055-01//fragments.tsv.gz",    
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//6_210401_snATAC_F62_HP-21062-01//fragments.tsv.gz",    
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//7_210401_snATAC_F7a_HP-21070-01//fragments.tsv.gz",    
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//9_210628_snATAC_F9a_HP-21079-01//fragments.tsv.gz",    
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//10_210628_snATAC_F10a_HP-21086-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//11_210714_snATAC_F11a_HP-21089-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//12_210714_snATAC_F12a_HP-21100-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//13_211208_snATAC_F13_HP-21216-01//fragments.tsv.gz",   
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//14_211208_snATAC_F14_HP-21232-01//fragments.tsv.gz",   
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//15_220303_snATAC_F15a_HP-21328-01//fragments.tsv.gz",  
#                  "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//16_220630_Fahd_snATAC16_HP-22021-01//fragments.tsv.gz")

# create a vector with all the new paths, in the correct order for your list of fragment objects
# Switch between paths, windows or linux
new.paths <- windows_paths
for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(hm.integrated.dfree) <- frags # assign updated list back to the object
Fragments(hm.integrated.dfree)

# Create RNA slot based on macs2
gene.activities.macs2 <- GeneActivity(hm.integrated.dfree)
hm.integrated.dfree[['RNA_macs2']] <- CreateAssayObject(counts = gene.activities.macs2)

# Normalize
DefaultAssay(hm.integrated.dfree) <- "RNA_macs2"
hm.integrated.dfree <- NormalizeData(
  object = hm.integrated.dfree,
  assay = 'RNA_macs2',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated.dfree$nCount_RNA)
)

############################ STAGE ############################
############################  14   ############################
# NORMALIZATION AND BATCH CORRECTION
# Run TFDIF  
DefaultAssay(hm.integrated.dfree) <- "macs2"
hm.integrated.dfree <- RunTFIDF(hm.integrated.dfree)
hm.integrated.dfree <- FindTopFeatures(hm.integrated.dfree, min.cutoff = 20)

# View UMAP
DimPlot(hm.integrated.dfree, group.by = "predicted.id", reduction = "umap", label = TRUE) + ggtitle("Predicted annotation")
DimPlot(hm.integrated.dfree, group.by = "ATAC_snn_res.0.8", reduction = "umap", label = TRUE) + ggtitle("Res")# + nolegend()

# Observing cells
Idents(hm.integrated.dfree) <- "celltype"
DimPlot(hm.integrated.dfree, 
        #split.by = "ancestry_sex", 
        #group.by = "celltype", 
        label = TRUE, 
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

#Peak calling
CoveragePlot(
  object = hm.integrated.dfree,
  assay = "macs2",
  group.by = "celltype",
  region = c("INS"
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
hm.integrated.dfree[['macs2']]
granges(hm.integrated.dfree)
table(hm.integrated.dfree@assays[["macs2"]]@annotation@seqinfo@genome)

############################ STAGE ############################
############################  16   ############################
# MOTIF INFORMATION ADDITION
# Adding Motifs
# Get a list of motif position frequency matrices from the JASPAR database
# Change back to peak data
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
hm.integrated.dfree <- AddMotifs(
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "macs2" # note appropriate assay
)

#Compute motif activity
# THIS STEP IS CRAZYYYYYYYYYYYY IT TAKES AGES OR BSOD's PC IF PARALLEL IS NOT TURNED OFF
register(SerialParam()) # VERY IMPORTANT TRANSITION FROM PARALLEL TO SERIAL COMPUTING
hm.integrated.dfree <- RunChromVAR(
  new.assay.name = "chromvar_macs2",
  object = hm.integrated.dfree,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "macs2" # note appropriate assay
)

DefaultAssay(hm.integrated.dfree) <- 'chromvar_macs2'

# SAVE YOUR FILE NOW
qsave(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.qs)")
hm.integrated.dfree <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.qs)")

# look at the activity of PDX1
FeaturePlot(
  object = hm.integrated.dfree,
  features = "MA0132.2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

# look at the activity of Arx
FeaturePlot(
  object = hm.integrated.dfree,
  features = "MA0874.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

############################ STAGE ############################
############################  17   ############################
# convert to CellDataSet format and make the cicero object
DefaultAssay(hm.integrated.dfree) <- 'macs2'
Idents(hm.integrated.dfree) <- "celltype"
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


############################ STAGE ############################
############################  17   ############################
table(hm.integrated.dfree@meta.data[["celltype"]])
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype"

# DE testing
#plan()
#plan("multisession", workers = 20)
#options(future.globals.maxSize = 20 * 1024^3)

{
beta_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "beta", 
  ident.2 = c("alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-308
write.csv(beta_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\beta_peaks.csv)")


alpha_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "alpha", 
  ident.2 = c("beta", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-308
write.csv(alpha_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\alpha_peaks.csv)")

delta_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "delta", 
  ident.2 = c("beta", "alpha", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-308
write.csv(delta_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\delta_peaks.csv)")

gamma_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "gamma", 
  ident.2 = c("beta", "alpha", "delta", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-308
write.csv(gamma_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\gamma_peaks.csv)")

acinar_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "acinar", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-308
write.csv(acinar_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\acinar_peaks.csv)")

ductal_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "ductal", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-308
write.csv(ductal_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\ductal_peaks.csv)")

activatedstellate_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "activated_stellate", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
activatedstellate_peaks$p_val_adj[activatedstellate_peaks$p_val_adj == 0] <- 2e-308
write.csv(activatedstellate_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\activatedstellate_peaks.csv)")

quiescentstellate_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "quiescent_stellate", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "endothelial", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
quiescentstellate_peaks$p_val_adj[quiescentstellate_peaks$p_val_adj == 0] <- 2e-308
write.csv(quiescentstellate_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\quiescentstellate_peaks.csv)")

endothelial_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "endothelial", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophage", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
endothelial_peaks$p_val_adj[endothelial_peaks$p_val_adj == 0] <- 2e-308
write.csv(endothelial_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\endothelial_peaks.csv)")

macrophage_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "macrophage", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "lymphocyte"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
macrophage_peaks$p_val_adj[macrophage_peaks$p_val_adj == 0] <- 2e-308
write.csv(macrophage_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\macrophage_peaks.csv)")

lymphocyte_peaks <- FindMarkers(
  object = hm.integrated.dfree,
  logfc.threshold = 0.263034, #1.2
  densify = TRUE,
  only.pos = TRUE,
  ident.1 = "lymphocyte", 
  ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_RNA_macs2'
)
lymphocyte_peaks$p_val_adj[lymphocyte_peaks$p_val_adj == 0] <- 2e-308
write.csv(lymphocyte_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\lymphocyte_peaks.csv)")
}


## TEST ACROSS SEX #
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype_sex"
sum(table(hm.integrated.dfree@meta.data[["celltype_sex"]]))
table(hm.integrated.dfree@meta.data[["celltype_sex"]])

# DE testing
#plan()
#plan("multisession", workers = 20)
#options(future.globals.maxSize = 20 * 1024^3)

{
  beta_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "beta_male", 
    ident.2 = c("beta_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  beta_mvf_peaks$p_val_adj[beta_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(beta_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\beta_mvf_peaks.csv)")
  
  
  alpha_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "alpha_male", 
    ident.2 = c("alpha_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  alpha_mvf_peaks$p_val_adj[alpha_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(alpha_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\alpha_mvf_peaks.csv)")
  
  delta_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "delta_male", 
    ident.2 = c("delta_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  delta_mvf_peaks$p_val_adj[delta_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(delta_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\delta_mvf_peaks.csv)")
  
  gamma_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "gamma_male", 
    ident.2 = c("gamma_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  gamma_mvf_peaks$p_val_adj[gamma_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(gamma_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\gamma_mvf_peaks.csv)")
  
  acinar_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "acinar_male", 
    ident.2 = c("acinar_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  acinar_mvf_peaks$p_val_adj[acinar_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(acinar_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\acinar_mvf_peaks.csv)")
  
  ductal_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "ductal_male", 
    ident.2 = c("ductal_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  ductal_mvf_peaks$p_val_adj[ductal_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(ductal_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\ductal_mvf_peaks.csv)")
  
  activatedstellate_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "activated_stellate_male", 
    ident.2 = c("activated_stellate_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  activatedstellate_mvf_peaks$p_val_adj[activatedstellate_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(activatedstellate_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\activatedstellate_mvf_peaks.csv)")
  
  quiescentstellate_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "quiescent_stellate_male", 
    ident.2 = c("quiescent_stellate_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  quiescentstellate_mvf_peaks$p_val_adj[quiescentstellate_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(quiescentstellate_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\quiescentstellate_mvf_peaks.csv)")
  
  endothelial_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "endothelial_male", 
    ident.2 = c("endothelial_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  endothelial_mvf_peaks$p_val_adj[endothelial_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(endothelial_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\endothelial_mvf_peaks.csv)")
  
  macrophage_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "macrophage_male", 
    ident.2 = c("macrophage_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  macrophage_mvf_peaks$p_val_adj[macrophage_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(macrophage_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\macrophage_mvf_peaks.csv)")
  
  lymphocyte_mvf_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034, #1.2
    densify = TRUE,
    only.pos = FALSE,
    ident.1 = "lymphocyte_male", 
    ident.2 = c("lymphocyte_female"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  lymphocyte_mvf_peaks$p_val_adj[lymphocyte_mvf_peaks$p_val_adj == 0] <- 2e-308
  write.csv(lymphocyte_mvf_peaks, file = r"(C:\Users\QadirMirzaMuhammadFa\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\lymphocyte_mvf_peaks.csv)")
}



#### END ####

















CoveragePlot(hm.integrated.dfree, region = c("CFTR", "PRSS2", "ESM1", "SDS"), 
                    window = 100,
                    ymax = 100,
                    links = TRUE,
                    #tile = TRUE,
                    extend.upstream = 50000,
                    extend.downstream = 50000,
             ncol = 4) & scale_fill_manual(values = c("dodgerblue3",      #beta
                                                      "chartreuse3",      #delta
                                                      "lightseagreen",    #alpha
                                                      "springgreen4",     #gamma
                                                      "salmon3",          #acinar
                                                      "darkorange2",      #ductal
                                                      "salmon",           #quiescent_stellate
                                                      "orange",           #activated_setallate
                                                      "red",              #endothelial
                                                      "orchid1",          #lymphocytes
                                                      "magenta3"))
