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
# Macs2 component run in linux
library(future)
plan()

plan("multicore", workers = 14)
plan()
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

hm.integrated.dfree$celltype_sex_ancestry <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, sep = "_")
table(hm.integrated.dfree$celltype_sex_ancestry)
sum(table(hm.integrated.dfree$celltype_sex_ancestry))
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
# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree@meta.data$celltype_sex_ancestry <- factor(x = hm.integrated.dfree@meta.data$celltype_sex_ancestry, levels = my_levels2)
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry"
table(hm.integrated.dfree$celltype_sex_ancestry)
sum(table(hm.integrated.dfree$celltype_sex_ancestry))
sum(table(hm.integrated.dfree$celltype_sex_ancestry_lib))

# Correct fragmnet paths
frags <- Fragments(hm.integrated.dfree)
Fragments(hm.integrated.dfree) <- NULL  # remove fragment information from assay

# We need to find the correct order of all paths
# First store old paths so we can switch back later
windows_paths <- c(frags[[1]]@path, frags[[2]]@path, frags[[3]]@path, frags[[4]]@path, frags[[5]]@path, 
                   frags[[6]]@path, frags[[7]]@path, frags[[8]]@path, frags[[9]]@path, frags[[10]]@path,
                   frags[[11]]@path, frags[[12]]@path, frags[[13]]@path, frags[[14]]@path, frags[[15]]@path)
windows_paths

# Store paths for linux
linux_paths <- c("//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//1_220628_Fahd_snATAC1_HP-20228-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//2_220701_Fahd_snATAC2_SAMN15877725//fragments.tsv.gz", 
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//3_220701_Fahd_snATAC3_HP-20240-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//4_220630_Fahd_snATAC4_HP-20314-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//5_220303_snATAC_F52_HP-21055-01//fragments.tsv.gz",    
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//6_210401_snATAC_F62_HP-21062-01//fragments.tsv.gz",    
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//7_210401_snATAC_F7a_HP-21070-01//fragments.tsv.gz",    
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//9_210628_snATAC_F9a_HP-21079-01//fragments.tsv.gz",    
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//10_210628_snATAC_F10a_HP-21086-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//11_210714_snATAC_F11a_HP-21089-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//12_210714_snATAC_F12a_HP-21100-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//13_211208_snATAC_F13_HP-21216-01//fragments.tsv.gz",   
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//14_211208_snATAC_F14_HP-21232-01//fragments.tsv.gz",   
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//15_220303_snATAC_F15a_HP-21328-01//fragments.tsv.gz",  
                 "//media//mirza//DATA2//1.SexbasedStudyrawdata//Cellranger_raw_data//snATACseq//16_220630_Fahd_snATAC16_HP-22021-01//fragments.tsv.gz")

# create a vector with all the new paths, in the correct order for your list of fragment objects
# Switch between paths, windows or linux
new.paths <- linux_paths
for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(hm.integrated.dfree) <- frags # assign updated list back to the object
Fragments(hm.integrated.dfree)


peaks <- CallPeaks(
  object = hm.integrated.dfree,
  group.by = "celltype_sex_ancestry"
)

# Filter peaks
peakwidths <- width(peaks)
peaks <- peaks[peakwidths  < 10000 & peakwidths > 20]
peaks

Idents(hm.integrated.dfree) <- "celltype"
CoveragePlot(
  object = hm.integrated.dfree,
  region = "INS",
  ranges = peaks,
  ranges.title = "MACS2"
)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(hm.integrated.dfree),
  features = peaks,
  cells = colnames(hm.integrated.dfree)
)

hm.integrated.dfree[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(hm.integrated.dfree),
annotation = Annotation(hm.integrated.dfree)
)

qsave(hm.integrated.dfree, "//media//mirza//DATA2//2.SexbasedStudyCurrent//QS files//hm.integrated.dfree.macs2.qs")
qsave(peaks, "//media//mirza//DATA2//2.SexbasedStudyCurrent//QS files//macs2_peaks.qs")
qsave(macs2_counts, "//media//mirza//DATA2//2.SexbasedStudyCurrent//QS files//macs2_counts.qs")
qsave(windows_paths, "//media//mirza//DATA2//2.SexbasedStudyCurrent//QS files//windows_paths.qs")
qsave(linux_paths, "//media//mirza//DATA2//2.SexbasedStudyCurrent//QS files//linux_paths.qs")

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


> sessionInfo()
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DESeq2_1.36.0                     ggseqlogo_0.1                     BiocParallel_1.32.0               cicero_1.3.9                     
 [5] Gviz_1.42.1                       SeuratWrappers_0.3.1              chromVAR_1.5.0                    motifmatchr_1.18.0               
 [9] biovizBase_1.44.0                 BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.67.1                   rtracklayer_1.56.1               
[13] Biostrings_2.66.0                 XVector_0.38.0                    TFBSTools_1.34.0                  JASPAR2020_0.99.10               
[17] qs_0.25.5                         R.utils_2.12.2                    R.oo_1.25.0                       R.methodsS3_1.8.2                
[21] devtools_2.4.5                    usethis_2.1.6                     ggVennDiagram_1.2.2               ggvenn_0.1.10                    
[25] DropletUtils_1.16.0               Nebulosa_1.6.0                    scCustomize_1.1.3                 circlize_0.4.15                  
[29] ComplexHeatmap_2.12.1             viridis_0.6.4                     viridisLite_0.4.1                 EnrichmentBrowser_2.26.0         
[33] graph_1.74.0                      escape_1.6.0                      dittoSeq_1.8.1                    DOSE_3.22.1                      
[37] clusterProfiler_4.4.4             MeSHDbi_1.32.0                    AnnotationHub_3.4.0               BiocFileCache_2.4.0              
[41] dbplyr_2.3.3                      org.Hs.eg.db_3.15.0               GOSemSim_2.22.0                   glmGamPoi_1.8.0                  
[45] EnhancedVolcano_1.14.0            DoubletFinder_2.0.3               future_1.29.0                     patchwork_1.1.2                  
[49] clustree_0.5.0                    ggraph_2.1.0                      plotly_4.10.1                     EnsDb.Hsapiens.v86_2.99.0        
[53] ensembldb_2.20.2                  AnnotationFilter_1.20.0           GenomicFeatures_1.48.4            AnnotationDbi_1.58.0             
[57] scDblFinder_1.10.0                Signac_1.8.0                      harmony_0.1.1                     monocle3_1.2.9                   
[61] SingleCellExperiment_1.20.0       SummarizedExperiment_1.28.0       GenomicRanges_1.50.1              GenomeInfoDb_1.32.4              
[65] IRanges_2.32.0                    S4Vectors_0.36.0                  MatrixGenerics_1.10.0             matrixStats_0.62.0               
[69] Biobase_2.58.0                    BiocGenerics_0.44.0               SeuratObject_4.1.3                Seurat_4.3.0.1                   
[73] reticulate_1.26                   data.table_1.14.4                 lubridate_1.9.2                   forcats_1.0.0                    
[77] purrr_1.0.1                       readr_2.1.4                       tidyr_1.3.0                       tibble_3.2.1                     
[81] tidyverse_2.0.0                   dplyr_1.1.2                       ggridges_0.5.4                    Matrix_1.5-1                     
[85] cowplot_1.1.1                     Rcpp_1.0.9                        SoupX_1.6.2                       hdf5r_1.3.8                      
[89] stringr_1.5.0                     leiden_0.4.3                      ggrepel_0.9.2                     ggplot2_3.4.2                    

loaded via a namespace (and not attached):
  [1] terra_1.6-17                  graphlayouts_0.8.3            pbapply_1.5-0                 lattice_0.20-45               GSVA_1.44.5                  
  [6] vctrs_0.6.3                   blob_1.2.3                    survival_3.4-0                spatstat.data_3.0-1           later_1.3.0                  
 [11] nloptr_2.0.3                  DBI_1.1.3                     rappdirs_0.3.3                uwot_0.1.14                   dqrng_0.3.0                  
 [16] jpeg_0.1-10                   zlibbioc_1.44.0               rgeos_0.5-9                   htmlwidgets_1.5.4             mvtnorm_1.1-3                
 [21] GlobalOptions_0.1.2           parallel_4.2.2                scater_1.24.0                 irlba_2.3.5.1                 tidygraph_1.2.2              
 [26] KernSmooth_2.23-20            DT_0.26                       promises_1.2.0.1              DelayedArray_0.24.0           limma_3.54.0                 
 [31] pkgload_1.3.1                 Hmisc_5.1-0                   RcppParallel_5.1.7            fs_1.5.2                      fastmatch_1.1-3              
 [36] digest_0.6.30                 png_0.1-7                     bluster_1.6.0                 sctransform_0.3.5             scatterpie_0.1.8             
 [41] here_1.0.1                    pkgconfig_2.0.3               GO.db_3.15.0                  spatstat.random_3.1-5         DelayedMatrixStats_1.20.0    
 [46] ggbeeswarm_0.6.0              iterators_1.0.14              minqa_1.2.5                   beeswarm_0.4.0                xfun_0.39                    
 [51] GetoptLong_1.0.5              zoo_1.8-11                    tidyselect_1.2.0              reshape2_1.4.4                ica_1.0-3                    
 [56] pkgbuild_1.3.1                rlang_1.1.1                   RVenn_1.1.0                   glue_1.6.2                    RColorBrewer_1.1-3           
 [61] CNEr_1.32.0                   httpuv_1.6.6                  BiocNeighbors_1.16.0          seqLogo_1.62.0                DO.db_2.9                    
 [66] annotate_1.74.0               jsonlite_1.8.7                bit_4.0.4                     mime_0.12                     gridExtra_2.3                
 [71] Rsamtools_2.14.0              stringi_1.7.8                 processx_3.8.0                RcppRoll_0.3.0                spatstat.sparse_3.0-2        
 [76] scattermore_0.8               spatstat.explore_3.2-1        yulab.utils_0.0.5             bitops_1.0-7                  cli_3.6.1                    
 [81] rhdf5filters_1.10.0           RSQLite_2.2.18                pheatmap_1.0.12               KEGGgraph_1.56.0              timechange_0.2.0             
 [86] rstudioapi_0.14               GenomicAlignments_1.32.1      nlme_3.1-160                  qvalue_2.28.0                 scran_1.24.1                 
 [91] ggprism_1.0.4                 locfit_1.5-9.8                janitor_2.2.0                 ks_1.14.0                     VariantAnnotation_1.42.1     
 [96] listenv_0.8.0                 miniUI_0.1.1.1                gridGraphics_0.5-1            urlchecker_1.0.1              sessioninfo_1.2.2            
[101] lifecycle_1.0.3               munsell_0.5.0                 caTools_1.18.2                codetools_0.2-18              vipor_0.4.5                  
[106] lmtest_0.9-40                 msigdbr_7.5.1                 htmlTable_2.4.1               xtable_1.8-4                  ROCR_1.0-11                  
[111] BiocManager_1.30.19           abind_1.4-5                   farver_2.1.1                  parallelly_1.32.1             RANN_2.6.1                   
[116] aplot_0.1.8                   poweRlaw_0.70.6               ggtree_3.7.1                  BiocIO_1.6.0                  RcppAnnoy_0.0.20             
[121] goftest_1.2-3                 dichromat_2.0-0.1             profvis_0.3.7                 cluster_2.1.4                 future.apply_1.10.0          
[126] tidytree_0.4.1                ellipsis_0.3.2                prettyunits_1.1.1             mclust_6.0.0                  igraph_1.3.5                 
[131] fgsea_1.22.0                  remotes_2.4.2                 paletteer_1.5.0               spatstat.utils_3.0-3          htmltools_0.5.3              
[136] yaml_2.3.6                    utf8_1.2.2                    interactiveDisplayBase_1.34.0 XML_3.99-0.12                 foreign_0.8-83               
[141] withr_2.5.0                   scuttle_1.8.0                 fitdistrplus_1.1-8            bit64_4.0.5                   xgboost_1.7.5.1              
[146] foreach_1.5.2                 ProtGenerics_1.28.0           progressr_0.11.0              rsvd_1.0.5                    ScaledMatrix_1.6.0           
[151] VGAM_1.1-7                    evaluate_0.18                 memoise_2.0.1                 RApiSerialize_0.1.2           geneplotter_1.74.0           
[156] tzdb_0.4.0                    callr_3.7.3                   ps_1.7.2                      curl_4.3.3                    fansi_1.0.3                  
[161] GSEABase_1.58.0               tensor_1.5                    edgeR_3.38.4                  checkmate_2.1.0               cachem_1.0.6                 
[166] interp_1.1-4                  deldir_1.0-6                  babelgene_22.9                metapod_1.4.0                 rjson_0.2.21                 
[171] clue_0.3-64                   rprojroot_2.0.3               tools_4.2.2                   magrittr_2.0.3                RCurl_1.98-1.9               
[176] TFMPvalue_0.0.9               ape_5.6-2                     ggplotify_0.1.0               xml2_1.3.3                    rmarkdown_2.17               
[181] httr_1.4.6                    boot_1.3-28                   globals_0.16.1                R6_2.5.1                      Rhdf5lib_1.20.0              
[186] nnet_7.3-18                   genefilter_1.78.0             DirichletMultinomial_1.38.0   progress_1.2.2                KEGGREST_1.36.3              
[191] treeio_1.23.0                 gtools_3.9.3                  shape_1.4.6                   statmod_1.5.0                 beachmat_2.14.0              
[196] BiocVersion_3.15.2            rematch2_2.1.2                HDF5Array_1.26.0              BiocSingular_1.14.0           ggrastr_1.0.1                
[201] rhdf5_2.42.0                  splines_4.2.2                 snakecase_0.11.0              ggfun_0.0.8                   colorspace_2.0-3             
[206] generics_0.1.3                base64enc_0.1-3               pracma_2.4.2                  pillar_1.9.0                  Rgraphviz_2.40.0             
[211] tweenr_2.0.2                  sp_1.5-1                      GenomeInfoDbData_1.2.9        plyr_1.8.7                    gtable_0.3.1                 
[216] stringfish_0.15.8             restfulr_0.0.15               latticeExtra_0.6-30           knitr_1.40                    shadowtext_0.1.2             
[221] biomaRt_2.52.0                fastmap_1.1.0                 doParallel_1.0.17             broom_1.0.5                   UCell_2.0.1                  
[226] scales_1.2.1                  filelock_1.0.2                backports_1.4.1               lme4_1.1-31                   enrichplot_1.16.2            
[231] hms_1.1.2                     ggforce_0.4.1                 Rtsne_0.16                    shiny_1.7.3                   polyclip_1.10-4              
[236] lazyeval_0.2.2                Formula_1.2-5                 crayon_1.5.2                  MASS_7.3-58.1                 downloader_0.4               
[241] sparseMatrixStats_1.10.0      rpart_4.1.19                  compiler_4.2.2                spatstat.geom_3.2-4  
