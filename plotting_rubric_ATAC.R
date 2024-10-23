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
install.packages("rlang")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

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

BiocManager::install("clusterProfiler")

# Pando compatible
remotes::install_version(package = 'Seurat', version = package_version('4.3.0.1'))
devtools::install_version("SeuratObject", version = "4.1.4", repos = "http://cran.us.r-project.org")
devtools::install_version("Signac", version = "1.11.0", repos = "http://cran.us.r-project.org")

BiocManager::install("scDblFinder", dependencies = T, force = TRUE)
BiocManager::install("HDO.db")

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
install.packages("ggseqlogo", force = TRUE)
install.packages('fastmap', force = TRUE)
devtools::install_github('quadbiolab/Pando')
usethis::create_github_token()

#devtools::install_github("stuart-lab/signac", ref = "develop", force = TRUE)
BiocManager::install(c("GenomicRanges", force = TRUE))
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

# Run the following code once you have Seurat installed
suppressPackageStartupMessages({
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
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(motifmatchr)
  library(chromVAR)
  #library(SeuratWrappers)
  library(cicero)
  library(BiocParallel)
  #library(ggseqlogo)
  library(DESeq2)
  library(Pando)
  # Set global environment parameter par-proc
  # options(future.globals.maxSize = 8000 * 1024^2)
  set.seed(1234)
})
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
packageVersion("scDblFinder")

# Load object
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")
#combined_atac <- qread(file = r"(E:\2.SexbasedStudyCurrent\QS files\snATACfiles_earlierpartsofworkflow\combined_atac.qs)")
saveRDS(hm.integrated.dfree, file = r"(E:\2.SexbasedStudyCurrent\RDS files\hm.integrated.dfree.rds)")

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma", 
               "acinar", "ductal",
               "quiescent_stellate", "activated_stellate", "endothelial",
               "lymphocyte", "macrophage")

table(hm.integrated.dfree$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
hm.integrated.dfree$celltype <- factor(x = hm.integrated.dfree$celltype, levels = my_levels)
table(unique(hm.integrated.dfree$celltype))

#QC
frac_reads_in_promoters <- combined_atac@meta.data$FRiP
quantile(frac_reads_in_promoters)
metadata <- combined_atac@meta.data
ggplot(data=metadata, mapping = aes(x=frac_reads_in_promoters)) +  geom_density(alpha = 0.2, fill= 'lightpink', color="pink") + 
  theme_linedraw() + geom_vline(xintercept=c(0.22,0.2,0.38), colour=c("blue", "red", "purple"),linetype = "longdash")

frac_reads_in_peaks <- combined_atac@meta.data$pct_reads_in_peaks
quantile(frac_reads_in_peaks)
ggplot(data=metadata, mapping = aes(x=frac_reads_in_peaks)) +  geom_density(alpha = 0.2, color="green", fill="lightgreen") + 
  theme_linedraw() + geom_vline(xintercept=c(30), colour=c("red"),linetype = "longdash")

TSS <- combined_atac@meta.data$TSS.enrichment
quantile(TSS)
ggplot(data=metadata, mapping = aes(x=TSS)) +  geom_density(alpha = 0.2, color="green", fill="lightgreen") + 
  theme_linedraw() + geom_vline(xintercept=c(2), colour=c("red"),linetype = "longdash")


DensityScatter(hm.integrated.dfree, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

Idents(hm.integrated.dfree) <- "sex"
Idents(combined_atac) <- "sex"
VlnPlot(
  object = hm.integrated.dfree,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', "FRiP"),
  pt.size = 0,
  ncol = 3
)

VlnPlot(
  object = combined_atac,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', "FRiP"),
  pt.size = 0,
  ncol = 3
)

# Dimplot
unique(hm.integrated.dfree@meta.data[["celltype"]])
DimPlot(hm.integrated.dfree, #switch here to plot
        #split.by = "Diabetes Status", 
        group.by = "celltype", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.05,
        cols = c("dodgerblue3",      #beta
                 "chartreuse3",      #delta
                 "lightseagreen",    #alpha
                 "springgreen4",     #gamma
                 "salmon3",          #acinar
                 "darkorange2",      #ductal
                 "salmon",           #quiescent_stellate
                 "orange",           #activated_setallate
                 "red",              #endothelial
                 "orchid1",          #lymphocytes
                 "magenta3"         #macrophages
        )
)

#
hm.integrated.dfree$sex_ancestry_sample <- paste(hm.integrated.dfree$sex, hm.integrated.dfree$ancestry, hm.integrated.dfree$sample, sep = "_")
dittoBarPlot(hm.integrated.dfree, "celltype", 
             retain.factor.levels = TRUE,
             scale = "percent",
             color.panel = c("dodgerblue3",      #beta
                             "chartreuse3",      #delta
                             "lightseagreen",    #alpha
                             "springgreen4",     #gamma
                             "salmon3",          #acinar
                             "darkorange2",      #ductal
                             "salmon",           #quiescent_stellate
                             "orange",           #activated_setallate
                             "red",              #endothelial
                             "orchid1",          #lymphocytes
                             "magenta3"),          
             group.by = "sex_ancestry_sample") + coord_flip()


# Umap of sex
DimPlot(hm.integrated.dfree, 
        #split.by = "Diabetes Status", http://127.0.0.1:42565/graphics/plot_zoom_png?width=1160&height=900
        group.by = "sex", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("red4", #ND
                 "deepskyblue3"         #T2D  
        ))
# Umap of ancestry_sex
DimPlot(hm.integrated.dfree, 
        #split.by = "Diabetes Status", 
        group.by = "ancestry", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.1,
        shuffle = TRUE,
        cols = c("black", #ND
                 "deepskyblue3"
        ))

# Umap of Accessibility
Idents(hm.integrated.dfree) <- "celltype"
DefaultAssay(hm.integrated.dfree) <- "RNA_macs2"

FeaturePlot(
  object = hm.integrated.dfree,
  features = c("GCG"),
  cols = c("lightgrey", "red4"),
  pt.size = 0.1,
  max.cutoff = 2,
  ncol = 1
  #raster = TRUE,
  #order = TRUE
)

# #, "GCG", "SST", "PPY", "CFTR", "PRSS1", "ESM1", "SDS", "RGS5", "PDGFRA", "CD7"
# 
# # Corelation plot
# # # Pseudobulk
# # Idents(hm.integrated.dfree) <- "celltype"
# # hm.integrated.dfree$celltype_sex <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, sep = "_")
# # table(hm.integrated.dfree@meta.data$celltype_sex)
# # Idents(hm.integrated.dfree) <- "celltype_sex"
# # DefaultAssay(hm.integrated.dfree) <- "ATAC"
# # combined_processed_atac <- Seurat:::PseudobulkExpression(object = hm.integrated.dfree, 
# #                                                          pb.method = 'aggregate', 
# #                                                          return.seurat = TRUE,
# #                                                          slot = 'counts')
# # 
# # DefaultAssay(combined_processed_atac) <- "ATAC"
# # combined_processed_atac <- RunTFIDF(combined_processed_atac, assay = "ATAC")
# # 
# # {
# #   combined_processed_atac$celltype_sex <- combined_processed_atac@active.ident
# #   Idents(combined_processed_atac) <- 'celltype_sex'
# #   combined_processed_atac$celltype <- combined_processed_atac$orig.ident
# #   metadat <- combined_processed_atac@meta.data
# #   metadat <- metadat %>% 
# #     mutate(celltype_sex = str_replace(celltype_sex, "activated", "activated-stellate"))
# #   metadat <- metadat %>% 
# #     mutate(celltype_sex = str_replace(celltype_sex, "quiescent", "quiescent-stellate"))
# #   metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex, "_", -1)
# #   combined_processed_atac@meta.data = metadat
# # }
# # 
# # table(combined_processed_atac@meta.data[["celltype"]])
# # table(combined_processed_atac@meta.data[["sex"]])
# # 
# # # cluster re-assignment occurs, which re-assigns clustering in my_levels
# # my_levels <- c("delta", "beta", "alpha", "gamma",
# #                "ductal", "acinar",
# #                "activated", "quiescent", "endothelial",
# #                "lymphocyte", "macrophage") 
# # 
# # table(combined_processed_atac$celltype)
# # 
# # # Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
# # combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
# # table(combined_processed_atac$celltype)
# # Idents(combined_processed_atac) <- "celltype"
# # DefaultAssay(combined_processed_atac) <- "ATAC"
# 
# # split object by cluster, take counts and aggregate by sum (pseudobulking)
# processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")
# hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
# 
# # Idents(processed_rna) <- "Sex"
# # unique(processed_rna$Sex)
# # processed_rna$sex <- plyr::mapvalues(
# #   x= processed_rna$Sex,
# #   from = c("F",
# #            "M"),
# #   to = c("female",
# #          "male"))
# # table(processed_rna@meta.data[["Sex"]])
# # table(processed_rna@meta.data[["sex"]])
# # Idents(processed_rna) <- "celltype_qadir"
# # processed_rna$celltype_sex <- paste(processed_rna$celltype_qadir, processed_rna$sex, sep = "_")
# Idents(processed_rna) <- "tissue_source"
# processed_rna <- subset(processed_rna, idents = c("Tulane"))
# Idents(processed_rna) <- "celltype_qadir"
# gc()
# gc()
# 
# #processed_rna <- subset(processed_rna, features = rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features))
# pblist_seu <- lapply(unique(processed_rna$celltype_qadir), function(x) {
#   rowSums(SeuratObject::GetAssayData(processed_rna[,processed_rna$celltype_qadir == x], "data"))
# })
# names(pblist_seu) = unique(processed_rna$celltype_qadir)
# 
# # DefaultAssay(hm.integrated.dfree) <- "RNA"
# # Idents(hm.integrated.dfree) <- "celltype"
# # hm.integrated.dfree$celltype_sex <- paste(hm.integrated.dfree$celltype, hm.integrated.dfree$sex, sep = "_")
# DefaultAssay(hm.integrated.dfree) <- "predicted"
# Idents(hm.integrated.dfree) <- "celltype"
# #hm.integrated.dfree <- subset(hm.integrated.dfree, features = rownames(processed_rna@assays[["RNA"]]@meta.features))
# pblist_sig <- lapply(unique(hm.integrated.dfree$celltype), function(x) {
#   rowSums(SeuratObject::GetAssayData(hm.integrated.dfree[,hm.integrated.dfree$celltype == x], "data"))
# })
# names(pblist_sig) = unique(hm.integrated.dfree$celltype)
# 
# # Subset all lists
# # Create a intersect of gene/peak names
# common_genes <- intersect(intersect(rownames(processed_rna@assays[["RNA"]]@meta.features), rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features)), processed_rna@assays[["RNA"]]@var.features)
# common_genes <- intersect(rownames(processed_rna@assays[["RNA"]]@meta.features), rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features))
# common_genes <- c("INS", "GCG", "SST", "PPY")
# common_genes <- processed_rna@assays[["RNA"]]@var.features
# common_genes <- intersect(rownames(hm.integrated.dfree@assays[["RNA"]]@meta.features), processed_rna@assays[["RNA"]]@var.features)
# 
# # Convert list of matrices to DFs
# ##### ATAC ###
# #####      ###
# pblist_sig <- lapply(
#   X = pblist_sig,
#   FUN = function(x){as.data.frame(x)})
# 
# # Initialize an empty list
# subset_dfs <- list()
# 
# # Iterate over each data frame
# for (i in 1:length(pblist_sig)) {
#   # Create subset based on common genes
#   subset_df <- pblist_sig[[i]][common_genes, ]
# 
#   # Append the subset data frame to the list
#   subset_dfs[[i]] <- subset_df
# }
# 
# # Add names
# names(subset_dfs) = unique(hm.integrated.dfree$celltype)
# pblist_sig <- subset_dfs
# 
# pblist_sig <- lapply(
#   X = pblist_sig,
#   FUN = function(x){as.matrix(x)})
# 
# 
# ##### scRNA ###
# #####       ###
# # Convert list of matrices to DFs
# pblist_seu <- lapply(
#   X = pblist_seu,
#   FUN = function(x){as.data.frame(x)})
# 
# # Initialize an empty list
# subset_dfs <- list()
# 
# # Iterate over each data frame
# for (i in 1:length(pblist_seu)) {
#   # Create subset based on common genes
#   subset_df <- pblist_seu[[i]][common_genes, ]
# 
#   # Append the subset data frame to the list
#   subset_dfs[[i]] <- subset_df
# }
# 
# # Add names
# names(subset_dfs) = unique(processed_rna$celltype_qadir)
# pblist_seu <- subset_dfs
# 
# pblist_seu <- lapply(
#   X = pblist_seu,
#   FUN = function(x){as.matrix(x)})
# 
# # create pairwise combinations
# ids = expand.grid(unique(processed_rna$celltype_qadir), unique(hm.integrated.dfree$celltype))
# 
# # iterate over combinations to extract correlations
# cors = apply(ids, 1, function(x) cor(pblist_seu[[x[1]]], pblist_sig[[x[2]]]))
# 
# pblist_sig[[2]]
# # fold the correlations vector into a matrix
# cmat = matrix(cors, nrow = length(unique(processed_rna$celltype_qadir)), byrow = TRUE)
# 
# # add cluster names. This implies that the names are matching, otherwise you need to add them separately for RNA and ATAC. This in turn depends on the order in which you expanded the grid above.
# rownames(cmat) = unique(processed_rna$celltype_qadir)
# colnames(cmat) = unique(hm.integrated.dfree$celltype)
# 
# # plot however you want
# pheatmap(cmat, clustering_distance_rows = "euclidean", 
#          cluster_rows = FALSE,cluster_cols = FALSE)
# 
# rna.seurat <- processed_rna[, sample(colnames(processed_rna), size =500, replace=F)]
# atac.signac <- hm.integrated.dfree[, sample(colnames(hm.integrated.dfree), size =500, replace=F)]
# qsave(rna.seurat, r"(E:\downsampling_collab\rna.seurat.qs)")
# qsave(atac.signac, r"(E:\downsampling_collab\atac.signac.qs)")



# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "acinar", "ductal",
               "quiescent_stellate", "activated_stellate", "endothelial",
               "lymphocyte", "macrophage")

table(hm.integrated.dfree$celltype)




# Confusion Matrix
predictions <- table(hm.integrated.dfree$celltype, hm.integrated.dfree$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
head(predictions)
predictions <- as.data.frame(predictions)
# predictions <- subset(predictions, Var2 == c("beta", "delta", "alpha", "gamma", 
#                                              "acinar", "ductal",
#                                              "quiescent_stellate", "activated_stellate", "endothelial",
#                                              "lymphocyte", "macrophage"))
head(predictions)
predicted_scale <- as.numeric(predictions$Freq)
predictions$Freq2 <- scale(predicted_scale, center = TRUE, scale = TRUE)

ggplot(predictions, aes(x=factor(Var1, level=c("delta", "beta", "alpha", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophage", "endothelial", "lymphocyte")), 
                        y=factor(Var2, level=c("delta", "beta+delta", "beta", "beta+alpha", "alpha", "cycling_endo", "gamma", "epsilon", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophages", "endothelial", "lymphocyte", "schwann", "mast")), 
                        fill = Freq2)) + 
  geom_tile() + scale_fill_gradient2(low = "dodgerblue3", high = "darkred", 
                                     midpoint = 0, limit = c(-5,5), space = "Lab", 
                                     name="Scaled Cell Contribution") +
  xlab("Cell type annotation (ATAC)") + ylab("Predicted cell type label (RNA)") +
  theme_minimal()+ 
  geom_text(aes(label = round(Freq, 1))) +
  ggtitle("Normalised Confusion Matrix") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))


# Plotting Chromatin accessibility
# when copying off the heatmap make sure to adjust dashes
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
qsave(conns, r"(E:\2.SexbasedStudyCurrent\QS files\conns.qs)")
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
#qsave(hm.integrated.dfree, r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")



DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype"
Idents(hm.integrated.dfree) <- "sex"
male_set <- subset(hm.integrated.dfree, idents = (c("male")))
female_set <- subset(hm.integrated.dfree, idents = (c("female")))
Idents(male_set) <- "celltype"
Idents(female_set) <- "celltype"
XIST_f <-  CoveragePlot(female_set, region = c("XIST"), 
                    assay = "macs2",
                    window = 100,
                    ymax = 100,
             links = FALSE,
             #tile = TRUE,
             extend.upstream = 20000,
             extend.downstream = 20000
             ) & scale_fill_manual(values = c("dodgerblue3",      #beta
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
XIST_m + XIST_f
kdm_m + kdm_f 
+ kdmd_m + kdmd_f
p1
#p1+p2

link_plot <- LinkPlot(
  object = hm.integrated.dfree,
  assay = "macs2",
  region = c("SDS"),
  scale.linewidth = FALSE,
  min.cutoff = 0,
  extend.upstream = 20000,
  extend.downstream = 20000
) + scale_color_gradient2(low = "white", # whatever colors you want here
                          mid = "white",
                          high = "darkred",
                          na.value = "grey50")

expr_plot <- ExpressionPlot(
  object = hm.integrated.dfree,
  features = "SDS",
  assay = "predicted"
) & scale_fill_manual(values = c("dodgerblue3",      #beta
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

CombineTracks(
  plotlist = list(p1, link_plot),
  expression.plot = expr_plot,
  heights = c(10, 2),
  widths = c(10, 1)
)


# Plotting heatmap for all significant DA regions
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")

# All genes
beta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\beta_peaks.csv)", sep = ",", row.names = 1)
alpha_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\alpha_peaks.csv)", sep = ",", row.names = 1)
delta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\delta_peaks.csv)", sep = ",", row.names = 1)
gamma_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\gamma_peaks.csv)", sep = ",", row.names = 1)
acinar_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\acinar_peaks.csv)", sep = ",", row.names = 1)
ductal_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\ductal_peaks.csv)", sep = ",", row.names = 1)
quiescentstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\quiescentstellate_peaks.csv)", sep = ",", row.names = 1)
activatedstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\activatedstellate_peaks.csv)", sep = ",", row.names = 1)
macrophage_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\macrophage_peaks.csv)", sep = ",", row.names = 1)
lymphocyte_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\lymphocyte_peaks.csv)", sep = ",", row.names = 1)
endothelial_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\endothelial_peaks.csv)", sep = ",", row.names = 1)

{
  beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-302
  beta_peaks <- dplyr::filter(beta_peaks, p_val_adj < 1e-5) 
  open_beta_peaks <- rownames(beta_peaks[beta_peaks$avg_log2FC > 1, ])
  
  alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-302
  alpha_peaks <- dplyr::filter(alpha_peaks, p_val_adj < 1e-5) 
  open_alpha_peaks <- rownames(alpha_peaks[alpha_peaks$avg_log2FC > 1, ])
  
  delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-302
  delta_peaks <- dplyr::filter(delta_peaks, p_val_adj < 1e-5) 
  open_delta_peaks <- rownames(delta_peaks[delta_peaks$avg_log2FC > 1, ])
  
  gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-302
  gamma_peaks <- dplyr::filter(gamma_peaks, p_val_adj < 1e-5) 
  open_gamma_peaks <- rownames(gamma_peaks[gamma_peaks$avg_log2FC > 1, ])
  
  acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-302
  acinar_peaks <- dplyr::filter(acinar_peaks, p_val_adj < 1e-5) 
  open_acinar_peaks <- rownames(acinar_peaks[acinar_peaks$avg_log2FC > 1, ])
  
  ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-302
  ductal_peaks <- dplyr::filter(ductal_peaks, p_val_adj < 1e-5) 
  open_ductal_peaks <- rownames(ductal_peaks[ductal_peaks$avg_log2FC > 1, ])
  
  quiescentstellate_peaks$p_val_adj[quiescentstellate_peaks$p_val_adj == 0] <- 2e-302
  quiescentstellate_peaks <- dplyr::filter(quiescentstellate_peaks, p_val_adj < 1e-5) 
  open_quiescentstellate_peaks <- rownames(quiescentstellate_peaks[quiescentstellate_peaks$avg_log2FC > 1, ])
  
  activatedstellate_peaks$p_val_adj[activatedstellate_peaks$p_val_adj == 0] <- 2e-302
  activatedstellate_peaks <- dplyr::filter(activatedstellate_peaks, p_val_adj < 1e-5) 
  open_activatedstellate_peaks <- rownames(activatedstellate_peaks[activatedstellate_peaks$avg_log2FC > 1, ])
  
  macrophage_peaks$p_val_adj[macrophage_peaks$p_val_adj == 0] <- 2e-302
  macrophage_peaks <- dplyr::filter(macrophage_peaks, p_val_adj < 1e-5) 
  open_macrophage_peaks <- rownames(macrophage_peaks[macrophage_peaks$avg_log2FC > 1, ])
  
  lymphocyte_peaks$p_val_adj[lymphocyte_peaks$p_val_adj == 0] <- 2e-302
  lymphocyte_peaks <- dplyr::filter(lymphocyte_peaks, p_val_adj < 1e-5) 
  open_lymphocyte_peaks <- rownames(lymphocyte_peaks[lymphocyte_peaks$avg_log2FC > 1, ])
  
  endothelial_peaks$p_val_adj[endothelial_peaks$p_val_adj == 0] <- 2e-302
  endothelial_peaks <- dplyr::filter(endothelial_peaks, p_val_adj < 1e-5) 
  open_endothelial_peaks <- rownames(endothelial_peaks[endothelial_peaks$avg_log2FC > 1, ])
  
  cgenes_beta <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks)
  cgenes_alpha <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks)
  cgenes_delta <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks)
  cgenes_gamma <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks)
  cgenes_acinar <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks)
  cgenes_ductal <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks)
  cgenes_qstel <- ClosestFeature(hm.integrated.dfree, regions = open_quiescentstellate_peaks)
  cgenes_astel <- ClosestFeature(hm.integrated.dfree, regions = open_activatedstellate_peaks)
  cgenes_macro <- ClosestFeature(hm.integrated.dfree, regions = open_macrophage_peaks)
  cgenes_lympho <- ClosestFeature(hm.integrated.dfree, regions = open_lymphocyte_peaks)
  cgenes_endo <- ClosestFeature(hm.integrated.dfree, regions = open_endothelial_peaks)
  
  cgenes_beta <- dplyr::filter(cgenes_beta, distance < 100000) 
  cgenes_alpha <- dplyr::filter(cgenes_alpha, distance < 100000) 
  cgenes_delta <- dplyr::filter(cgenes_delta, distance < 100000) 
  cgenes_gamma <- dplyr::filter(cgenes_gamma, distance < 100000) 
  cgenes_acinar <- dplyr::filter(cgenes_acinar, distance < 100000) 
  cgenes_ductal <- dplyr::filter(cgenes_ductal, distance < 100000) 
  cgenes_qstel <- dplyr::filter(cgenes_qstel, distance < 100000) 
  cgenes_astel <- dplyr::filter(cgenes_astel, distance < 100000) 
  cgenes_macro <- dplyr::filter(cgenes_macro, distance < 100000) 
  cgenes_lympho <- dplyr::filter(cgenes_lympho, distance < 100000) 
  cgenes_endo <- dplyr::filter(cgenes_endo, distance < 100000)
}

allregions <- c(
  as.character(cgenes_beta$query_region),
  as.character(cgenes_delta$query_region),
  as.character(cgenes_alpha$query_region),
  as.character(cgenes_gamma$query_region),
  as.character(cgenes_ductal$query_region),
  as.character(cgenes_acinar$query_region),
  as.character(cgenes_astel$query_region),
  as.character(cgenes_qstel$query_region),
  as.character(cgenes_endo$query_region),
  as.character(cgenes_lympho$query_region),
  as.character(cgenes_macro$query_region)
)

cgenes_beta <- distinct(cgenes_beta, gene_name, .keep_all = TRUE)
cgenes_alpha <- distinct(cgenes_alpha, gene_name, .keep_all = TRUE)
cgenes_delta <- distinct(cgenes_delta, gene_name, .keep_all = TRUE)
cgenes_gamma <- distinct(cgenes_gamma, gene_name, .keep_all = TRUE)
cgenes_acinar <- distinct(cgenes_acinar, gene_name, .keep_all = TRUE)
cgenes_ductal <- distinct(cgenes_ductal, gene_name, .keep_all = TRUE)
cgenes_qstel <- distinct(cgenes_qstel, gene_name, .keep_all = TRUE)
cgenes_astel <- distinct(cgenes_astel, gene_name, .keep_all = TRUE)
cgenes_macro <- distinct(cgenes_macro, gene_name, .keep_all = TRUE)
cgenes_lympho <- distinct(cgenes_lympho, gene_name, .keep_all = TRUE)
cgenes_endo <- distinct(cgenes_endo, gene_name, .keep_all = TRUE)


allgenes <- c(
  as.character(cgenes_beta$gene_name),
  as.character(cgenes_delta$gene_name),
  as.character(cgenes_alpha$gene_name),
  as.character(cgenes_gamma$gene_name),
  as.character(cgenes_ductal$gene_name),
  as.character(cgenes_acinar$gene_name),
  as.character(cgenes_astel$gene_name),
  as.character(cgenes_qstel$gene_name),
  as.character(cgenes_endo$gene_name),
  as.character(cgenes_lympho$gene_name),
  as.character(cgenes_macro$gene_name)
)

# Sex differences
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")
DefaultAssay(hm.integrated.dfree) <- "macs2"
beta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\beta_mvf_peaks.csv)", sep = ",", row.names = 1)
alpha_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\alpha_mvf_peaks.csv)", sep = ",", row.names = 1)
delta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\delta_mvf_peaks.csv)", sep = ",", row.names = 1)
gamma_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\gamma_mvf_peaks.csv)", sep = ",", row.names = 1)
acinar_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\acinar_mvf_peaks.csv)", sep = ",", row.names = 1)
ductal_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\ductal_mvf_peaks.csv)", sep = ",", row.names = 1)
qstell_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\quiescentstellate_mvf_peaks.csv)", sep = ",", row.names = 1)
astell_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\activatedstellate_mvf_peaks.csv)", sep = ",", row.names = 1)
macro_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\macrophage_mvf_peaks.csv)", sep = ",", row.names = 1)
lympho_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\lymphocyte_mvf_peaks.csv)", sep = ",", row.names = 1)
endo_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\endothelial_mvf_peaks.csv)", sep = ",", row.names = 1)

# Plot ranges for Ruth
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype"
x_cells <- subset(hm.integrated.dfree, idents = c("delta"))
x_cells <- subset(hm.integrated.dfree, idents = c("endothelial"))
x_cells <- subset(hm.integrated.dfree, idents = c("acinar"))

ranges.show <- StringToGRanges("chr2-121035441-121037137")
CoveragePlot(ductal_cells, region = c("chr2-121035441-121037137"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr19-19627168-19629130")
CoveragePlot(lymphocytes, region = c("chr19-19627168-19629130"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr5-68433822-68434744")
CoveragePlot(x_cells, region = c("chr5-68433822-68434744"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr7-44111586-44113624")
CoveragePlot(x_cells, region = c("chr7-44111586-44113624"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr12-121128766-121129441")
CoveragePlot(x_cells, region = c("chr12-121128766-121129441"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr11-3159900-3161041")
CoveragePlot(x_cells, region = c("chr11-3159900-3161041"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr10-69657288-69657771")
CoveragePlot(x_cells, region = c("chr10-69657288-69657771"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)

ranges.show <- StringToGRanges("chr8-117376094-117376998")
CoveragePlot(x_cells, region = c("chr8-117376094-117376998"), 
             extend.upstream = 50000,
             extend.downstream = 50000, links = TRUE, split.by = "sex", region.highlight = ranges.show)



beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-302
beta_peaks <- dplyr::filter(beta_peaks, p_val_adj < 5e-2) 
open_beta_peaks_male <- rownames(beta_peaks[beta_peaks$avg_log2FC > 0 & beta_peaks$p_val_adj < 5e-2, ])
open_beta_peaks_female <- rownames(beta_peaks[beta_peaks$avg_log2FC < 0 & beta_peaks$p_val_adj < 5e-2, ])
length(open_beta_peaks_male)
length(open_beta_peaks_female)

alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-302
alpha_peaks <- dplyr::filter(alpha_peaks, p_val_adj < 5e-2) 
open_alpha_peaks_male <- rownames(alpha_peaks[alpha_peaks$avg_log2FC > 0 & alpha_peaks$p_val_adj < 5e-2, ])
open_alpha_peaks_female <- rownames(alpha_peaks[alpha_peaks$avg_log2FC < 0 & alpha_peaks$p_val_adj < 5e-2, ])
length(open_alpha_peaks_male)
length(open_alpha_peaks_female)

cgenes_beta_male <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_male)
cgenes_beta_female <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_female)

cgenes_beta_male <- dplyr::filter(cgenes_beta_male, distance < 100000) 
cgenes_beta_female <- dplyr::filter(cgenes_beta_female, distance < 100000) 

{
beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-302
beta_peaks <- dplyr::filter(beta_peaks, p_val_adj < 5e-2) 
open_beta_peaks_male <- rownames(beta_peaks[beta_peaks$avg_log2FC > 1 & beta_peaks$p_val_adj < 5e-2, ])
open_beta_peaks_female <- rownames(beta_peaks[beta_peaks$avg_log2FC < -1 & beta_peaks$p_val_adj < 5e-2, ])

alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-302
alpha_peaks <- dplyr::filter(alpha_peaks, p_val_adj < 5e-2) 
open_alpha_peaks_male <- rownames(alpha_peaks[alpha_peaks$avg_log2FC > 1 & alpha_peaks$p_val_adj < 5e-2, ])
open_alpha_peaks_female <- rownames(alpha_peaks[alpha_peaks$avg_log2FC < -1 & alpha_peaks$p_val_adj < 5e-2, ])

delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-302
delta_peaks <- dplyr::filter(delta_peaks, p_val_adj < 5e-2) 
open_delta_peaks_male <- rownames(delta_peaks[delta_peaks$avg_log2FC > 1 & delta_peaks$p_val_adj < 5e-2, ])
open_delta_peaks_female <- rownames(delta_peaks[delta_peaks$avg_log2FC < -1 & delta_peaks$p_val_adj < 5e-2, ])

gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-302
gamma_peaks <- dplyr::filter(gamma_peaks, p_val_adj < 5e-2) 
open_gamma_peaks_male <- rownames(gamma_peaks[gamma_peaks$avg_log2FC > 1 & gamma_peaks$p_val_adj < 5e-2, ])
open_gamma_peaks_female <- rownames(gamma_peaks[gamma_peaks$avg_log2FC < -1 & gamma_peaks$p_val_adj < 5e-2, ])

ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-302
ductal_peaks <- dplyr::filter(ductal_peaks, p_val_adj < 5e-2) 
open_ductal_peaks_male <- rownames(ductal_peaks[ductal_peaks$avg_log2FC > 1 & ductal_peaks$p_val_adj < 5e-2, ])
open_ductal_peaks_female <- rownames(ductal_peaks[ductal_peaks$avg_log2FC < -1 & ductal_peaks$p_val_adj < 5e-2, ])

acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-302
acinar_peaks <- dplyr::filter(acinar_peaks, p_val_adj < 5e-2) 
open_acinar_peaks_male <- rownames(acinar_peaks[acinar_peaks$avg_log2FC > 1 & acinar_peaks$p_val_adj < 5e-2, ])
open_acinar_peaks_female <- rownames(acinar_peaks[acinar_peaks$avg_log2FC < -1 & acinar_peaks$p_val_adj < 5e-2, ])

qstell_peaks$p_val_adj[qstell_peaks$p_val_adj == 0] <- 2e-302
qstell_peaks <- dplyr::filter(qstell_peaks, p_val_adj < 5e-2) 
open_qstell_peaks_male <- rownames(qstell_peaks[qstell_peaks$avg_log2FC > 1 & qstell_peaks$p_val_adj < 5e-2, ])
open_qstell_peaks_female <- rownames(qstell_peaks[qstell_peaks$avg_log2FC < -1 & qstell_peaks$p_val_adj < 5e-2, ])

astell_peaks$p_val_adj[astell_peaks$p_val_adj == 0] <- 2e-302
astell_peaks <- dplyr::filter(astell_peaks, p_val_adj < 5e-2) 
open_astell_peaks_male <- rownames(astell_peaks[astell_peaks$avg_log2FC > 1 & astell_peaks$p_val_adj < 5e-2, ])
open_astell_peaks_female <- rownames(astell_peaks[astell_peaks$avg_log2FC < -1 & astell_peaks$p_val_adj < 5e-2, ])

macro_peaks$p_val_adj[macro_peaks$p_val_adj == 0] <- 2e-302
macro_peaks <- dplyr::filter(macro_peaks, p_val_adj < 5e-2) 
open_macro_peaks_male <- rownames(macro_peaks[macro_peaks$avg_log2FC > 1 & macro_peaks$p_val_adj < 5e-2, ])
open_macro_peaks_female <- rownames(macro_peaks[macro_peaks$avg_log2FC < -1 & macro_peaks$p_val_adj < 5e-2, ])

lympho_peaks$p_val_adj[lympho_peaks$p_val_adj == 0] <- 2e-302
lympho_peaks <- dplyr::filter(lympho_peaks, p_val_adj < 5e-2) 
open_lympho_peaks_male <- rownames(lympho_peaks[lympho_peaks$avg_log2FC > 1 & lympho_peaks$p_val_adj < 5e-2, ])
open_lympho_peaks_female <- rownames(lympho_peaks[lympho_peaks$avg_log2FC < -1 & lympho_peaks$p_val_adj < 5e-2, ])

endo_peaks$p_val_adj[endo_peaks$p_val_adj == 0] <- 2e-302
endo_peaks <- dplyr::filter(endo_peaks, p_val_adj < 5e-2) 
open_endo_peaks_male <- rownames(endo_peaks[endo_peaks$avg_log2FC > 1 & endo_peaks$p_val_adj < 5e-2, ])
open_endo_peaks_female <- rownames(endo_peaks[endo_peaks$avg_log2FC < -1 & endo_peaks$p_val_adj < 5e-2, ])

cgenes_beta_male <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_male)
cgenes_beta_female <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks_female)
cgenes_alpha_male <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks_male)
cgenes_alpha_female <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks_female)
cgenes_delta_male <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks_male)
cgenes_delta_female <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks_female)
cgenes_gamma_male <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks_male)
cgenes_gamma_female <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks_female)
cgenes_ductal_male <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks_male)
cgenes_ductal_female <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks_female)
cgenes_acinar_male <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks_male)
cgenes_acinar_female <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks_female)
cgenes_qstell_male <- ClosestFeature(hm.integrated.dfree, regions = open_qstell_peaks_male)
cgenes_qstell_female <- ClosestFeature(hm.integrated.dfree, regions = open_qstell_peaks_female)
cgenes_astell_male <- ClosestFeature(hm.integrated.dfree, regions = open_astell_peaks_male)
cgenes_astell_female <- ClosestFeature(hm.integrated.dfree, regions = open_astell_peaks_female)
cgenes_macro_male <- ClosestFeature(hm.integrated.dfree, regions = open_macro_peaks_male)
cgenes_macro_female <- ClosestFeature(hm.integrated.dfree, regions = open_macro_peaks_female)
cgenes_lympho_male <- ClosestFeature(hm.integrated.dfree, regions = open_lympho_peaks_male)
#cgenes_lympho_female <- ClosestFeature(hm.integrated.dfree, regions = open_lympho_peaks_female)
cgenes_endo_male <- ClosestFeature(hm.integrated.dfree, regions = open_endo_peaks_male)
cgenes_endo_female <- ClosestFeature(hm.integrated.dfree, regions = open_endo_peaks_female)


cgenes_beta_male <- dplyr::filter(cgenes_beta_male, distance < 100000) 
cgenes_beta_female <- dplyr::filter(cgenes_beta_female, distance < 100000) 
cgenes_alpha_male <- dplyr::filter(cgenes_alpha_male, distance < 100000) 
cgenes_alpha_female <- dplyr::filter(cgenes_alpha_female, distance < 100000) 
cgenes_delta_male <- dplyr::filter(cgenes_delta_male, distance < 100000) 
cgenes_delta_female <- dplyr::filter(cgenes_delta_female, distance < 100000)
cgenes_gamma_male <- dplyr::filter(cgenes_gamma_male, distance < 100000) 
cgenes_gamma_female <- dplyr::filter(cgenes_gamma_female, distance < 100000)
cgenes_ductal_male <- dplyr::filter(cgenes_ductal_male, distance < 100000) 
cgenes_ductal_female <- dplyr::filter(cgenes_ductal_female, distance < 100000)
cgenes_acinar_male <- dplyr::filter(cgenes_acinar_male, distance < 100000) 
cgenes_acinar_female <- dplyr::filter(cgenes_acinar_female, distance < 100000)
cgenes_qstell_male <- dplyr::filter(cgenes_qstell_male, distance < 100000) 
cgenes_qstell_female <- dplyr::filter(cgenes_qstell_female, distance < 100000)
cgenes_astell_male <- dplyr::filter(cgenes_astell_male, distance < 100000) 
cgenes_astell_female <- dplyr::filter(cgenes_astell_female, distance < 100000)
cgenes_macro_male <- dplyr::filter(cgenes_macro_male, distance < 100000) 
cgenes_macro_female <- dplyr::filter(cgenes_macro_female, distance < 100000)
cgenes_lympho_male <- dplyr::filter(cgenes_lympho_male, distance < 100000) 
# cgenes_lympho_female <- dplyr::filter(cgenes_lympho_female, distance < 100000)
cgenes_endo_male <- dplyr::filter(cgenes_endo_male, distance < 100000) 
cgenes_endo_female <- dplyr::filter(cgenes_endo_female, distance < 100000)
}

allregions <- c(
  as.character(cgenes_beta_male$query_region),
  as.character(cgenes_beta_female$query_region),
  as.character(cgenes_alpha_male$query_region),
  as.character(cgenes_alpha_female$query_region),
  as.character(cgenes_delta_male$query_region),
  as.character(cgenes_delta_female$query_region),
  as.character(cgenes_gamma_male$query_region),
  as.character(cgenes_gamma_female$query_region),
  as.character(cgenes_ductal_male$query_region),
  as.character(cgenes_ductal_female$query_region),
  as.character(cgenes_acinar_male$query_region),
  as.character(cgenes_acinar_female$query_region),
  as.character(cgenes_astell_male$query_region),
  as.character(cgenes_astell_female$query_region),
  as.character(cgenes_qstell_male$query_region),
  as.character(cgenes_qstell_female$query_region),
  as.character(cgenes_macro_male$query_region),
  as.character(cgenes_macro_female$query_region),
  as.character(cgenes_lympho_male$query_region),
  #as.character(cgenes_lympho_female$query_region),
  as.character(cgenes_endo_male$query_region),
  as.character(cgenes_endo_female$query_region)
)

cgenes_beta_male <- distinct(cgenes_beta_male, gene_name, .keep_all = TRUE)
cgenes_beta_female <- distinct(cgenes_beta_female, gene_name, .keep_all = TRUE)
cgenes_alpha_male <- distinct(cgenes_alpha_male, gene_name, .keep_all = TRUE)
cgenes_alpha_female <- distinct(cgenes_alpha_female, gene_name, .keep_all = TRUE)
cgenes_delta_male <- distinct(cgenes_delta_male, gene_name, .keep_all = TRUE)
cgenes_delta_female <- distinct(cgenes_delta_female, gene_name, .keep_all = TRUE)
cgenes_gamma_male <- distinct(cgenes_gamma_male, gene_name, .keep_all = TRUE)
cgenes_gamma_female <- distinct(cgenes_gamma_female, gene_name, .keep_all = TRUE)
cgenes_ductal_male <- distinct(cgenes_ductal_male, gene_name, .keep_all = TRUE)
cgenes_ductal_female <- distinct(cgenes_ductal_female, gene_name, .keep_all = TRUE)
cgenes_acinar_male <- distinct(cgenes_acinar_male, gene_name, .keep_all = TRUE)
cgenes_acinar_female <- distinct(cgenes_acinar_female, gene_name, .keep_all = TRUE)
cgenes_astell_male <- distinct(cgenes_astell_male, gene_name, .keep_all = TRUE)
cgenes_astell_female <- distinct(cgenes_astell_female, gene_name, .keep_all = TRUE)
cgenes_qstell_male <- distinct(cgenes_qstell_male, gene_name, .keep_all = TRUE)
cgenes_qstell_female <- distinct(cgenes_qstell_female, gene_name, .keep_all = TRUE)
cgenes_macro_male <- distinct(cgenes_macro_male, gene_name, .keep_all = TRUE)
cgenes_macro_female <- distinct(cgenes_macro_female, gene_name, .keep_all = TRUE)
cgenes_lympho_male <- distinct(cgenes_lympho_male, gene_name, .keep_all = TRUE)
#cgenes_lympho_female <- distinct(cgenes_lympho_female, gene_name, .keep_all = TRUE)
cgenes_endo_male <- distinct(cgenes_endo_male, gene_name, .keep_all = TRUE)
cgenes_endo_female <- distinct(cgenes_endo_female, gene_name, .keep_all = TRUE)


allgenes <- c(
  as.character(cgenes_beta_male$gene_name),
  as.character(cgenes_beta_female$gene_name),
  as.character(cgenes_alpha_male$gene_name),
  as.character(cgenes_alpha_female$gene_name),
  as.character(cgenes_delta_male$gene_name),
  as.character(cgenes_delta_female$gene_name),
  as.character(cgenes_gamma_male$gene_name),
  as.character(cgenes_gamma_female$gene_name),
  as.character(cgenes_ductal_male$gene_name),
  as.character(cgenes_ductal_female$gene_name),
  as.character(cgenes_acinar_male$gene_name),
  as.character(cgenes_acinar_female$gene_name),
  as.character(cgenes_astell_male$gene_name),
  as.character(cgenes_astell_female$gene_name),
  as.character(cgenes_qstell_male$gene_name),
  as.character(cgenes_qstell_female$gene_name),
  as.character(cgenes_macro_male$gene_name),
  as.character(cgenes_macro_female$gene_name),
  as.character(cgenes_lympho_male$gene_name),
# as.character(cgenes_lympho_female$gene_name),
  as.character(cgenes_endo_male$gene_name),
  as.character(cgenes_endo_female$gene_name)
)

allregions
allgenes

# Concatenate and remove dupliates
uniquegenes <- unique(allgenes)
uniqueregions <- unique(allregions)

plotgenesnow <- intersect(combined_processed_atac@assays[["RNA_macs2"]]@counts@Dimnames[[1]], uniquegenes)

# Heatmap
# Pseudobulk
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry_lib"
combined_processed_atac <- AverageExpression(hm.integrated.dfree, 
                                               assays = c("macs2", "RNA_macs2"), 
                                               features = NULL, return.seurat = TRUE,  
                                               group.by = "celltype_sex_ancestry_lib",
                                               slot = "counts", verbose = FALSE)

combined_processed_atac$celltype_sex_ancestry_lib <- Cells(combined_processed_atac) #6892 Seurat

{
  Idents(combined_processed_atac) <- 'celltype_sex_ancestry_lib'
  combined_processed_atac$celltype <- combined_processed_atac@meta.data[["orig.ident"]]
  metadat <- combined_processed_atac@meta.data
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "quiescent_stellate", "quiescent-stellate"))
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex_ancestry_lib, "_", -3)
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -2)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -1)
  combined_processed_atac@meta.data = metadat
}


table(combined_processed_atac@meta.data[["celltype"]])
table(combined_processed_atac@meta.data[["sex"]])

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "ductal", "acinar",
               "activated", "quiescent", "endothelial",
               "lymphocyte", "macrophage")

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
table(combined_processed_atac$celltype)
Idents(combined_processed_atac) <- "celltype"


# Plot heatmap
genes.to.plot <- uniquegenes
label_genes <- c("INS", "GCG", "SST") #uniquegenes
all_genes_inobj <- rownames(combined_processed_atac@assays[["RNA"]])
genes.to.plot.intr <- intersect(all_genes_inobj, uniquegenes) # to get genes make a seurat object by combining data over counts
xistcheck <- intersect(combined_processed_atac@assays[["RNA_macs2"]]@counts@Dimnames[[1]], c("SNX5")) # to get genes make a seurat object by combining data over counts

regions.to.plot <- uniqueregions
label_genes <- c("INS", "GCG", "SST") #uniquegenes


write.csv(regions.to.plot, R"(C:\Users\mqadir\Desktop\regions.csv)")
DefaultAssay(combined_processed_atac) <- "RNA_macs2" #RNA
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 8,
    height = 6)
dittoHeatmap(
  object = combined_processed_atac,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = c("XIST"), #c("chr1-181244-181600", "chr1-184156-184547"), #uniqueregions, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("celltype", "sex", "ancestry"),
  order.by = c("sex", "celltype"),
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
                                        "activated" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophage" = "magenta3",
                                        "quiescent" = "salmon"),
                           sex = c("female" = "red4",
                                   "male" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black")),
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
  show_rownames = TRUE,
  # scale = "row",
  cluster_row = TRUE,
  cluster_cols = FALSE,
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
  #gaps_col = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150),
  #gaps_row = c(2980, 4760, 7790, 8090, 18900, 24850, 39540, 42385, 45649, 46290),
  #gaps_row = c(2985, 4760, 7790, 8090, 18900, 24850, 39540),
  #gaps_row = c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)
  #gaps_row = c(2804, 4584, 7464, 7799, 17833, 23493, 36665, 39029, 41669, 42089)
) 
dev.off()
dev.off()

+ rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_atac[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1
)
)


# GO plot This uses allgenes that was derived from above
# Compare

# create gene list
# Sort out top 1000 accessible sites
cgenes_beta <- cgenes_beta %>% slice_max(cgenes_beta, n = 1000)
cgenes_delta <- cgenes_delta %>% slice_max(cgenes_delta, n = 1000)
cgenes_alpha <- cgenes_alpha %>% slice_max(cgenes_alpha, n = 1000)
cgenes_gamma <- cgenes_gamma %>% slice_max(cgenes_gamma, n = 1000)
cgenes_ductal <- cgenes_ductal %>% slice_max(cgenes_ductal, n = 1000)
cgenes_acinar <- cgenes_acinar %>% slice_max(cgenes_acinar, n = 1000)
cgenes_astel <- cgenes_astel %>% slice_max(cgenes_astel, n = 1000)
cgenes_qstel <- cgenes_qstel %>% slice_max(cgenes_qstel, n = 1000)
cgenes_endo <- cgenes_endo %>% slice_max(cgenes_endo, n = 1000)
cgenes_lympho <- cgenes_lympho %>% slice_max(cgenes_lympho, n = 1000)
cgenes_macro <- cgenes_macro %>% slice_max(cgenes_macro, n = 1000)

# For accessibility remember to run the code again generating the _beta files and subset 100kb window peaks
cgenes_beta <- as.character(cgenes_beta$gene_name)
cgenes_delta <- as.character(cgenes_delta$gene_name)
cgenes_alpha <- as.character(cgenes_alpha$gene_name)
cgenes_gamma <- as.character(cgenes_gamma$gene_name)
cgenes_ductal <- as.character(cgenes_ductal$gene_name)
cgenes_acinar <- as.character(cgenes_acinar$gene_name)
cgenes_astel <- as.character(cgenes_astel$gene_name)
cgenes_qstel <- as.character(cgenes_qstel$gene_name)
cgenes_endo <- as.character(cgenes_endo$gene_name)
cgenes_lympho <- as.character(cgenes_lympho$gene_name)
cgenes_macro <- as.character(cgenes_macro$gene_name)

#across sex
beta_male_ac <- as.character(cgenes_beta_male$gene_name)
beta_female_ac <- as.character(cgenes_beta_female$gene_name)
alpha_male_ac <- as.character(cgenes_alpha_male$gene_name)
alpha_female_ac <- as.character(cgenes_alpha_female$gene_name)
delta_male_ac <- as.character(cgenes_delta_male$gene_name)
delta_female_ac <- as.character(cgenes_delta_female$gene_name)
gamma_male_ac <- as.character(cgenes_gamma_male$gene_name)
gamma_female_ac <-as.character(cgenes_gamma_female$gene_name)
ductal_male_ac <- as.character(cgenes_ductal_male$gene_name)
ductal_female_ac <- as.character(cgenes_ductal_female$gene_name)
acinar_male_ac <- as.character(cgenes_acinar_male$gene_name)
acinar_female_ac <- as.character(cgenes_acinar_female$gene_name)
astell_male_ac <- as.character(cgenes_astell_male$gene_name)
astell_female_ac <- as.character(cgenes_astell_female$gene_name)
macro_male_ac <- as.character(cgenes_macro_male$gene_name)
macro_female_ac <- as.character(cgenes_macro_female$gene_name)
endo_male_ac <- as.character(cgenes_endo_male$gene_name)
endo_female_ac <- as.character(cgenes_endo_female$gene_name)


gene.list <- list("beta" = cgenes_beta, "alpha" = cgenes_alpha, "delta" = cgenes_delta, "gamma" = cgenes_gamma, 
                  "acinar" = cgenes_acinar, "ductal" = cgenes_ductal, "qstel" = cgenes_qstel, "astel" = cgenes_astel,
                  "macro" = cgenes_macro, "lympho" = cgenes_lympho, "endo" = cgenes_endo)

gene.list <- list("beta_male_ac" = beta_male_ac, "alpha_male_ac" = alpha_male_ac, "delta_male_ac" = delta_male_ac, "gamma_male_ac" = gamma_male_ac, "ductal_male_ac" = ductal_male_ac, "acinar_male_ac" = acinar_male_ac, "astell_male_ac"= astell_male_ac, "macro_male_ac" = macro_male_ac, "endo_male_ac" = endo_male_ac,
                  "beta_female_ac" = beta_female_ac,  "alpha_female_ac" = alpha_female_ac,  "delta_female_ac" = delta_female_ac,  "gamma_female_ac" = gamma_female_ac, "ductal_female_ac" = ductal_female_ac, "acinar_female_ac" = acinar_female_ac, "astell_female_ac" = astell_female_ac, "macro_female_ac"=macro_female_ac, "endo_female_ac" = endo_female_ac)

# Unified genes for 
gene.list <- unique(allgenes)
gene.list <- intersect(combined_processed_atac@assays[["RNA"]]@counts@Dimnames[[1]], uniquegenes)

ck <- compareCluster(geneCluster = gene.list, # list of genes
                     fun = enrichGO, 
                     #universe = rownames(processed_rna@assays[["RNA"]]@counts), 
                     keyType = "SYMBOL", #keytypes(org.Hs.eg.db)
                     OrgDb = org.Hs.eg.db, 
                     ont = c("ALL"), 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1, #if not set default is at 0.05
                     readable = TRUE) 
head(ck) 
cluster_summary <- data.frame(ck)
ck <- ck[ck@compareClusterResult[["qvalue"]] < 0.1, asis=T]
dotplot(ck, showCategory = 8) + #coord_flip() + 
  scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
dotplot(ck, showCategory = 1)
ck.save <- ck@compareClusterResult
write.csv(ck.save, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\ORA\sex_accessibility.save.csv)")

dotplot(ck, 
        showCategory = c("insulin secretion", "calcium ion homeostasis", "protein localization to plasma membrane", "glucose homeostasis",
                             "regulation of G protein-coupled receptor signaling pathway", "pancreas development", "hormone metabolic process",
                             "glutamate receptor signaling pathway", "synaptic transmission, GABAergic", "peptide secretion",
                             "hormone transport",
                             "developmental growth involved in morphogenesis", "sodium ion transport", "regulation of actin filament-based process",
                             "regulation of Wnt signaling pathway", "integrin-mediated signaling pathway",
                             "muscle contraction", "establishment of endothelial barrier", "cellular response to fibroblast growth factor stimulus", "extracellular matrix organization",
                             "activation of immune response", "response to interferon-gamma", "cytokine-mediated signaling pathway", "leukocyte chemotaxis", "leukocyte degranulation", "interleukin-1 beta production", "phagocytosis",
                             "alpha-beta T cell activation",
                             "endothelial cell development", "cellular response to angiotensin"),
        # showCategory = c("histone lysine demethylation", "protein dealkylation", "histone modification", "sequestering of actin monomers", "gonad development", "sex determination", 
        #                  "androgen receptor signaling pathway", 
        #                  "dosage compensation", "regulation of gene expression, epigenetic", "chromatin remodeling"),
        font.size=14) + coord_flip() + 
  scale_x_discrete(limits=rev) + scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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


# Unified genes for plotting dotplot of genes intersecting DARs
unique.de.genes.all # Data comes from Line 530-535 this is a list of unique DE genes 

gene.list <- list("beta_male_ac" = beta_male_ac, "alpha_male_ac" = alpha_male_ac, "delta_male_ac" = delta_male_ac, "gamma_male_ac" = gamma_male_ac, "ductal_male_ac" = ductal_male_ac, "acinar_male_ac" = acinar_male_ac, "astell_male_ac"= astell_male_ac, "macro_male_ac" = macro_male_ac, "endo_male_ac" = endo_male_ac,
                  "beta_female_ac" = beta_female_ac,  "alpha_female_ac" = alpha_female_ac,  "delta_female_ac" = delta_female_ac,  "gamma_female_ac" = gamma_female_ac, "ductal_female_ac" = ductal_female_ac, "acinar_female_ac" = acinar_female_ac, "astell_female_ac" = astell_female_ac, "macro_female_ac"=macro_female_ac, "endo_female_ac" = endo_female_ac)
gene.list <- unique(allgenes)
allDEandDAR <- intersect(unique.de.genes.all, gene.list)

# Go back up and plot heatmap note that this is RNAseq data
dittoHeatmap(
  nd.combined,
  genes = allDEandDAR, #genes.to.plot.intr,
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
  use_raster = FALSE,
  raster_quality = 5
)


# Motif analysis an plotting
# For accessibility remember to run the code again generating the _beta files and subset 100kb window peaks
allregions <- c(
  as.character(cgenes_beta$query_region),
  as.character(cgenes_delta$query_region),
  as.character(cgenes_alpha$query_region),
  as.character(cgenes_gamma$query_region),
  as.character(cgenes_ductal$query_region),
  as.character(cgenes_acinar$query_region),
  as.character(cgenes_astel$query_region),
  as.character(cgenes_qstel$query_region),
  as.character(cgenes_endo$query_region),
  as.character(cgenes_lympho$query_region),
  as.character(cgenes_macro$query_region)
)

uniqueregions <- unique(allregions)

peaks_beta <- as.character(cgenes_beta$query_region)
peaks_delta <- as.character(cgenes_delta$query_region)
peaks_alpha <- as.character(cgenes_alpha$query_region)
peaks_gamma <- as.character(cgenes_gamma$query_region)
peaks_ductal <- as.character(cgenes_ductal$query_region)
peaks_acinar <- as.character(cgenes_acinar$query_region)
peaks_astel <- as.character(cgenes_astel$query_region)
peaks_qstel <- as.character(cgenes_qstel$query_region)
peaks_endo <- as.character(cgenes_endo$query_region)
peaks_lympho <- as.character(cgenes_lympho$query_region)
peaks_macro <- as.character(cgenes_macro$query_region)

# Motif analysis
Idents(hm.integrated.dfree) <- "celltype"
open.peaks <- AccessiblePeaks(hm.integrated.dfree, idents = unique(as.character(hm.integrated.dfree$celltype)))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(hm.integrated.dfree, assay = "macs2", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[uniqueregions, ],
  n = 50000
)

# Motif tests
# This is not super accurate its better to run motif analysis using chromvar
{
enriched.motifs.beta <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_beta)
enriched.motifs.delta <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_delta)
enriched.motifs.alpha <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_alpha)
enriched.motifs.gamma <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_gamma)
enriched.motifs.ductal <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_ductal)
enriched.motifs.acinar <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_acinar)
enriched.motifs.astel <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_astel)
enriched.motifs.qstel <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_qstel)
enriched.motifs.endo <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_endo)
enriched.motifs.lympho <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_lympho)
enriched.motifs.macro <- FindMotifs(object = hm.integrated.dfree, background=peaks.matched, features = peaks_macro)

# write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.beta.csv)")
# write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.delta.csv)")
# write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.alpha.csv)")
# write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.gamma.csv)")
# write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.ductal.csv)")
# write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.acinar.csv)")
# write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.astel.csv)")
# write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.qstel.csv)")
# write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.endo.csv)")
# write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.lympho.csv)")
# write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\motifs\allsex\enriched.motifs.macro.csv)")

write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.beta.csv)")
write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.delta.csv)")
write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.alpha.csv)")
write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.gamma.csv)")
write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.ductal.csv)")
write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.acinar.csv)")
write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.astel.csv)")
write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.qstel.csv)")
write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.endo.csv)")
write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.lympho.csv)")
write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\Motif_ave\enriched.motifs.macro.csv)")
}


# Add Motif scores to pseduo bulk and make a heat map
# Heatmap
# Pseudobulk
DefaultAssay(hm.integrated.dfree) <- "chromvar_macs2"
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry_lib"
combined_processed_atac <- AverageExpression(hm.integrated.dfree, 
                                             assays = c("chromvar_macs2", "macs2", "RNA"), 
                                             features = NULL, return.seurat = TRUE,  
                                             group.by = "celltype_sex_ancestry_lib",
                                             slot = "data", verbose = FALSE)

combined_processed_atac$celltype_sex_ancestry_lib <- Cells(combined_processed_atac) #6892 Seurat

{
  Idents(combined_processed_atac) <- 'celltype_sex_ancestry_lib'
  combined_processed_atac$celltype <- combined_processed_atac@meta.data[["orig.ident"]]
  metadat <- combined_processed_atac@meta.data
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "activated_stellate", "activated-stellate"))
  metadat <- metadat %>% 
    mutate(celltype_sex_ancestry_lib = str_replace(celltype_sex_ancestry_lib, "quiescent_stellate", "quiescent-stellate"))
  metadat$sex <- metadat[c('sex')] <- str_split_i(metadat$celltype_sex_ancestry_lib, "_", -3)
  metadat$ancestry <- metadat[c('ancestry')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -2)
  metadat$lib <- metadat[c('lib')] <- str_split_i(metadat$celltype_sex_ancestry_lib, '_', -1)
  combined_processed_atac@meta.data = metadat
}


table(combined_processed_atac@meta.data[["celltype"]])
table(combined_processed_atac@meta.data[["sex"]])

# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "delta", "alpha", "gamma",
               "ductal", "acinar",
               "activated", "quiescent", "endothelial",
               "lymphocyte", "macrophage")

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
combined_processed_atac$celltype <- factor(x = combined_processed_atac$celltype, levels = my_levels)
table(combined_processed_atac$celltype)
Idents(combined_processed_atac) <- "celltype"


# DE testing Motifs
DefaultAssay(hm.integrated.dfree) <- "chromvar_macs2"
Idents(hm.integrated.dfree) <- "celltype"

enriched.cvar.motifs.beta <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'beta',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.alpha <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'alpha',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.delta <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'delta',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.gamma <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'gamma',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.astel <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'activated_stellate',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.qstel <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'quiescent_stellate',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.ductal <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'ductal',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.acinar <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'acinar',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.macro <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'macrophage',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.lympho <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'lymphocyte',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.endo <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'endothelial',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.cvar.motifs.beta
enriched.cvar.motifs.alpha
enriched.cvar.motifs.delta
enriched.cvar.motifs.gamma
enriched.cvar.motifs.astel
enriched.cvar.motifs.qstel
enriched.cvar.motifs.ductal
enriched.cvar.motifs.acinar
enriched.cvar.motifs.macro
enriched.cvar.motifs.lympho
enriched.cvar.motifs.endo

# Adjust for incalculable pvals
enriched.cvar.motifs.beta$p_val_adj[enriched.cvar.motifs.beta$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.alpha$p_val_adj[enriched.cvar.motifs.alpha$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.delta$p_val_adj[enriched.cvar.motifs.delta$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.gamma$p_val_adj[enriched.cvar.motifs.gamma$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.astel$p_val_adj[enriched.cvar.motifs.astel$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.qstel$p_val_adj[enriched.cvar.motifs.qstel$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.ductal$p_val_adj[enriched.cvar.motifs.ductal$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.acinar$p_val_adj[enriched.cvar.motifs.acinar$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.macro$p_val_adj[enriched.cvar.motifs.macro$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.lympho$p_val_adj[enriched.cvar.motifs.lympho$p_val_adj == 0] <- 2e-307
enriched.cvar.motifs.endo$p_val_adj[enriched.cvar.motifs.endo$p_val_adj == 0] <- 2e-307

# Translating Motifs
motif_id <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\motif.ID.csv)", row.names = 1)
enriched.motifs.beta <- merge(motif_id, enriched.cvar.motifs.beta, by = 'row.names')
enriched.motifs.delta <- merge(motif_id, enriched.cvar.motifs.delta, by = 'row.names')
enriched.motifs.alpha <- merge(motif_id, enriched.cvar.motifs.alpha, by = 'row.names')
enriched.motifs.gamma <- merge(motif_id, enriched.cvar.motifs.gamma, by = 'row.names')
enriched.motifs.astel <- merge(motif_id, enriched.cvar.motifs.astel, by = 'row.names')
enriched.motifs.qstel <- merge(motif_id, enriched.cvar.motifs.qstel, by = 'row.names')
enriched.motifs.ductal <- merge(motif_id, enriched.cvar.motifs.ductal, by = 'row.names')
enriched.motifs.acinar <- merge(motif_id, enriched.cvar.motifs.acinar, by = 'row.names')
enriched.motifs.macro <- merge(motif_id, enriched.cvar.motifs.macro, by = 'row.names')
enriched.motifs.lympho <- merge(motif_id, enriched.cvar.motifs.lympho, by = 'row.names')
enriched.motifs.endo <- merge(motif_id, enriched.cvar.motifs.endo, by = 'row.names')

enriched.motifs.beta <- enriched.motifs.beta %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.delta <- enriched.motifs.delta %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.alpha <- enriched.motifs.alpha %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.gamma <- enriched.motifs.gamma %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.astel <- enriched.motifs.astel %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.qstel <- enriched.motifs.qstel %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.ductal <- enriched.motifs.ductal %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.acinar <- enriched.motifs.acinar %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.macro <- enriched.motifs.macro %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.lympho <- enriched.motifs.lympho %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.endo <- enriched.motifs.endo %>% remove_rownames %>% column_to_rownames(var="tf")

# write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.beta.csv)")
# write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.alpha.csv)")
# write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.delta.csv)")
# write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.gamma.csv)")
# write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.activated.stellate.csv)")
# write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.quiescent.stellate.csv)")
# write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.ductal.csv)")
# write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.acinar.csv)")
# write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.macro.csv)")
# write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.lympho.csv)")
# write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.endo.csv)")

# Load in motifs From 
#Manually add ordering
enriched.motifs.beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.beta.csv)", row.names = 1)
enriched.motifs.alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.alpha.csv)", row.names = 1)
enriched.motifs.delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.delta.csv)", row.names = 1)
enriched.motifs.gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.gamma.csv)", row.names = 1)
enriched.motifs.ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.ductal.csv)", row.names = 1)
enriched.motifs.acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.acinar.csv)", row.names = 1)
enriched.motifs.astel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.activated.stellate.csv)", row.names = 1)
enriched.motifs.qstel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.quiescent.stellate.csv)", row.names = 1)
enriched.motifs.endo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.endo.csv)", row.names = 1)
enriched.motifs.lympho <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.lympho.csv)", row.names = 1)
enriched.motifs.macro <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.macro.csv)", row.names = 1)

# Order data based on pval #beta_motifs <- as.character(rownames(enriched.motifs.beta[enriched.motifs.beta$p_val_adj < 1e-50,]))
enriched.motifs.beta <- enriched.motifs.beta[order(enriched.motifs.beta$ordering),] 
enriched.motifs.delta <- enriched.motifs.delta[order(enriched.motifs.delta$ordering),] 
enriched.motifs.alpha <- enriched.motifs.alpha[order(enriched.motifs.alpha$ordering),] 
enriched.motifs.gamma <- enriched.motifs.gamma[order(enriched.motifs.gamma$ordering),]
enriched.motifs.astel <- enriched.motifs.astel[order(enriched.motifs.astel$ordering),] 
enriched.motifs.qstel <- enriched.motifs.qstel[order(enriched.motifs.qstel$ordering),]
enriched.motifs.ductal <- enriched.motifs.ductal[order(enriched.motifs.ductal$ordering),] 
enriched.motifs.acinar <- enriched.motifs.acinar[order(enriched.motifs.acinar$ordering),]
enriched.motifs.macro <- enriched.motifs.macro[order(enriched.motifs.macro$ordering),] 
enriched.motifs.lympho <- enriched.motifs.lympho[order(enriched.motifs.lympho$ordering),]
enriched.motifs.endo <- enriched.motifs.endo[order(enriched.motifs.endo$ordering),]

# Add divided data
enriched.motifs.beta$foldenrich <- enriched.motifs.beta$pct.1 / enriched.motifs.beta$pct.2
enriched.motifs.delta$foldenrich <- enriched.motifs.delta$pct.1 / enriched.motifs.delta$pct.2
enriched.motifs.alpha$foldenrich <- enriched.motifs.alpha$pct.1 / enriched.motifs.alpha$pct.2
enriched.motifs.gamma$foldenrich <- enriched.motifs.gamma$pct.1 / enriched.motifs.gamma$pct.2
enriched.motifs.astel$foldenrich <- enriched.motifs.astel$pct.1 / enriched.motifs.astel$pct.2
enriched.motifs.qstel$foldenrich <- enriched.motifs.qstel$pct.1 / enriched.motifs.qstel$pct.2
enriched.motifs.ductal$foldenrich <- enriched.motifs.ductal$pct.1 / enriched.motifs.ductal$pct.2
enriched.motifs.acinar$foldenrich <- enriched.motifs.acinar$pct.1 / enriched.motifs.acinar$pct.2
enriched.motifs.macro$foldenrich <- enriched.motifs.macro$pct.1 / enriched.motifs.macro$pct.2
enriched.motifs.lympho$foldenrich <- enriched.motifs.lympho$pct.1 / enriched.motifs.lympho$pct.2
enriched.motifs.endo$foldenrich <- enriched.motifs.endo$pct.1 / enriched.motifs.endo$pct.2

# Plotting heatmap from ->
enriched.motifs.beta <- dplyr::filter(enriched.motifs.beta, p_val_adj < 5e-2) 
enriched.motifs.delta <- dplyr::filter(enriched.motifs.delta, p_val_adj < 5e-2) 
enriched.motifs.alpha <- dplyr::filter(enriched.motifs.alpha, p_val_adj < 5e-2) 
enriched.motifs.gamma <- dplyr::filter(enriched.motifs.gamma, p_val_adj < 5e-2) 
enriched.motifs.astel <- dplyr::filter(enriched.motifs.astel, p_val_adj < 5e-2) 
enriched.motifs.qstel <- dplyr::filter(enriched.motifs.qstel, p_val_adj < 5e-2) 
enriched.motifs.ductal <- dplyr::filter(enriched.motifs.ductal, p_val_adj < 5e-2) 
enriched.motifs.acinar <- dplyr::filter(enriched.motifs.acinar, p_val_adj < 5e-2) 
enriched.motifs.macro <- dplyr::filter(enriched.motifs.macro, p_val_adj < 5e-2) 
enriched.motifs.lympho <- dplyr::filter(enriched.motifs.lympho, p_val_adj < 5e-2)
enriched.motifs.endo <- dplyr::filter(enriched.motifs.endo, p_val_adj < 5e-2)

# Select top 50 motifs
enriched.motifs.beta <- top_n(enriched.motifs.beta, -50, ordering)
enriched.motifs.delta <- top_n(enriched.motifs.delta, -50, ordering)
enriched.motifs.alpha <- top_n(enriched.motifs.alpha, -50, ordering)
enriched.motifs.gamma <- top_n(enriched.motifs.gamma, -50, ordering)
enriched.motifs.astel <- top_n(enriched.motifs.astel, -50, ordering)
enriched.motifs.qstel <- top_n(enriched.motifs.qstel, -50, ordering)
enriched.motifs.ductal <- top_n(enriched.motifs.ductal, -50, ordering)
enriched.motifs.acinar <- top_n(enriched.motifs.acinar, -50, ordering)
enriched.motifs.macro <- top_n(enriched.motifs.macro, -50, ordering)
enriched.motifs.lympho <- top_n(enriched.motifs.lympho, -50, ordering)
enriched.motifs.endo <- top_n(enriched.motifs.endo, -50, ordering)

# Order data #beta_motifs <- as.character(rownames(enriched.motifs.beta[enriched.motifs.beta$p_val_adj < 1e-50,]))
enriched.motifs.beta <- enriched.motifs.beta[order(enriched.motifs.beta$foldenrich),] 
enriched.motifs.delta <- enriched.motifs.delta[order(enriched.motifs.delta$foldenrich),] 
enriched.motifs.alpha <- enriched.motifs.alpha[order(enriched.motifs.alpha$foldenrich),] 
enriched.motifs.gamma <- enriched.motifs.gamma[order(enriched.motifs.gamma$foldenrich),]
enriched.motifs.astel <- enriched.motifs.astel[order(enriched.motifs.astel$foldenrich),] 
enriched.motifs.qstel <- enriched.motifs.qstel[order(enriched.motifs.qstel$foldenrich),]
enriched.motifs.ductal <- enriched.motifs.ductal[order(enriched.motifs.ductal$foldenrich),] 
enriched.motifs.acinar <- enriched.motifs.acinar[order(enriched.motifs.acinar$foldenrich),]
enriched.motifs.macro <- enriched.motifs.macro[order(enriched.motifs.macro$foldenrich),] 
enriched.motifs.lympho <- enriched.motifs.lympho[order(enriched.motifs.lympho$foldenrich),]
enriched.motifs.endo <- enriched.motifs.endo[order(enriched.motifs.endo$foldenrich),]

# -log10adjpval
enriched.motifs.beta$neglogpval <- -log10(enriched.motifs.beta$p_val_adj)
enriched.motifs.delta$neglogpval <- -log10(enriched.motifs.delta$p_val_adj)
enriched.motifs.alpha$neglogpval <- -log10(enriched.motifs.alpha$p_val_adj)
enriched.motifs.gamma$neglogpval <- -log10(enriched.motifs.gamma$p_val_adj)
enriched.motifs.astel$neglogpval <- -log10(enriched.motifs.astel$p_val_adj)
enriched.motifs.qstel$neglogpval <- -log10(enriched.motifs.qstel$p_val_adj)
enriched.motifs.ductal$neglogpval <- -log10(enriched.motifs.ductal$p_val_adj)
enriched.motifs.acinar$neglogpval <- -log10(enriched.motifs.acinar$p_val_adj)
enriched.motifs.macro$neglogpval <- -log10(enriched.motifs.macro$p_val_adj)
enriched.motifs.lympho$neglogpval <- -log10(enriched.motifs.lympho$p_val_adj)
enriched.motifs.endo$neglogpval <- -log10(enriched.motifs.endo$p_val_adj)

# Plot
plot_celldata <- top_n(enriched.motifs.macro, -50, foldenrich)
plot_celldata <- tibble::rownames_to_column(plot_celldata, "motif.name")
level_order <- plot_celldata$motif.name
mid <- mean(plot_celldata$foldenrich)
plot_celldata %>%
  arrange(foldenrich) %>%
  mutate(name=factor(motif.name, levels=motif.name)) %>%
  ggplot(aes(x = factor(motif.name, level = level_order),
             y = foldenrich, 
             color = neglogpval, 
             size = pct.1)) + 
  #geom_point(alpha=0.5, colour = "springgreen4") + #use when all neglog pvals are same
  geom_point(alpha=0.5) + #dont geom twice
  ylab("") + 
  xlab("") + 
  ggtitle("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3, hjust=1, size =6, face = "plain", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =6, face = "plain", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "plain"),
        legend.title=element_text(size=8, face = "plain"), 
        legend.text=element_text(size=8, face = "plain")) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  coord_flip() +
  scale_size_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)) + #only till here when neglogs are the same for all
  #scale_x_discrete(limits=rev) + #reverses order 
  scale_color_gradient(low = "white", high="magenta3", space ="Lab")

## Motif analysis across sex
# DE testing Motifs
DefaultAssay(hm.integrated.dfree) <- "chromvar_macs2"
colnames(hm.integrated.dfree@meta.data)
Idents(hm.integrated.dfree) <- "celltype_sex"
unique(hm.integrated.dfree@active.ident)

enriched.sex.motifs.beta <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'beta_male',
  ident.2 = 'beta_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.alpha <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'alpha_male',
  ident.2 = 'alpha_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.delta <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'delta_male',
  ident.2 = 'delta_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.gamma <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'gamma_male',
  ident.2 = 'gamma_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.astel <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'activated_stellate_male',
  ident.2 = 'activated_stellate_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.qstel <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'quiescent_stellate_male',
  ident.2 = 'quiescent_stellate_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.ductal <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'ductal_male',
  ident.2 = 'ductal_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.acinar <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'acinar_male',
  ident.2 = 'acinar_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.macro <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'macrophage_male',
  ident.2 = 'macrophage_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.lympho <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'lymphocyte_male',
  ident.2 = 'lymphocyte_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.endo <- FindMarkers(
  object = hm.integrated.dfree,
  ident.1 = 'endothelial_male',
  ident.2 = 'endothelial_female',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  min.pct = 0.1,
  logfc.threshold = 0
)

enriched.sex.motifs.beta
enriched.sex.motifs.alpha
enriched.sex.motifs.delta
enriched.sex.motifs.gamma
enriched.sex.motifs.astel
enriched.sex.motifs.qstel
enriched.sex.motifs.ductal
enriched.sex.motifs.acinar
enriched.sex.motifs.macro
enriched.sex.motifs.lympho
enriched.sex.motifs.endo

# Adjust for incalculable pvals
enriched.sex.motifs.beta$p_val_adj[enriched.sex.motifs.beta$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.alpha$p_val_adj[enriched.sex.motifs.alpha$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.delta$p_val_adj[enriched.sex.motifs.delta$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.gamma$p_val_adj[enriched.sex.motifs.gamma$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.astel$p_val_adj[enriched.sex.motifs.astel$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.qstel$p_val_adj[enriched.sex.motifs.qstel$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.ductal$p_val_adj[enriched.sex.motifs.ductal$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.acinar$p_val_adj[enriched.sex.motifs.acinar$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.macro$p_val_adj[enriched.sex.motifs.macro$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.lympho$p_val_adj[enriched.sex.motifs.lympho$p_val_adj == 0] <- 2e-307
enriched.sex.motifs.endo$p_val_adj[enriched.sex.motifs.endo$p_val_adj == 0] <- 2e-307

# Translating Motifs
motif_id <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\motif.ID.csv)", row.names = 1)
enriched.motifs.beta <- merge(motif_id, enriched.sex.motifs.beta, by = 'row.names')
enriched.motifs.delta <- merge(motif_id, enriched.sex.motifs.delta, by = 'row.names')
enriched.motifs.alpha <- merge(motif_id, enriched.sex.motifs.alpha, by = 'row.names')
enriched.motifs.gamma <- merge(motif_id, enriched.sex.motifs.gamma, by = 'row.names')
enriched.motifs.astel <- merge(motif_id, enriched.sex.motifs.astel, by = 'row.names')
enriched.motifs.qstel <- merge(motif_id, enriched.sex.motifs.qstel, by = 'row.names')
enriched.motifs.ductal <- merge(motif_id, enriched.sex.motifs.ductal, by = 'row.names')
enriched.motifs.acinar <- merge(motif_id, enriched.sex.motifs.acinar, by = 'row.names')
enriched.motifs.macro <- merge(motif_id, enriched.sex.motifs.macro, by = 'row.names')
enriched.motifs.lympho <- merge(motif_id, enriched.sex.motifs.lympho, by = 'row.names')
enriched.motifs.endo <- merge(motif_id, enriched.sex.motifs.endo, by = 'row.names')

enriched.motifs.beta <- enriched.motifs.beta %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.delta <- enriched.motifs.delta %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.alpha <- enriched.motifs.alpha %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.gamma <- enriched.motifs.gamma %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.astel <- enriched.motifs.astel %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.qstel <- enriched.motifs.qstel %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.ductal <- enriched.motifs.ductal %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.acinar <- enriched.motifs.acinar %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.macro <- enriched.motifs.macro %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.lympho <- enriched.motifs.lympho %>% remove_rownames %>% column_to_rownames(var="tf")
enriched.motifs.endo <- enriched.motifs.endo %>% remove_rownames %>% column_to_rownames(var="tf")

# write.csv(enriched.motifs.beta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.beta.csv)")
# write.csv(enriched.motifs.alpha, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.alpha.csv)")
# write.csv(enriched.motifs.delta, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.delta.csv)")
# write.csv(enriched.motifs.gamma, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.gamma.csv)")
# write.csv(enriched.motifs.astel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.activated.stellate.csv)")
# write.csv(enriched.motifs.qstel, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.quiescent.stellate.csv)")
# write.csv(enriched.motifs.ductal, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.ductal.csv)")
# write.csv(enriched.motifs.acinar, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.acinar.csv)")
# write.csv(enriched.motifs.macro, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.macro.csv)")
# write.csv(enriched.motifs.lympho, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.lympho.csv)")
# write.csv(enriched.motifs.endo, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.endo.csv)")

# Load in motifs From 
#Manually add ordering
enriched.motifs.beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.beta.csv)", row.names = 1)
enriched.motifs.alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.alpha.csv)", row.names = 1)
enriched.motifs.delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.delta.csv)", row.names = 1)
enriched.motifs.gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.gamma.csv)", row.names = 1)
enriched.motifs.ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.ductal.csv)", row.names = 1)
enriched.motifs.acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.acinar.csv)", row.names = 1)
enriched.motifs.astel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.activated.stellate.csv)", row.names = 1)
enriched.motifs.qstel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.quiescent.stellate.csv)", row.names = 1)
enriched.motifs.endo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.endo.csv)", row.names = 1)
enriched.motifs.lympho <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.lympho.csv)", row.names = 1)
enriched.motifs.macro <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.macro.csv)", row.names = 1)

# # Order data based on pval #beta_motifs <- as.character(rownames(enriched.motifs.beta[enriched.motifs.beta$p_val_adj < 1e-50,]))
# enriched.motifs.beta <- enriched.motifs.beta[order(enriched.motifs.beta$ordering),] 
# enriched.motifs.delta <- enriched.motifs.delta[order(enriched.motifs.delta$ordering),] 
# enriched.motifs.alpha <- enriched.motifs.alpha[order(enriched.motifs.alpha$ordering),] 
# enriched.motifs.gamma <- enriched.motifs.gamma[order(enriched.motifs.gamma$ordering),]
# enriched.motifs.astel <- enriched.motifs.astel[order(enriched.motifs.astel$ordering),] 
# enriched.motifs.qstel <- enriched.motifs.qstel[order(enriched.motifs.qstel$ordering),]
# enriched.motifs.ductal <- enriched.motifs.ductal[order(enriched.motifs.ductal$ordering),] 
# enriched.motifs.acinar <- enriched.motifs.acinar[order(enriched.motifs.acinar$ordering),]
# enriched.motifs.macro <- enriched.motifs.macro[order(enriched.motifs.macro$ordering),] 
# enriched.motifs.lympho <- enriched.motifs.lympho[order(enriched.motifs.lympho$ordering),]
# enriched.motifs.endo <- enriched.motifs.endo[order(enriched.motifs.endo$ordering),]

# Add divided data
enriched.motifs.beta$foldenrich <- enriched.motifs.beta$pct.1 / enriched.motifs.beta$pct.2
enriched.motifs.delta$foldenrich <- enriched.motifs.delta$pct.1 / enriched.motifs.delta$pct.2
enriched.motifs.alpha$foldenrich <- enriched.motifs.alpha$pct.1 / enriched.motifs.alpha$pct.2
enriched.motifs.gamma$foldenrich <- enriched.motifs.gamma$pct.1 / enriched.motifs.gamma$pct.2
enriched.motifs.astel$foldenrich <- enriched.motifs.astel$pct.1 / enriched.motifs.astel$pct.2
enriched.motifs.qstel$foldenrich <- enriched.motifs.qstel$pct.1 / enriched.motifs.qstel$pct.2
enriched.motifs.ductal$foldenrich <- enriched.motifs.ductal$pct.1 / enriched.motifs.ductal$pct.2
enriched.motifs.acinar$foldenrich <- enriched.motifs.acinar$pct.1 / enriched.motifs.acinar$pct.2
enriched.motifs.macro$foldenrich <- enriched.motifs.macro$pct.1 / enriched.motifs.macro$pct.2
enriched.motifs.lympho$foldenrich <- enriched.motifs.lympho$pct.1 / enriched.motifs.lympho$pct.2
enriched.motifs.endo$foldenrich <- enriched.motifs.endo$pct.1 / enriched.motifs.endo$pct.2

# # Lets merge into 1
# MyMerge<- function(x, y){
#   df<- merge(x, y, by= "row.names", all.x= F, all.y= F)
#   rownames(df)<- df$Row.names
#   df$Row.names<- NULL
#   return(df)
# }
# mergedmotifs <- Reduce(MyMerge, list(enriched.motifs.beta, enriched.motifs.delta, enriched.motifs.alpha, enriched.motifs.gamma, enriched.motifs.astel, enriched.motifs.qstel, 
#                                       enriched.motifs.ductal, enriched.motifs.acinar, enriched.motifs.macro, enriched.motifs.lympho, enriched.motifs.endo))
# 
# mergedmotifs <- do.call("rbind", list(enriched.motifs.beta, enriched.motifs.delta, enriched.motifs.alpha, enriched.motifs.gamma, enriched.motifs.astel, enriched.motifs.qstel, 
#                                       enriched.motifs.ductal, enriched.motifs.acinar, enriched.motifs.macro, enriched.motifs.lympho, enriched.motifs.endo))

# Plotting heatmap from ->
# enriched.motifs.beta <- dplyr::filter(enriched.motifs.beta, p_val_adj < 5e-2) 
# enriched.motifs.delta <- dplyr::filter(enriched.motifs.delta, p_val_adj < 5e-2) 
# enriched.motifs.alpha <- dplyr::filter(enriched.motifs.alpha, p_val_adj < 5e-2) 
# enriched.motifs.gamma <- dplyr::filter(enriched.motifs.gamma, p_val_adj < 5e-2) 
# enriched.motifs.astel <- dplyr::filter(enriched.motifs.astel, p_val_adj < 5e-2) 
# enriched.motifs.qstel <- dplyr::filter(enriched.motifs.qstel, p_val_adj < 5e-2) 
# enriched.motifs.ductal <- dplyr::filter(enriched.motifs.ductal, p_val_adj < 5e-2) 
# enriched.motifs.acinar <- dplyr::filter(enriched.motifs.acinar, p_val_adj < 5e-2) 
# enriched.motifs.macro <- dplyr::filter(enriched.motifs.macro, p_val_adj < 5e-2) 
# enriched.motifs.lympho <- dplyr::filter(enriched.motifs.lympho, p_val_adj < 5e-2)
# enriched.motifs.endo <- dplyr::filter(enriched.motifs.endo, p_val_adj < 5e-2)

# Plot violins
volcanodat <- enriched.motifs.beta

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(volcanodat))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(volcanodat))

# modify keyvals for variables with fold change > 1 note for T2D padj < 0.01 #keyvals[which(volcanodat$log2FoldChange > 1 & volcanodat$padj < 0.01)]
keyvals[which(volcanodat$foldenrich > 1 & volcanodat$p_val_adj < 0.05)] <- 'royalblue'
names(keyvals)[which(volcanodat$foldenrich > 1 & volcanodat$p_val_adj < 0.05)] <- 'high'

# modify keyvals for variables with fold change < 1
keyvals[which(volcanodat$foldenrich < 1 & volcanodat$p_val_adj < 0.05)] <- 'red'
names(keyvals)[which(volcanodat$foldenrich < 1 & volcanodat$p_val_adj < 0.05)] <- 'low'

unique(names(keyvals))

unique(keyvals)
keyvals[1:20]
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(volcanodat,
                lab = rownames(volcanodat),
                x = 'foldenrich',
                y = 'p_val_adj',
                #selectLab = TRUE,
                #selectLab = rownames(volcanodat)[which(names(keyvals) %in% c('high', 'low'))],
                # selectLab = c("NFYA", "NFYB", "NFYC", "ELK3", "ETS2",
                #               "ETV3", "THAP11", "FEV", "FLI1", "ELK1", "Klf12", 
                #               "ETS1", "CTCF", "ERG", "ERF", "ETV2", "Ddit3::Cebpa", 
                #               "DBP", "ZNF143", "RFX2",
                #               
                #               "IRF1", "JUND", "JUNB", "FOSL2::JUNB",
                #               "EGR1", "GATA5", "FOS", "KLF4", "Sox3",
                #               "JUN(var.2)", "FOS::JUN", "FOSL1::JUNB", "FOSL2::JUND", "FOSB::JUNB",
                #               "FOS::JUNB", "BATF::JUN", "FOSL1", "BATF3", "BATF", "FOSL2::JUN" #beta
                selectLab = c("FOSL2", "FOS::JUN", "Smad2::Smad3", "BATF3", "FOSL1::JUND",
                              "FOS::JUND", "BATF", "JUN(var.2)", "FOS", "FOSL1::JUNB", 
                              "FOSL1::JUN", "FOS::JUNB", "NFE2", "FOSL2::JUN", "FOSL1", 
                              "FOSB::JUNB", "BATF::JUN", "BACH2", "FOSL2::JUND", "JUNB",

                              "Ptf1a", "TAL1::TCF3", "HAND2", "Wt1", "Arid3a", 
                              "Ascl2", "EWSR1-FLI1", "Tcf12", "NEUROD1", "Tcf21", 
                              "ZBTB18", "MYOG", "PRDM1", "BHLHA15(var.2)", "NEUROG2(var.2)",
                              "EGR1", "KLF15", "ZNF263", "ATOH7", "ATOH1(var.2)" #alpha
                ), # use this for labelling genes on plot
                #encircle = c('VAMP3'),
                boxedLabels = FALSE,
                #xlim = c(0.1,3),
                #ylim = c(0,30),
                xlab = bquote('Fold Enrichment'),
                title = 'Custom colour over-ride',
                pCutoff = 0.05,
                FCcutoff = c(1), # you can make multiple lines
                pointSize = c(ifelse((volcanodat$foldenrich > 1 & volcanodat$p_val_adj < 0.05) | (volcanodat$foldenrich < 1 & volcanodat$p_val_adj < 0.05), 3, 2)), #changed from T2D 0.01
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


# Plot difference in motif accessibility for select motifs across sex
DefaultAssay(hm.integrated.dfree) <- "chromvar_macs2"
Idents(hm.integrated.dfree) <- "celltype_sex"

#Subset out only alpha and beta cells
alpha_beta <- subset(hm.integrated.dfree, idents = c("beta_female", "beta_male", "alpha_female", "alpha_male"))

VlnPlot(
  alpha_beta,
  features = "MA0007.3", #Ar
  pt.size = 0)

DotPlot(alpha_beta,  
        dot.scale = 20,
        col.min = -1, #minimum level
        col.max = 1,  #maximum level
        features = c("MA0007.3"),
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
  guides(color = guide_colorbar(title = 'Average Expression')) +
  scale_size_continuous(limits = c(10, 60)) +  # specify your desired range
  coord_flip() 

DotPlot(hm.integrated.dfree,  
        dot.scale = 20,
        col.min = -1, #minimum level
        col.max = 1,  #maximum level
        features = c("chr19-19627168-19629130"),
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
  guides(color = guide_colorbar(title = 'Average Expression')) +
  scale_size_continuous(limits = c(10, 60)) +  # specify your desired range
  coord_flip() 

# Plot DE motifs
# Select top 50 motifs
enriched.motifs.beta <- top_n(enriched.motifs.beta, -50, ordering)
enriched.motifs.delta <- top_n(enriched.motifs.delta, -50, ordering)
enriched.motifs.alpha <- top_n(enriched.motifs.alpha, -50, ordering)
enriched.motifs.gamma <- top_n(enriched.motifs.gamma, -50, ordering)
enriched.motifs.astel <- top_n(enriched.motifs.astel, -50, ordering)
enriched.motifs.qstel <- top_n(enriched.motifs.qstel, -50, ordering)
enriched.motifs.ductal <- top_n(enriched.motifs.ductal, -50, ordering)
enriched.motifs.acinar <- top_n(enriched.motifs.acinar, -50, ordering)
enriched.motifs.macro <- top_n(enriched.motifs.macro, -50, ordering)
enriched.motifs.lympho <- top_n(enriched.motifs.lympho, -50, ordering)
enriched.motifs.endo <- top_n(enriched.motifs.endo, -50, ordering)

# Order data #beta_motifs <- as.character(rownames(enriched.motifs.beta[enriched.motifs.beta$p_val_adj < 1e-50,]))
enriched.motifs.beta <- enriched.motifs.beta[order(enriched.motifs.beta$foldenrich),] 
enriched.motifs.delta <- enriched.motifs.delta[order(enriched.motifs.delta$foldenrich),] 
enriched.motifs.alpha <- enriched.motifs.alpha[order(enriched.motifs.alpha$foldenrich),] 
enriched.motifs.gamma <- enriched.motifs.gamma[order(enriched.motifs.gamma$foldenrich),]
enriched.motifs.astel <- enriched.motifs.astel[order(enriched.motifs.astel$foldenrich),] 
enriched.motifs.qstel <- enriched.motifs.qstel[order(enriched.motifs.qstel$foldenrich),]
enriched.motifs.ductal <- enriched.motifs.ductal[order(enriched.motifs.ductal$foldenrich),] 
enriched.motifs.acinar <- enriched.motifs.acinar[order(enriched.motifs.acinar$foldenrich),]
enriched.motifs.macro <- enriched.motifs.macro[order(enriched.motifs.macro$foldenrich),] 
enriched.motifs.lympho <- enriched.motifs.lympho[order(enriched.motifs.lympho$foldenrich),]
enriched.motifs.endo <- enriched.motifs.endo[order(enriched.motifs.endo$foldenrich),]

# -log10adjpval
enriched.motifs.beta$neglogpval <- -log10(enriched.motifs.beta$p_val_adj)
enriched.motifs.delta$neglogpval <- -log10(enriched.motifs.delta$p_val_adj)
enriched.motifs.alpha$neglogpval <- -log10(enriched.motifs.alpha$p_val_adj)
enriched.motifs.gamma$neglogpval <- -log10(enriched.motifs.gamma$p_val_adj)
enriched.motifs.astel$neglogpval <- -log10(enriched.motifs.astel$p_val_adj)
enriched.motifs.qstel$neglogpval <- -log10(enriched.motifs.qstel$p_val_adj)
enriched.motifs.ductal$neglogpval <- -log10(enriched.motifs.ductal$p_val_adj)
enriched.motifs.acinar$neglogpval <- -log10(enriched.motifs.acinar$p_val_adj)
enriched.motifs.macro$neglogpval <- -log10(enriched.motifs.macro$p_val_adj)
enriched.motifs.lympho$neglogpval <- -log10(enriched.motifs.lympho$p_val_adj)
enriched.motifs.endo$neglogpval <- -log10(enriched.motifs.endo$p_val_adj)

# Plot
plot_celldata <- top_n(enriched.motifs.macro, -50, foldenrich)
plot_celldata <- tibble::rownames_to_column(plot_celldata, "motif.name")
level_order <- plot_celldata$motif.name
mid <- mean(plot_celldata$foldenrich)
plot_celldata %>%
  arrange(foldenrich) %>%
  mutate(name=factor(motif.name, levels=motif.name)) %>%
  ggplot(aes(x = factor(motif.name, level = level_order),
             y = foldenrich, 
             color = neglogpval, 
             size = pct.1)) + 
  #geom_point(alpha=0.5, colour = "springgreen4") + #use when all neglog pvals are same
  geom_point(alpha=0.5) + #dont geom twice
  ylab("") + 
  xlab("") + 
  ggtitle("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.3, hjust=1, size =6, face = "plain", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =6, face = "plain", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "plain"),
        legend.title=element_text(size=8, face = "plain"), 
        legend.text=element_text(size=8, face = "plain")) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  coord_flip() +
  scale_size_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)) + #only till here when neglogs are the same for all
  #scale_x_discrete(limits=rev) + #reverses order 
  scale_color_gradient(low = "white", high="magenta3", space ="Lab")

# TF Foot printing
# gather the footprinting information for sets of motifs
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")
DefaultAssay(hm.integrated.dfree) <- "macs2"
Idents(hm.integrated.dfree) <- "celltype"
hm.integrated.dfree <- Footprint(
  object = hm.integrated.dfree,
  motif.name = c("PDX1", "MAFA", "Arx", "NKX6-1", "Isl1", "GBX2", "GATA5", "GLIS3"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

PlotFootprint(hm.integrated.dfree, features = c("GATA5"))

# Correlation plots
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
hm.integrated.dfree

# ATAC
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")
Idents(hm.integrated.dfree) <- "celltype_sex_ancestry"
av.exp <- AverageExpression(hm.integrated.dfree)$RNA
DefaultAssay(hm.integrated.dfree) <- "RNA"
hm.integrated.dfree <- FindVariableFeatures(hm.integrated.dfree, selection.method = "vst", nfeatures = 2000)
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["RNA"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")

av.exp <- AverageExpression(hm.integrated.dfree)$ATAC
DefaultAssay(hm.integrated.dfree) <- "ATAC"
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["ATAC"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")

av.exp <- AverageExpression(hm.integrated.dfree)$predicted
DefaultAssay(hm.integrated.dfree) <- "predicted"
hm.integrated.dfree <- FindVariableFeatures(hm.integrated.dfree, selection.method = "vst", nfeatures = 2000)
av.exp <- av.exp[rownames(av.exp) %in% hm.integrated.dfree@assays[["predicted"]]@var.features, ]
cor.exp <- as.data.frame(round(cor(av.exp, method = "pearson"),2))
cor.exp$x <- rownames(cor.exp)
datanames_toplot <- unique(hm.integrated.dfree@meta.data[["celltype_sex_ancestry"]])
datanames_toplot <- as.character(datanames_toplot[!is.na(datanames_toplot)])
cor.df <- tidyr::gather(data = cor.exp, y, correlation, datanames_toplot)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(#color = "black"
  ) + theme_bw() + theme(axis.text.x = element_text(angle=70,vjust = 1, hjust=1)) + scale_fill_gradient2(low = "dodgerblue4", high = "red4", mid = "white", 
                                                                                                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                                                                                                         name="Pearson\nCorrelation")


# Load in motifs From 
enriched.motifs.beta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.beta.csv)")
enriched.motifs.delta <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.alpha.csv)")
enriched.motifs.alpha <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.delta.csv)")
enriched.motifs.gamma <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.gamma.csv)")
enriched.motifs.ductal <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.ductal.csv)")
enriched.motifs.acinar <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.acinar.csv)")
enriched.motifs.astel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.astel.csv)")
enriched.motifs.qstel <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.qstel.csv)")
enriched.motifs.endo <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.endo.csv)")
enriched.motifs.lympho <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.lympho.csv)")
enriched.motifs.macro <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\chromvar\enriched.motifs.macro.csv)")

motifs.to.plot <- rownames(hm.integrated.dfree@assays[["chromvar_macs2"]])
combined_processed_atac <- FindVariableFeatures(combined_processed_atac, selection.method = "vst", nfeatures = 800) #changed to 800 from 500 for Figure 4

DefaultAssay(combined_processed_atac) <- "chromvar_macs2"
combined_processed_atac$celltype_sex <- paste(combined_processed_atac$'celltype', combined_processed_atac$'sex', sep = "_")
Idents(combined_processed_atac) <- "celltype_sex" #change to celltype or celltype_sex
pdf(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\R_imageouts\temp_files\figure.pdf)",
    width = 10,
    height = 9)
genes.to.plot <- combined_processed_atac@assays[["chromvar_macs2"]]@var.features #complete geneset
label_genes <- c("MA0132.2", #PDX1
                 "MA0148.4", #FOXA1
                 "MA0874.1", #Arx
                 "MA1608.1", #ISL1
                 "MA0674.1", #NKX6-1
                 "MA1645.1", #NKX2-2
                 "MA0117.2", #MAFB
                 "MA0077.1", #SOX9
                 "MA0068.2", #PAX4
                 "MA1618.1", #PTF1A
                 "MA0114.4", #HNF4A
                 "MA0046.2", #HNF1A
                 "MA0084.1", #SRY
                 "MA0007.3", #Ar
                 "MA0098.3", #ETS1
                 "MA1484.1"  #ETS2
                 )
heatmap <- dittoHeatmap(
  object = combined_processed_atac,#(subset(combined_processed_rna, idents = c("alpha"))),
  genes = genes.to.plot, #this is a compete gene set
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("celltype", "sex", "ancestry"),
  order.by = c("celltype", "sex"),
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
                                        "activated" = "orange",
                                        "alpha"= "lightseagreen",
                                        "beta" = "dodgerblue3",
                                        "delta" = "chartreuse3",
                                        "ductal" = "darkorange2",
                                        "endothelial" = "red",
                                        "gamma" = "springgreen4",
                                        "lymphocyte" = "orchid1",
                                        "macrophage" = "magenta3",
                                        "quiescent" = "salmon"),
                           sex = c("female" = "red4",
                                   "male" = "deepskyblue3"),
                           ancestry = c("white" = "deepskyblue3",
                                        "black" = "black")),
  # # data.out = FALSE,
 #highlight.features = c("MA0132.2", #PDX1
  #                      "MA0148.4",  #FOXA1
   #                     "MA0874.1",  #Arx
    #                     "MA0674.1" #NKX6-1
                         #gamma
                         #ductal
                         #acinar
                         #activated
                         #quiescent
                         #endo
                         #lymphocyte
   #                      ), #macrophages
  #right_annotations = rowAnnotation(foo = anno_mark(at = c(1), labels = c("HHEX"))),
  # show_colnames = isBulk(object),
  show_rownames = FALSE,
  # scale = "row",
  cluster_row = TRUE,
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
  #gaps_col = c(15, 30, 45, 60, 75, 90, 105, 120, 135, 150),
  #gaps_row = c(1819, 3205, 5170, 5405, 10850, 14650, 21250, 23500, 25750, 26200)
  #gaps_row = c(2804, 4584, 7464, 7799, 17833, 23493, 36665, 39029, 41669, 42089)
) + rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(combined_processed_atac[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=0.3),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1))

heatmap

dev.off()
dev.off()


# Plot motifs
# look at the activity of Mef2c
DefaultAssay(hm.integrated.dfree) <- "chromvar_macs2"
FeaturePlot(
  object = hm.integrated.dfree,
  features = "MA0077.1",
  min.cutoff = 0,
  max.cutoff = 1,
  cols = c("lightgrey", "red4"),
  #order = TRUE,
  raster = TRUE,
  pt.size = 1,
  raster.dpi = c(1024, 1024)
)

FeaturePlot(
  object = processed_rna,
  features = "SOX9",
  min.cutoff = 0,
  max.cutoff = 1,
  cols = c("lightgrey", "red4"),
  #order = TRUE,
  raster = TRUE,
  pt.size = 1,
  raster.dpi = c(1024, 1024)
)

MotifPlot(
  object = hm.integrated.dfree,
  motifs = "MA0132.2",
  assay = 'macs2'
)

# Dotplots
# Selected genes
DefaultAssay(hm.integrated.dfree) <- "chromvar"
markers.to.plot <- c("MA0132.2", #PDX1
                     "MA0148.4",  #FOXA1
                     "MA0874.1",  #Arx
                     "MA1608.1", #ISL1
                     "MA0674.1", #NKX6-1
                     "MA1645.1", #NKX2-2
                     "MA0117.2", #MAFB
                     "MA0077.1", #SOX9
                     "MA0068.2", #PAX4
                     "MA1618.1", #PTF1A
                     "MA0114.4", #HNF4A
                     "MA0046.2", #HNF1A
                     "MA0084.1", #SRY
                     "MA0007.3", #Ar
                     "MA0098.3", #ETS1
                     "MA1484.1" #ETS2
                     )

# Dotplot
DotPlot(hm.integrated.dfree,  
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

# Finding enriched motifs
enriched.motifs.xist <- FindMotifs(
  object = hm.integrated.dfree,
  features = c("chrX-73849537-73853078")
)

write.csv(enriched.motifs.xist, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\ORA\Motifs\xist.csv)")

# Plot
MotifPlot(
  object = hm.integrated.dfree,
  motifs = c("GLIS2")
)

# Alternative to calculating 
bone <- Footprint(
  object = bone,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg19
)


# Translating Motifs
motif_id <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\allsex\motif.ID.csv)", row.names = 1)
table1.df <- merge(motif_id, enriched.cvar.motifs.beta,
                   by = 'row.names', all = TRUE)
enriched.cvar.motifs.beta

# Swict frag paths
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.qs)")

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

# Save file
qsave(hm.integrated.dfree, r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.qs)")

# Stop

# All genes
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.qs)")

# DE on this PC because running cicero on the other
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
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "beta", 
    ident.2 = c("alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-308
  write.csv(beta_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\beta_peaks.csv)")
  
  
  alpha_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "alpha", 
    ident.2 = c("beta", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-308
  write.csv(alpha_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\alpha_peaks.csv)")
  
  delta_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "delta", 
    ident.2 = c("beta", "alpha", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-308
  write.csv(delta_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\delta_peaks.csv)")
  
  gamma_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "gamma", 
    ident.2 = c("beta", "alpha", "delta", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-308
  write.csv(gamma_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\gamma_peaks.csv)")
  
  acinar_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "acinar", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-308
  write.csv(acinar_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\acinar_peaks.csv)")
  
  ductal_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "ductal", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-308
  write.csv(ductal_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\ductal_peaks.csv)")
  
  activatedstellate_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "activated_stellate", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "quiescent_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  activatedstellate_peaks$p_val_adj[activatedstellate_peaks$p_val_adj == 0] <- 2e-308
  write.csv(activatedstellate_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\activatedstellate_peaks.csv)")
  
  quiescentstellate_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "quiescent_stellate", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "endothelial", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  quiescentstellate_peaks$p_val_adj[quiescentstellate_peaks$p_val_adj == 0] <- 2e-308
  write.csv(quiescentstellate_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\quiescentstellate_peaks.csv)")
  
  endothelial_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "endothelial", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "macrophage", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  endothelial_peaks$p_val_adj[endothelial_peaks$p_val_adj == 0] <- 2e-308
  write.csv(endothelial_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\endothelial_peaks.csv)")
  
  macrophage_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "macrophage", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "lymphocyte"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  macrophage_peaks$p_val_adj[macrophage_peaks$p_val_adj == 0] <- 2e-308
  write.csv(macrophage_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\macrophage_peaks.csv)")
  
  lymphocyte_peaks <- FindMarkers(
    object = hm.integrated.dfree,
    logfc.threshold = 0.263034406, #1.2
    densify = TRUE,
    only.pos = TRUE,
    ident.1 = "lymphocyte", 
    ident.2 = c("beta", "alpha", "delta", "gamma", "acinar", "ductal", "activated_stellate", "quiescent_stellate", "endothelial", "macrophage"),
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_RNA_macs2'
  )
  lymphocyte_peaks$p_val_adj[lymphocyte_peaks$p_val_adj == 0] <- 2e-308
  write.csv(lymphocyte_peaks, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\allsex\lymphocyte_peaks.csv)")
}




beta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\beta_mvf_peaks.csv)", sep = ",", row.names = 1)
alpha_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\alpha_mvf_peaks.csv)", sep = ",", row.names = 1)
delta_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\delta_mvf_peaks.csv)", sep = ",", row.names = 1)
gamma_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\gamma_mvf_peaks.csv)", sep = ",", row.names = 1)
acinar_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\acinar_mvf_peaks.csv)", sep = ",", row.names = 1)
ductal_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\ductal_mvf_peaks.csv)", sep = ",", row.names = 1)
quiescentstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\quiescentstellate_mvf_peaks.csv)", sep = ",", row.names = 1)
activatedstellate_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\activatedstellate_mvf_peaks.csv)", sep = ",", row.names = 1)
macrophage_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\macrophage_mvf_peaks.csv)", sep = ",", row.names = 1)
lymphocyte_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\lymphocyte_mvf_peaks.csv)", sep = ",", row.names = 1)
endothelial_peaks <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_peaks_macs2\sex\endothelial_mvf_peaks.csv)", sep = ",", row.names = 1)
{
  beta_peaks$p_val_adj[beta_peaks$p_val_adj == 0] <- 2e-302
  beta_peaks <- dplyr::filter(beta_peaks, p_val_adj < 1e-5) 
  open_beta_peaks <- rownames(beta_peaks[beta_peaks$avg_log2FC > 1, ])
  
  alpha_peaks$p_val_adj[alpha_peaks$p_val_adj == 0] <- 2e-302
  alpha_peaks <- dplyr::filter(alpha_peaks, p_val_adj < 1e-5) 
  open_alpha_peaks <- rownames(alpha_peaks[alpha_peaks$avg_log2FC > 1, ])
  
  delta_peaks$p_val_adj[delta_peaks$p_val_adj == 0] <- 2e-302
  delta_peaks <- dplyr::filter(delta_peaks, p_val_adj < 1e-5) 
  open_delta_peaks <- rownames(delta_peaks[delta_peaks$avg_log2FC > 1, ])
  
  gamma_peaks$p_val_adj[gamma_peaks$p_val_adj == 0] <- 2e-302
  gamma_peaks <- dplyr::filter(gamma_peaks, p_val_adj < 1e-5) 
  open_gamma_peaks <- rownames(gamma_peaks[gamma_peaks$avg_log2FC > 1, ])
  
  acinar_peaks$p_val_adj[acinar_peaks$p_val_adj == 0] <- 2e-302
  acinar_peaks <- dplyr::filter(acinar_peaks, p_val_adj < 1e-5) 
  open_acinar_peaks <- rownames(acinar_peaks[acinar_peaks$avg_log2FC > 1, ])
  
  ductal_peaks$p_val_adj[ductal_peaks$p_val_adj == 0] <- 2e-302
  ductal_peaks <- dplyr::filter(ductal_peaks, p_val_adj < 1e-5) 
  open_ductal_peaks <- rownames(ductal_peaks[ductal_peaks$avg_log2FC > 1, ])
  
  quiescentstellate_peaks$p_val_adj[quiescentstellate_peaks$p_val_adj == 0] <- 2e-302
  quiescentstellate_peaks <- dplyr::filter(quiescentstellate_peaks, p_val_adj < 1e-5) 
  open_quiescentstellate_peaks <- rownames(quiescentstellate_peaks[quiescentstellate_peaks$avg_log2FC > 1, ])
  
  activatedstellate_peaks$p_val_adj[activatedstellate_peaks$p_val_adj == 0] <- 2e-302
  activatedstellate_peaks <- dplyr::filter(activatedstellate_peaks, p_val_adj < 1e-5) 
  open_activatedstellate_peaks <- rownames(activatedstellate_peaks[activatedstellate_peaks$avg_log2FC > 1, ])
  
  macrophage_peaks$p_val_adj[macrophage_peaks$p_val_adj == 0] <- 2e-302
  macrophage_peaks <- dplyr::filter(macrophage_peaks, p_val_adj < 1e-5) 
  open_macrophage_peaks <- rownames(macrophage_peaks[macrophage_peaks$avg_log2FC > 1, ])
  
  lymphocyte_peaks$p_val_adj[lymphocyte_peaks$p_val_adj == 0] <- 2e-302
  lymphocyte_peaks <- dplyr::filter(lymphocyte_peaks, p_val_adj < 1e-5) 
  open_lymphocyte_peaks <- rownames(lymphocyte_peaks[lymphocyte_peaks$avg_log2FC > 1, ])
  
  endothelial_peaks$p_val_adj[endothelial_peaks$p_val_adj == 0] <- 2e-302
  endothelial_peaks <- dplyr::filter(endothelial_peaks, p_val_adj < 1e-5) 
  open_endothelial_peaks <- rownames(endothelial_peaks[endothelial_peaks$avg_log2FC > 1, ])
  
  cgenes_beta <- ClosestFeature(hm.integrated.dfree, regions = open_beta_peaks)
  cgenes_alpha <- ClosestFeature(hm.integrated.dfree, regions = open_alpha_peaks)
  cgenes_delta <- ClosestFeature(hm.integrated.dfree, regions = open_delta_peaks)
  cgenes_gamma <- ClosestFeature(hm.integrated.dfree, regions = open_gamma_peaks)
  cgenes_acinar <- ClosestFeature(hm.integrated.dfree, regions = open_acinar_peaks)
  cgenes_ductal <- ClosestFeature(hm.integrated.dfree, regions = open_ductal_peaks)
  cgenes_qstel <- ClosestFeature(hm.integrated.dfree, regions = open_quiescentstellate_peaks)
  cgenes_astel <- ClosestFeature(hm.integrated.dfree, regions = open_activatedstellate_peaks)
  cgenes_macro <- ClosestFeature(hm.integrated.dfree, regions = open_macrophage_peaks)
  cgenes_lympho <- ClosestFeature(hm.integrated.dfree, regions = open_lymphocyte_peaks)
  cgenes_endo <- ClosestFeature(hm.integrated.dfree, regions = open_endothelial_peaks)
  
  cgenes_beta <- dplyr::filter(cgenes_beta, distance < 100000) 
  cgenes_alpha <- dplyr::filter(cgenes_alpha, distance < 100000) 
  cgenes_delta <- dplyr::filter(cgenes_delta, distance < 100000) 
  cgenes_gamma <- dplyr::filter(cgenes_gamma, distance < 100000) 
  cgenes_acinar <- dplyr::filter(cgenes_acinar, distance < 100000) 
  cgenes_ductal <- dplyr::filter(cgenes_ductal, distance < 100000) 
  cgenes_qstel <- dplyr::filter(cgenes_qstel, distance < 100000) 
  cgenes_astel <- dplyr::filter(cgenes_astel, distance < 100000) 
  cgenes_macro <- dplyr::filter(cgenes_macro, distance < 100000) 
  cgenes_lympho <- dplyr::filter(cgenes_lympho, distance < 100000) 
  cgenes_endo <- dplyr::filter(cgenes_endo, distance < 100000)
}


head(beta_peaks)

# FirP
Identshm.integrated.dfree
Idents
 plotingvln <- VlnPlot(
  hm.integrated.dfree,
  features = "FRiP",
  cols = NULL,
  pt.size = 0.1,
  #alpha = 1,
  #idents = ,
  sort = FALSE,
  assay = NULL,
  #group.by = "sex",
  split.by = "celltype_sex_ancestry")

############################ STAGE ############################
############################   XX  ############################
# Pando GRN
# Load object
hm.integrated.dfree <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\hm.integrated.dfree.macs2.chromvarmacs2.conns.qs)")
Idents(hm.integrated.dfree) <- "sex"
male.atac <- subset(hm.integrated.dfree, idents = c("male"))
female.atac <- subset(hm.integrated.dfree, idents = c("female"))
 
# Run seperatly
# Initiate Pando GRN
# hm.integrated.dfree <- initiate_grn(hm.integrated.dfree, 
#                                     regions = union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38), 
#                                     peak_assay = "macs2", 
#                                     rna_assay = "RNA_macs2", 
#                                     #exclude_exons = TRUE
#                                     )

male.atac <- initiate_grn(male.atac, 
                          regions = union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38), 
                          peak_assay = "macs2", 
                          rna_assay = "RNA_macs2", 
                          #exclude_exons = TRUE
)

female.atac <- initiate_grn(female.atac, 
                            regions = union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38), 
                            peak_assay = "macs2", 
                            rna_assay = "RNA_macs2", 
                            #exclude_exons = TRUE
)

# Find TF binding sites
# hm.integrated.dfree <- find_motifs(
#   hm.integrated.dfree,
#   pfm = motifs,
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )

male.atac <- find_motifs(
  male.atac,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

female.atac <- find_motifs(
  female.atac,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

gc()
gc()

# Find variable genes in RNA_macs2
#DefaultAssay(hm.integrated.dfree) <- "RNA_macs2"
DefaultAssay(male.atac) <- "RNA_macs2"
DefaultAssay(female.atac) <- "RNA_macs2"


#hm.integrated.dfree <- FindVariableFeatures(hm.integrated.dfree, selection.method = "vst", nfeatures = 2000)
male.atac <- FindVariableFeatures(male.atac, selection.method = "vst", nfeatures = 2000)
female.atac <- FindVariableFeatures(female.atac, selection.method = "vst", nfeatures = 2000)


#DefaultAssay(hm.integrated.dfree) <- "macs2"
DefaultAssay(male.atac) <- "macs2"
DefaultAssay(female.atac) <- "macs2"

jaspartfdb <- read.csv(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\TFs\motif.ID.jaspar2020.csv)")
intersecttfs <- intersect(hm.integrated.dfree@assays[["RNA_macs2"]]@counts@Dimnames[[1]], jaspartfdb$X1.TranscripFact)
#qsave(intersecttfs, file = r"(E:\2.SexbasedStudyCurrent\QS files\intersecttfs.qs)")
intersecttfs <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\intersecttfs.qs)")

# Infer GRN
# hm.integrated.dfree <- infer_grn(
#   hm.integrated.dfree,
#   peak_to_gene_method = 'Signac',
#   #genes = hm.integrated.dfree@assays[["RNA_macs2"]]@var.features, #all var genes by default
#   genes = intersecttfs, #TFs for which motifs exist
#   method = 'glm'
# )

male.atac <- infer_grn(
  male.atac,
  peak_to_gene_method = 'Signac',
  #genes = male.atac@assays[["RNA_macs2"]]@var.features, #all var genes by default
  genes = intersecttfs, #TFs for which motifs exist
  method = 'glm'
)

female.atac <- infer_grn(
  female.atac,
  peak_to_gene_method = 'Signac',
  #genes = female.atac@assays[["RNA_macs2"]]@var.features, #all var genes by default
  genes = intersecttfs, #TFs for which motifs exist
  method = 'glm'
)

# Inspect coeff
coef(male.atac)
coef(female.atac)

# Mod extr
male.atac <- find_modules(male.atac)
modules_male <- NetworkModules(male.atac)
modules_male@meta

female.atac <- find_modules(female.atac)
modules_female <- NetworkModules(female.atac)
modules_female@meta

#qsave(male.atac, file = r"(E:\2.SexbasedStudyCurrent\QS files\male.atac.dfree.macs2.chromvarmacs2.conns.grn.qs)")
#qsave(female.atac, file = r"(E:\2.SexbasedStudyCurrent\QS files\female.atac.dfree.macs2.chromvarmacs2.conns.grn.qs)")
male.atac <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\male.atac.dfree.macs2.chromvarmacs2.conns.grn.qs)")
female.atac <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\female.atac.dfree.macs2.chromvarmacs2.conns.grn.qs)")
processed_rna <- qread(r"(E:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

# Visuals
plot_gof(male.atac, point_size=3)
plot_module_metrics(male.atac)

plot_gof(female.atac, point_size=3)
plot_module_metrics(female.atac)

# Graph
DefaultAssay(male.atac) <- "RNA_macs2"
DefaultAssay(female.atac) <- "RNA_macs2"
vargenes_male <- head(male.atac@assays[["RNA_macs2"]]@var.features, 2000)
vargenes_female <- head(female.atac@assays[["RNA_macs2"]]@var.features, 2000)
intersecttfs_vargenes_male <- intersect(intersecttfs, head(male.atac@assays[["RNA_macs2"]]@var.features, 2000))
intersecttfs_vargenes_female <- intersect(intersecttfs, head(female.atac@assays[["RNA_macs2"]]@var.features, 2000))

# GRNs
# Networks
male.atac <- get_network_graph(male.atac, 
                               graph_name='umap_graph',
                               rna_assay = "RNA_macs2",
                               features = intersecttfs
                               )
female.atac <- get_network_graph(female.atac, 
                                 graph_name='umap_graph',
                                 rna_assay = "RNA_macs2",
                                 features = intersecttfs
                                 )

# Graphs
subgraph_male <- NetworkGraph(male.atac, 
                         graph='umap_graph')
subgraph_female <- NetworkGraph(female.atac, 
                              graph='umap_graph')

beta_mvf_genes <- read.table(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\DETesting\DE_testing\alldata\beta.deseq.WaldTest.M_ND.vs.F_ND.tsv)")
beta_mvf_genes <- beta_mvf_genes[rownames(beta_mvf_genes) %in% as.character(intersect(rownames(beta_mvf_genes), as.character(hm.integrated.dfree@grn@networks[["glm_network"]]@features))), ]
beta_mvf_genes <- beta_mvf_genes[rownames(beta_mvf_genes) %in% ggplot_pando[["data"]][["name"]], ] #ggplot_pando is a ggplot object making a global pando GRN

# Plotting
library(ggraph)
library(RColorBrewer)
display.brewer.all()
p1 <- ggraph(subgraph_male, x=UMAP_1, y=UMAP_2) +
  geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) +
  geom_node_point(aes(size=centrality, fill=centrality), shape=21, stroke = 0.2) +
  scale_size(range = c(1,10)) +
  geom_node_text(aes(label=name), size=2, shape=21, repel=F) +
  #scale_fill_gradientn(colors=c('-1'='darkgrey', '1'='orange')) +
  scale_fill_viridis(aes(fill= centrality), option='magma', direction = -1) +
  coord_flip() +
  theme_void()

p2 <- ggraph(subgraph_female, x=UMAP_1, y=UMAP_2) +
  geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) +
  geom_node_point(aes(size=centrality, fill=centrality), shape=21, stroke = 0.2) +
  scale_size(range = c(1,10)) +
  geom_node_text(aes(label=name), size=2, shape=21, repel=F) +
  #scale_fill_gradientn(colors=c('-1'='darkgrey', '1'='orange')) +
  scale_fill_viridis(aes(fill= centrality), option='magma', direction = -1) +
  theme_void()

p1 + p2

plot_network_graph(hm.integrated.dfree,
                   node_size = c(1, 10),
                   graph='umap_graph')

# Plot TFs
# Make average seurat object
Idents(processed_rna) <- "disease_ancestry_lib_sex_source_celltype"
combined_processed_rna <- AverageExpression(processed_rna, return.seurat = TRUE, slot = 'data')
unique(Idents(combined_processed_rna))
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
processed_rna <- NULL

# Subset datasets
combined_processed_rna$celltype_sex_disease <- paste(combined_processed_rna$celltype, combined_processed_rna$sex, combined_processed_rna$disease, sep = '_')
Idents(combined_processed_rna) <- "celltype_sex_disease"
beta_alpha_combined <- subset(combined_processed_rna, idents = c("beta_M_ND", "beta_F_ND", "alpha_M_ND", "alpha_F_ND"))

Idents(beta_alpha_combined) <- "celltype"
beta_cell_genes <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\beta.csv)", row.names = 1, sep = ",")
beta_cell_genes <- beta_cell_genes[rownames(beta_cell_genes) %in% as.character(intersect(rownames(beta_cell_genes), as.character(hm.integrated.dfree@grn@networks[["glm_network"]]@features))), ]
beta_cell_genes <- beta_cell_genes[order(beta_cell_genes$avg_log2FC,decreasing = TRUE), ]
beta_cell_genes <- top_n(beta_cell_genes, 25, avg_log2FC)
beta_cell_genes <- rownames(beta_cell_genes)

alpha_cell_genes <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\scRNA\Conserved markers\DEtesting\DEseq2\alpha.csv)", row.names = 1, sep = ",")
alpha_cell_genes <- alpha_cell_genes[rownames(alpha_cell_genes) %in% as.character(intersect(rownames(alpha_cell_genes), as.character(hm.integrated.dfree@grn@networks[["glm_network"]]@features))), ]
alpha_cell_genes <- alpha_cell_genes[order(alpha_cell_genes$avg_log2FC,decreasing = TRUE), ]
alpha_cell_genes <- top_n(alpha_cell_genes, 25, avg_log2FC)
alpha_cell_genes <- rownames(alpha_cell_genes)
cell_genes <- c(beta_cell_genes, alpha_cell_genes)

markers.to.plot <- cell_genes
Idents(beta_alpha_combined) <- "celltype_sex_disease"
DotPlot(beta_alpha_combined,  
        group.by = "celltype_sex_disease",
        #split.by = "sex",
        dot.scale = 5,
        #scale.min = 0,
        #col.min = 0, #minimum level
        #col.max = 1,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =10, colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =10, colour = "black")) +
  theme(plot.title = element_text(size = 1),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10)) +
  scale_colour_gradient2(low =c("dodgerblue"), mid = c("white"), high =c("red4")) +
  guides(color = guide_colorbar(title = 'Average Expression')) + coord_flip()

# GRN Beta
Idents(hm.integrated.dfree) <- "celltype"
beta_macs <- subset(x = hm.integrated.dfree, idents = "beta")
beta_alpha_macs <- subset(x = hm.integrated.dfree, idents = c("beta", "alpha"))

Idents(beta_macs) <- "sex"
beta_macs_m <- subset(x = beta_macs, idents = "male")
beta_macs_f <- subset(x = beta_macs, idents = "female")

# Read in motifs for all beta cells
male_beta_motifs <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.beta.csv)", row.names = 1, sep = ",")
male_beta_motifs$foldenrich <- male_beta_motifs$pct.1 / male_beta_motifs$pct.2

row.names(male_beta_motifs$foldenrich < 1 & male_beta_motifs$p_val_adj < 0.05)
rownames(male_beta_motifs)[male_beta_motifs$foldenrich < 1 & male_beta_motifs$p_val_adj < 0.05]

male_alpha_motifs <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\Data Output\snATAC\DE_accessible_sites\DA_motifs\sex\chromvar\enriched.motifs.alpha.csv)", row.names = 1, sep = ",")
male_alpha_motifs$foldenrich <- male_alpha_motifs$pct.1 / male_alpha_motifs$pct.2

row.names(male_alpha_motifs$foldenrich < 1 & male_alpha_motifs$p_val_adj < 0.05)
rownames(male_alpha_motifs)[male_alpha_motifs$foldenrich < 1 & male_alpha_motifs$p_val_adj < 0.05]


DefaultAssay(beta_macs) <- "predicted"
VlnPlot(beta_macs, features = c("SRY", "AR", "XIST"), pt.size = 0)
VlnPlot(beta_macs, features = c("MA0084.1", "MA0007.3"), pt.size = 0.1)

beta_alpha_macs <- get_network_graph(beta_alpha_macs, 
                               graph_name='umap_graph',
                               rna_assay = "RNA_macs2",
                               # features = c("DDX3X", "KDM6A", "MAPT", "BARHL1", "MAPT-IT1", "SH3BP5", "AL672277.1", 
                               #              "XIST", "TSIX", "JPX",
                               #              
                               #              "RPS4Y1", "EIF1AY", "USP9Y", "DDX3Y", "TTTY14", "KDM5D", "UTY", "ZFY",
                               #              "AC244213.1", "NLGN4Y", "TMSB4Y", "AC006157.1", "SRY", "LINC00278", "TTTY10", "AC010889.2", "AC010889.1", "AC011297.1", 
                               #              "HAR1B", "CD99", "BCHE", "ODF3L1", "EFEMP2", "DENND1B", "SOCS6", "TWSG1", "ADGRG7"),
                               #features = c(rownames(male_beta_motifs)[male_beta_motifs$p_val_adj < 0.05], rownames(male_alpha_motifs)[male_alpha_motifs$p_val_adj < 0.05])
                               features = VariableFeatures(beta_alpha_macs)
)

plot_network_graph(beta_alpha_macs, graph='umap_graph')

DefaultAssay(beta_alpha_macs) <- "RNA_macs2"
beta_alpha_macs <- FindVariableFeatures(beta_alpha_macs, selection.method = "vst", nfeatures = 500)
top1000 <- head(VariableFeatures(hm.integrated.dfree), 1000)
vargenes <- head(processed_rna@assays[["RNA"]]@var.features, 1000)

#####

## Discover motifs in a gene across sex
Idents(hm.integrated.dfree) <- "sex"
male_atac <- subset(hm.integrated.dfree, idents = c("male"))

# Get all peaks linked to INS gene
peaks_linked_ins <- GetLinkedPeaks(hm.integrated.dfree, features = c("INS-IGF2", "INS", "IGF2", "TH"), assay = NULL, min.abs.score = 0)

CoveragePlot(hm.integrated.dfree, region = c("INS"), 
             extend.upstream = 10000,
             extend.downstream = 10000, links = TRUE)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(male_atac),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)

# example peak set we're interested in
peaks.use <- peaks_linked_ins
motif.use <- colnames(motif.matrix)[1]

# compute fraction of peaks containing a certain motif
sum(motif.matrix[peaks.use, motif.use]) / length(peaks.use)

#





















   
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














link_plot + scale_colour_viridis_d()

+ scale_fill_gradient2(low = "red",
                                 mid = "white",
                                 high = "blue",
                                 na.value = "grey50")





p1[["data"]][["Var2"]])
my_levels2 <- c("beta", "delta", "alpha", "gamma", 
                "acinar", "ductal",
                "quiescent_stellate", "activated_stellate", "endothelial",
                "lymphocyte", "macrophage")
processed_rna$disease_ancestry_lib_sex_source <- factor(x = processed_rna$disease_ancestry_lib_sex_source, levels = my_levels)

p1

predictions <- round(cor(predictions),2)
head(predictions)
melted_predictions <- melt(predictions)
head(melted_predictions)
ggplot(data = melted_predictions, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()



predictions <- t(predictions)
predictions <- as.data.frame(predictions)




predictions <- as.data.frame(cor(predictions, method = "spearman"))


predictions$x <- rownames(predictions)
predictions <- tidyr::gather(data = predictions, y, correlation, c('beta', 'alpha', 'delta', 'gamma', 'acinar', 'ductal', 'quiescent_stellate', 'activated_stellate', 'endothelial', 'macrophage'))
ggplot(predictions, aes(x, y, fill = correlation)) +
  geom_tile()

ggplot(predictions, aes(x= Var1, y = Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2

# END #
> sessionInfo()
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3          rstudioapi_0.14             jsonlite_1.8.4              magrittr_2.0.3             
  [5] spatstat.utils_3.0-3        BiocIO_1.10.0               fs_1.6.2                    zlibbioc_1.46.0            
  [9] vctrs_0.6.2                 ROCR_1.0-11                 Rsamtools_2.16.0            memoise_2.0.1              
 [13] spatstat.explore_3.2-1      RCurl_1.98-1.12             S4Arrays_1.0.4              htmltools_0.5.5            
 [17] usethis_2.1.6               sctransform_0.3.5.9002      parallelly_1.35.0           KernSmooth_2.23-21         
 [21] htmlwidgets_1.6.2           ica_1.0-3                   plyr_1.8.8                  plotly_4.10.1              
 [25] zoo_1.8-12                  cachem_1.0.8                GenomicAlignments_1.36.0    igraph_1.4.2               
 [29] mime_0.12                   lifecycle_1.0.3             pkgconfig_2.0.3             Matrix_1.5-4               
 [33] R6_2.5.1                    fastmap_1.1.1               MatrixGenerics_1.12.0       GenomeInfoDbData_1.2.10    
 [37] fitdistrplus_1.1-11         future_1.32.0               shiny_1.7.4                 digest_0.6.31              
 [41] colorspace_2.1-0            S4Vectors_0.38.1            patchwork_1.1.2             ps_1.7.5                   
 [45] Seurat_4.3.0.9002           tensor_1.5                  irlba_2.3.5.1               pkgload_1.3.2              
 [49] GenomicRanges_1.52.0        progressr_0.13.0            fansi_1.0.4                 spatstat.sparse_3.0-1      
 [53] httr_1.4.6                  polyclip_1.10-4             abind_1.4-5                 compiler_4.3.0             
 [57] remotes_2.4.2               BiocParallel_1.34.1         pkgbuild_1.4.0              MASS_7.3-60                
 [61] DelayedArray_0.26.2         sessioninfo_1.2.2           rjson_0.2.21                tools_4.3.0                
 [65] lmtest_0.9-40               httpuv_1.6.11               future.apply_1.11.0         goftest_1.2-3              
 [69] glue_1.6.2                  restfulr_0.0.15             callr_3.7.3                 nlme_3.1-162               
 [73] promises_1.2.0.1            grid_4.3.0                  Rtsne_0.16                  cluster_2.1.4              
 [77] reshape2_1.4.4              generics_0.1.3              gtable_0.3.3                spatstat.data_3.0-1        
 [81] tidyr_1.3.0                 data.table_1.14.8           XVector_0.40.0              sp_1.6-0                   
 [85] utf8_1.2.3                  BiocGenerics_0.46.0         spatstat.geom_3.2-1         RcppAnnoy_0.0.20           
 [89] ggrepel_0.9.3               RANN_2.6.1                  pillar_1.9.0                stringr_1.5.0              
 [93] later_1.3.1                 splines_4.3.0               dplyr_1.1.2                 lattice_0.21-8             
 [97] survival_3.5-5              rtracklayer_1.60.0          deldir_1.0-9                tidyselect_1.2.0           
[101] Biostrings_2.68.1           miniUI_0.1.1.1              pbapply_1.7-0               gridExtra_2.3              
[105] IRanges_2.34.0              SummarizedExperiment_1.30.1 scattermore_1.1             stats4_4.3.0               
[109] Biobase_2.60.0              devtools_2.4.5              matrixStats_0.63.0          stringi_1.7.12             
[113] yaml_2.3.7                  lazyeval_0.2.2              codetools_0.2-19            tibble_3.2.1               
[117] BiocManager_1.30.20         cli_3.6.1                   uwot_0.1.14                 xtable_1.8-4               
[121] reticulate_1.28             munsell_0.5.0               processx_3.8.1              GenomeInfoDb_1.36.0        
[125] Rcpp_1.0.10                 globals_0.16.2              spatstat.random_3.1-5       png_0.1-8                  
[129] XML_3.99-0.14               parallel_4.3.0              ellipsis_0.3.2              ggplot2_3.4.2              
[133] prettyunits_1.1.1           profvis_0.3.8               urlchecker_1.0.1            bitops_1.0-7               
[137] listenv_0.9.0               viridisLite_0.4.2           scales_1.2.1                ggridges_0.5.4             
[141] SeuratObject_4.1.3          leiden_0.4.3                purrr_1.0.1                 crayon_1.5.2               
[145] rlang_1.1.1                 cowplot_1.1.1              
Warning message:
replacing previous import ‘S4Arrays::read_block’ by ‘DelayedArray::read_block’ when loading ‘SummarizedExperiment’ 

