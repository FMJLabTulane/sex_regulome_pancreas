
# Sex-Specific Regulatory Architecture of Pancreatic Islets in Type 2 Diabetes

This repository contains the analysis scripts and data associated with the study "Sex-specific regulatory architecture of pancreatic islets from subjects with and without type 2 diabetes" by Qadir et al., published in *The EMBO Journal* in 2024.

## Table of Contents

- [Overview](#overview)
- [Data Access](#data-access)
  - [Downloading Data Files](#downloading-data-files)
  - [Data Structure](#data-structure)
- [Analysis Tools](#analysis-tools)
- [Setting Up the Environment](#setting-up-the-environment)
- [Computing Environment](#computing-environment)
- [Contributors](#contributors)
- [Contact Information](#contact-information)

## Overview

This repository provides the code and data used to analyze sex-specific differences in the regulatory architecture of pancreatic islets from individuals with and without type 2 diabetes (T2D). The study integrates single-cell RNA sequencing (scRNA-seq) and single-nucleus ATAC sequencing (snATAC-seq) to explore chromatin accessibility and gene expression patterns. The goal is to understand how biological sex influences β-cell function and the pathogenesis of T2D.

## Data Access

### Downloading Data Files

The raw and processed data supporting the findings of this study are available in public repositories. Please refer to the publication for specific links and access instructions.

### Data Structure

The data is organized into the following categories:

#### plotting_rubric

- **plotting rubric for scRNAseq figures

#### plotting_rubric_ATAC

- **plotting rubric for snATACseq figures

#### scRNAseq_analysis

- **Analysis rubric for scRNASeq

#### snATAC_analysis_cellranger

- **Analysis rubric for snATACseq

#### snATACanalysis_macs2

- **Analysis rubric for snATACseq MACS2 pipeline
  
## Analysis Tools

The following tools were used for data analysis and visualization:

### Cell Ranger

- **Description:** A set of analysis pipelines from 10x Genomics for processing scRNA-seq and snATAC-seq data.
- **Usage:** Employed for initial data processing, including demultiplexing, alignment, and generation of feature-barcode matrices.

### R and RStudio

- **Description:** R is a programming language for statistical computing; RStudio is an integrated development environment for R.
- **Usage:** Used for downstream data analysis and visualization.

### Seurat

- **Description:** An R package designed for the analysis and visualization of single-cell RNA-seq data.
- **Usage:** Utilized for clustering, differential expression analysis, and integration of scRNA-seq data.

## Setting Up the Environment

To replicate the analyses, ensure that the following software and packages are installed:

- **Cell Ranger:** Follow installation instructions from 10x Genomics.
- **R and RStudio:** Install the latest versions from CRAN and the RStudio website.
- **Seurat:** Install via CRAN or GitHub using `install.packages("Seurat")` or `devtools::install_github("satijalab/seurat")`.

Additional R packages required include:

- **dplyr**
- **ggplot2**
- **Matrix**
- **cowplot**

Install these packages using `install.packages("package_name")`.

## Computing Environment

Some analyses, particularly those involving large datasets, may require high-performance computing resources. Ensure that your system has sufficient memory and processing power to handle computational demand.

## Contributors

- **Mirza Muhammad Fahd Qadir**
- **Ruth M. Elgamal**

## Contact Information

For questions or further information, please contact:

- **Franck Mauvais-Jarvis**
- **Email:** [email protected]

---

This README is adapted from the AR-DHT repository's README file.
