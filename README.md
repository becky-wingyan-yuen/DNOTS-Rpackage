# On Data Normalization for Tumor Subtyping with microRNA data (DNOTS)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)

# Overview
Single-cell RNA-sequencing (scRNA-seq) technologies enable the measurement of the transcriptome of individual cells, which provides unprecedented opportunities to discover cell types and understand cellular heterogeneity. Despite their widespread applications, single-cell RNA-sequencing (scRNA-seq) experiments are still plagued by batch effects and dropout events.

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat and the surrogate variable analysis (SVA), are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, the R package `BUSseq` fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data (BUSseq)---to correct batch effects in the presence of unknown cell types. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate read count matrices obtained from different scRNA-seq platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions listed in our [manuscript](https://www.biorxiv.org/content/10.1101/533372v3).

# Repo Contents

- [R](./R): `R` code.
- [data](./data): the example data for the demo.
- [man](./man): help manual.

### Software dependencies

Before installing the `DNOTS` package, users should have installed `R` with version 3.5.0 or higher.


#### Package dependencies

Users should install the following packages prior to installing `DNOTS`, from an `R` session:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva")
BiocManager::install("vsn")
```

# Installation Guide

From an `R` session, type:

```
require(devtools)
install_github("becky-wing-yan-yuen/DNOTS-Rpackage") 
```


# Demo

Please check the user's guide for the detailed instructions on how to use the functions in the package. You may also run the following code in the `R` session:

```
data("example_data", package = "DNOTS")

uni_true_cluster = rep(NA,dim(example_data)[2])
uni_true_cluster[which(substr(colnames(example_data),7,7)=="E")] = 1
uni_true_cluster[which(substr(colnames(example_data),7,7)=="V")] = 2

results = cluster_other(dat = example_data, true_cluster = uni_true_cluster, 
						clust_method = "kmeans")
						
str(results)
```


