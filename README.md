# PACells: Phenotype-Associated Cells
PACells identifies clinical phenotype-associated cell states by integrating bulk and single-cell chromatin accessibility profiles

## Overview
We introduce PACells, a computational framework to identify cell states associated with clinical phenotypes from scATAC-seq data guided by bulk ATAC-seq data annotated with clinical phenotype information. PACells links phenotype information measured in bulk samples to individual cells in single-cell data and outputs the subset of cells most associated with the phenotype of interest (e.g., disease status, prognosis, treatment response), helping users interpret clinically relevant cellular heterogeneity. PACells also provides an optional extension for scRNA-seq guided by bulk RNA-seq data.

<p align="center">
<img  src="tutorial/PACells.png" width="800" height=auto > 
</p>

## Installation

PACells is provided as an R package and can be installed directly from GitHub.

### Install from GitHub

```R
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("jiaojhua/PACells")
```
## Verify installation

```R
library(PACells)
packageVersion("PACells")
```

## Quick Start
### scATAC-seq (bulk ATAC-seq + scATAC-seq)

```R
library(PACells)

# Prepare motif set (JASPAR or cisBP)
motifs <- getMotifs(database = "JASPAR")

# phenotype should be aligned to bulk samples (columns of bulk_dataset)
sc_res <- PACells(
  sc_dataset   = sc_dataset,
  bulk_dataset = bulk_dataset,
  phenotype    = phenotype,
  motifs       = motifs,
  family       = "binomial",  
  method       = "KL",    
  cutoff       = 0.1,
  screenRatio  = 0.8,
  batch        = "none"  
)

table(sc_res$PACells_label)
```

**Tip (batch effect)**: if cells are mainly separated by batchs, consider batch = "harmony".

### scRNA-seq (bulk RNA-seq + scRNA-seq)

```R
library(PACells)

# sc_mat: genes x cells; bulk_mat: genes x samples
sc_res_rna <- PACells.RNA(
  sc_dataset   = sc_mat,
  bulk_dataset = bulk_mat,
  phenotype    = phenotype,
  family       = "binomial",
  method       = "KL",
  cutoff       = 0.1,
  screenRatio  = 0.8,
  batch        = "none"
)

table(sc_res_rna$PACells_label)
```


# Tutorials

* For CLL datasets, please see [here](https://github.com/jiaojhua/PACells/blob/main/tutorial/Tutorial_CLL.ipynb), datasets are available at this [link](https://drive.google.com/drive/folders/1PpDxiRl8wv2JUdCtBd146cWche-o3U5e?usp=sharing).

* For AD datasets, please see [here](https://github.com/jiaojhua/PACells/blob/main/tutorial/Tutorial_AD.ipynb), datasets are available at this [link](https://drive.google.com/drive/folders/1PpDxiRl8wv2JUdCtBd146cWche-o3U5e?usp=sharing).

# Dependencies
- Seurat, Signac
- chromVAR, chromVARmotifs, motifmatchr
- SummarizedExperiment
- SGL, gam, survival
- snowfall

# Contact
If you have any suggestions or problems, please feel free to contact Jiao Hua (jhua@stu.hit.edu.cn).
