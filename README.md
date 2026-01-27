# PACells: Phenotype-Associated Cells
PACells identifies clinical phenotype-associated cell states by integrating bulk and single-cell chromatin accessibility profiles

## Overview
We introduce PACells, a computational framework to identify cell states associated with clinical phenotypes from scATAC-seq data guided by bulk ATAC-seq data annotated with clinical phenotype information. PACells links phenotype information measured in bulk samples to individual cells in single-cell data and outputs the subset of cells most associated with the phenotype of interest (e.g., **disease status**, **prognosis**, **treatment response**), helping users interpret clinically relevant cellular heterogeneity. PACells also provides an optional extension for **scRNA-seq** guided by bulk RNA-seq data.

<p align="center">
<img  src="tutorial/PACells.png" width="800" height=auto > 
</p>

# Installation
To run ``PACells`` R package, install from GitHub through ``devtools`` directly:
```R
install.packages('devtools')
library(devtools)
devtools::install_github("jiaojhua/PACells")
```

# Tutorials

* For CLL datasets, please see [here](https://github.com/jiaojhua/PACells/blob/main/tutorial/Tutorial_CLL.ipynb), datasets are available at this [link](https://drive.google.com/drive/folders/1PpDxiRl8wv2JUdCtBd146cWche-o3U5e?usp=sharing).

* For AD datasets, please see [here](https://github.com/jiaojhua/PACells/blob/main/tutorial/Tutorial_AD.ipynb), datasets are available at this [link](https://drive.google.com/drive/folders/1PpDxiRl8wv2JUdCtBd146cWche-o3U5e?usp=sharing).

# Dependencies
- Seurat
- chromVAR
- SGL
- chromVARmotifs
- motifmatchr
- Signac
- SummarizedExperiment
- gam
- survival
- snowfall

# Contact
If you have any suggestions or problems, please feel free to contact Jiao Hua (jhua@stu.hit.edu.cn).
