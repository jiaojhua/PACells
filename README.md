# PACells: Phenotype-Associated Cells
PACells identifies clinical phenotype-associated cell states by integrating bulk and single-cell chromatin accessibility profiles

## Overview
We introduce PACells, 

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
To run PACells, it needs to install a list of R packages, including: Seurat, chromVAR, SGL, chromVARmotifs, motifmatchr, Signac, SummarizedExperiment, gam, survival, snowfall.

# Contact
If you have any suggestions or problems, please contact Jiao Hua (jhua@stu.hit.edu.cn).
