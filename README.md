# Ciona 16

The R code underlying the analysis of a scientific paper is hosted here in the interests of reproducibility. The manuscript is undergoing peer review and this code might be refactored or changed. 

---
Ilsley, G. R., Suyama, R., Noda, T., Satoh, N. & Luscombe, N. M. Identifying developmentally important genes with single-cell RNA-seq from an embryo. bioRxiv 197699 (2017). https://doi.org/10.1101/197699

---

The code was tested on R 3.4.1 using Bioconductor 3.6 (development version). 

The project has the following directory structure:

| Directory     | Content     |
|:------------- |:------------|
| data-raw      | Data in csv format for importing into R.  | 
| data          | Count data in rda format (of SingleCellExperiment class).  | 
| exports       | Data exported from the DESeq analysis script.  |
| images        | Images generated for Figure 1. |
| scripts       | R scripts used for Figure 1, pattern scoring and DESeq comparison. |

