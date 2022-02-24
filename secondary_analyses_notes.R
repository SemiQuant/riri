---
title: "RiboSeek Anlaysis"
author: "SemiQuant"
date: "2/23/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r packages, include=FALSE, eval=FALSE}
install.packages(c("plotly", "tidyverse", "BiocManager", "DESeq2", "parallel", "xtail", "devtools"))
BiocManager::install(c("riboSeqR", "systemPipeR"))
devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
```

## Inspiration

- https://bioconductor.org/packages/release/bioc/html/riboSeqR.html

- https://github.com/xryanglab/xtail

- https://bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeRIBOseq.html

- https://github.com/LabTranslationalArchitectomics/riboWaltz

- https://github.com/ratschlab/ribodiff

- https://pubmed.ncbi.nlm.nih.gov/19213877/#&gid=article-figures&pid=fig-4-uid-3


## Read in processed data and check/clean