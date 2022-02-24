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
install.packages(c("plotly", "tidyverse", "BiocManager", "DESeq2", "parallel", "xtail", "devtools",
  "ggplot2", "heatmaply", "adegenet"), quietly = T)
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


## Read in processed data and metadata, and check/clean


### Deseq
```{r}
dds_1 <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                # tidy = T,
                                design = ~ Cell + TBstatus)
dds_1 <- estimateSizeFactors(dds_1)

vst_1_blind <- vst(dds_1, blind = T)
vst_1_blind_cor <- cor(assay(vst_1_blind))

# Heatmap of all transcripts color by groups ------------------------------------------------
hmap_QC_1 <- heatmaply(vst_1_blind_cor,
                    symm = T,
                    col_side_colors = meta$phen,
                    row_side_colors = meta$phen,
                    hclust_method = "complete",
                    dist_method = "manhattan",
                    fontsize_row = 5,
                    fontsize_col = 5,
                    # row_text_angle = 45,
                    plot_method = "plotly"
) 



PCA <- dudi.pca(df = mat, center = center,
                      scale = scale)

p1 <- qplot(data = PCA$l1, x = RS1, y = RS2, colour = pheno)
      
      
      if (circ == T) {
        p1 <- p1 + stat_ellipse(geom = "polygon", 
                                level = 0.95,
                                type = "norm",
                                alpha = 0.3, aes(fill = pheno)) +
          stat_ellipse(type = "euclid", linetype = 2, alpha = 0.3) # level = 0.2) +
      }
      
      p1 <- p1 + labs(col = ifelse(is.null(name_in), "Phenotype", name_in)) 
      
      if (!is.factor(pheno) & !is.character(pheno))
        p1 <- p1 + scale_colour_viridis_c()
      # else
      #   p1 <- p1 + scale_colour_viridis_d()
      p1 <- ggplotly(p1 +
                       #              + theme(
                       # panel.background = element_rect(fill = "#272b30"),
                       # plot.background = element_rect(fill = "#272b30"),
                       # legend.background = element_rect(fill = "#272b30"),
                       # panel.grid.major = element_line(colour = "#c8c8c8"),
                       # axis.title.x = element_text(colour = "#c8c8c8"),
                       # axis.title.y = element_text(colour = "#c8c8c8"),
                       # legend.title = element_text(colour = "#c8c8c8"),
                       # plot.title = element_text(colour = "#c8c8c8"),
                       # legend.text = element_text(colour = "#c8c8c8"),
                       # axis.text = element_text(colour = "#c8c8c8")
                     # ) +
                     labs(title = title) + #, subtitle = "A subtitle"
                       xlab(label = paste("PC1", "(", round(pca_sum$variance.percent[1], 0), "%)")) +
                       ylab(label =  paste("PC2", "(", round(pca_sum$variance.percent[2], 0), "%)")
                       )
      ) %>% 
        layout(
          legend = list(
            orientation = 'h',
            y=-0.15)
        )
      p1 


# PCA
PCA_vst_1 <- dudi.pca(df = t(assay(vst_1_blind)), center = T, scale = F, scannf = F, nf = 3)
pca_QC_1 <- pca_qplot(pheno = meta$phen,
                   type = "plotly",
                   title = "PCA",
                   pca_in = PCA_vst_1,
                   name_in = "Group",
                   circ = T,
                   scale = F,
                   # tsne = T,
                   center = T
)



# Corr heatmap ---------------
vst_1_blind_cor <- cor(assay(vst_1_blind_QC_1))

(hp_corr <- heatmaply(vst_1_blind_cor,
          symm = T,
          col_side_colors = meta_QC_1$phen,
          row_side_colors = meta_QC_1$phen,
          hclust_method = "complete",
          dist_method = "euclidean", #"manhattan",
          fontsize_row = 5,
          fontsize_col = 5,
          k_col=9,
          # k_row = 4,
          dendrogram = "both",
          # row_text_angle = 45,
          plot_method = "plotly"
))







run_once=0
if (run_once == 1){
  load("QC_done.Rimg")
  meta_QC_1$Cell <- relevel(meta_QC_1$Cell, "RM") # this is setting RM as the reference
  meta_QC_1$TBstatus <- relevel(meta_QC_1$TBstatus, "neg")
  dds_in <- DESeqDataSetFromMatrix(countData = counts_QC_1,
                                   colData = meta_QC_1,
                                   design = ~ Cell * TBstatus)
  
  dds_in <- DESeq(dds_in, parallel = T)
  save(dds_in, file = "dds_in.Rdat")
}
load("QC_done.Rimg")
load("dds_in.Rdat")












function(deSeq_results = NULL, adj_p = 0.05, num_annotate = 10, effect_sizeThresh = 2, makePlot = T){
  volc.plot <- data.frame(ID = rownames(deSeq_results),
                          P = deSeq_results$pvalue,
                          P.adj = deSeq_results$padj,
                          EffectSize = deSeq_results$log2FoldChange,
                          threshold = as.factor(deSeq_results$padj <= adj_p)
                          # comparison = c("TB vs No TB")
  )

  volc.plot$threshold <- as.factor(volc.plot$P.adj <= adj_p & abs(volc.plot$EffectSize) >= effect_sizeThresh)


  # rewrite this, in a rush
  effect_sizeThresh_2 <- abs(volc.plot[volc.plot$P.adj <= adj_p,]$EffectSize)
  effect_sizeThresh_2 <- effect_sizeThresh_2[order(abs(effect_sizeThresh_2), decreasing = T)][num_annotate]


  volc.plot$lab <- ifelse(volc.plot$P.adj <= adj_p & abs(volc.plot$EffectSize) >= effect_sizeThresh_2, as.character(volc.plot$ID), "")

  volc_DE_1 <- ggplot2::ggplot(data = volc.plot,
                               ggplot2::aes(x = EffectSize, y = -log10(P),
                                            colour = threshold, text = ID)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.75) +
    # xlim(c(-6, 6)) +
    ggplot2::xlab("Effect Size") + ggplot2::ylab("-log10 p-value") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::geom_text(ggplot2::aes(label = volc.plot$lab),
                       hjust = 0, vjust = 0, colour = "dark grey", size = 2, nudge_x = 0.2)

  volc_DE_1 <- volc_DE_1 + ggplot2::geom_hline(yintercept = -log10(min(volc.plot$P[volc.plot$P.adj <= adj_p])), colour = "grey", linetype = "longdash") +
    ggplot2::geom_vline(xintercept = c(-1*effect_sizeThresh, effect_sizeThresh), colour = "grey",linetype = "longdash")
  # ggplot2::geom_hline(yintercept = -log10(adj_p), colour = "grey", linetype = "longdash") +


  if (makePlot)
    return(volc_DE_1)
  else
    return(volc.plot)

}







    deSeq_assay <- deSeq_assay[rownames(deSeq_assay) %in% genes_in, grepl(grps_in, colnames(deSeq_assay))]

  p <- heatmaply::heatmaply(deSeq_assay,
                       col_side_colors = colnames(deSeq_assay),
                       # row_side_colors = meta$phen,
                       hclust_method = "complete",
                       dist_method = "manhattan",
                       # fontsize_row = 5,
                       # fontsize_col = 5,
                       # row_text_angle = 45,
                       plot_method = "plotly"
  )



```








