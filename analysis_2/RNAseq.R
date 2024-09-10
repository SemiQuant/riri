require(tidyverse)
require(DESeq2)
require(heatmaply)
require(factoextra)
require(clusterProfiler)
require(KEGGREST)


# setwd("analysis_2/")


meta <- readxl::read_excel("metadata.xlsx")
meta <- meta %>% 
  filter(!is.na(Group))


mn_fld <- "/Users/semiquant/Downloads/RNAseq_pipeline_20240906_223306/"
files_in <- list.files(paste0(mn_fld), "*featCount.counts$",
                       full.names = T, recursive = T)

genes <- read_tsv(files_in[1], skip = 1, col_select = 1)


counts <- files_in %>%
  map(read_tsv, skip = 1, col_select = 7, show_col_types = F) %>%
  purrr::reduce(cbind)

colnames(counts) <- dirname(colnames(counts))
counts <- cbind(genes, counts)


counts <- counts %>% 
  select(one_of(c("Geneid", meta$Filename)))

# filter those without any counts but do it after selecting groups etc
counts <- counts %>%
  filter(rowSums(select(., -1), na.rm = TRUE) >= 10)


gen_names <- counts$Geneid
rownames(counts) <- gen_names

counts <- counts %>%
  dplyr::select(-c(1))
genes <- NULL



meta <- meta %>% 
  filter(Filename %in% colnames(counts))

meta <- meta[match(meta$Filename, colnames(counts)), ]
table(meta$Filename == colnames(counts))
meta <- meta %>%
  mutate(Group = as.factor(Group)) %>% 
  mutate(Timepoint_min = as.factor(Timepoint_min)) %>% 
  mutate(Rep = as.factor(Rep))

dds_1 <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                # tidy = T,
                                design = ~ Timepoint_min) #Group + Batch
dds_1 <- estimateSizeFactors(dds_1)
dds_1 <- DESeq(dds_1, parallel = T)
res_1 <- DESeq2::results(dds_1)

vst_1_blind <- vst(dds_1, blind = T)
vst_1_blind_cor <- cor(assay(vst_1_blind))







# Heatmap of all transcripts color by groups ------------------------------------------------
hmap_QC_1 <- heatmaply(vst_1_blind_cor,
                       # symm = T,
                       col_side_colors = meta$Batch,
                       row_side_colors = meta$Timepoint_min,
                       hclust_method = "complete",
                       dist_method = "manhattan",
                       fontsize_row = 5,
                       fontsize_col = 5,
                       # row_text_angle = 45,
                       plot_method = "plotly"
)
# hmap_QC_1

vst_1_blind_t <- t(assay(vst_1_blind))
res.pca <- prcomp(vst_1_blind_t,  scale = T)
fviz_pca_ind(res.pca, 
             label="quanti.sup",
             habillage=meta$Batch,
             addEllipses=TRUE, ellipse.level=0.95)

"these are the two outliers"
DT::datatable(meta %>% 
                filter(Filename %in% c("NC_1", "NC_2"))
)

"with those removed"
meta_sub <- meta %>% 
    filter(!Filename %in% c("NC_1", "NC_2"))

vst_1_blind_t_sub <- vst_1_blind_t[!rownames(vst_1_blind_t) %in% c("NC_1", "NC_2"),]
res.pca_sub <- prcomp(vst_1_blind_t_sub,  scale = F)
fviz_pca_ind(res.pca_sub, 
             label="none",
             habillage=meta_sub$Batch,
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(res.pca_sub, 
             label="none",
             habillage=meta_sub$Timepoint_min,
             addEllipses=TRUE, ellipse.level=0.95)







vol_plot <- function(deSeq_results = NULL, adj_p = 0.05, num_annotate = 10, effect_sizeThresh = 2, makePlot = T){
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



resultsNames(dds_1)
res_2 <- DESeq2::results(dds_1, contrast = c("Timepoint_min","30","0"))
v1 <- ggplotly(
  vol_plot(deSeq_results = res_2, adj_p = 0.05, num_annotate = 10, effect_sizeThresh = 2, makePlot = T)
)

res_2b <- DESeq2::results(dds_1, contrast =c("Timepoint_min","360","0"))
v2 <- ggplotly(
  vol_plot(deSeq_results = res_2b, adj_p = 0.05, num_annotate = 10, effect_sizeThresh = 2, makePlot = T)
)

print("0 vs 30 | 0 vs 360")
subplot(v1, v2, shareY = T)



res_1a <- data.frame(ID = rownames(res_2),
                     rna_30 = res_2$log2FoldChange,
                     threshold_30 = res_2$padj <= 0.05
)

res_2a <- data.frame(ID = rownames(res_2b),
                     rna_360 = res_2b$log2FoldChange,
                     threshold_360 = res_2b$padj <= 0.05
)


rna_rib <- res_1a %>% 
  left_join(res_2a) %>% 
  drop_na


rna_rib$threshold <- ifelse(
  rna_rib$threshold_360 & rna_rib$threshold_30, "Both",
  ifelse(rna_rib$threshold_360, "360m",
         ifelse(rna_rib$threshold_30, "30m",
                "None"
         )))

rna_rib %>% 
  plot_ly(x = ~rna_30, y = ~rna_360, color = ~threshold, text = ~ID, legendgroup = ~ID,
          hovertemplate = paste(
            "<b>Gene: </b>%{text}<br>",
            "<b>30: </b>%{x}<br>",
            "<b>360: </b>%{y}<br>"
          )
  ) %>% 
  layout(legend=list(title=list(text='<b> Significant </b>')))




rna_rib %>% 
  data.frame() %>% 
  datatable(extensions = 'Buttons',
            options = list(
              paging = TRUE,
              searching = TRUE,
              fixedColumns = TRUE,
              autoWidth = TRUE,
              ordering = TRUE,
              dom = 'tBp',
              buttons = c('copy', 'csv', 'excel')
            ))










(kegg_enrich_both <- enrichKEGG(gene         = rna_rib[rna_rib$threshold == "Both", ]$ID,
                                organism     = 'mtu',
                                keyType      = 'kegg',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
)
# kegg_enrich_both

(kegg_enrich_30 <- enrichKEGG(gene         = rna_rib[rna_rib$threshold_30, ]$ID,
                              organism     = 'mtu',
                              keyType      = 'kegg',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)
)
# dotplot(kegg_enrich_30, showCategory = 10)

(kegg_enrich_360 <- enrichKEGG(gene         = rna_rib[rna_rib$threshold_360, ]$ID,
                               organism     = 'mtu',
                               keyType      = 'kegg',
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)
)
dotplot(kegg_enrich_360, showCategory = 10)









