# require(plotly)
# require(tidyverse)
# 
# files <- list.files("/Users/SemiQuant/Downloads/riboDelete/out_dir_test", full.names = T)
# 
# 
# data <- files %>%
#   map(read_tsv, col_names = F, show_col_types = FALSE) %>%
#   reduce(cbind) 
# colnames(data) <- gsub(".txt", "", basename(files))
# 
# data %>% 
#   mutate(n = 1:nrow(data)) %>% 
#   pivot_longer(!n, names_to = "Gene", values_to = "count") %>% 
#   plot_ly(x = ~n, y = ~count, type = 'scatter', mode = 'lines', color = ~Gene)
# 
# 
# 
# t() %>% 
#   data.frame() %>% 
#   
#   




require(RiboDiPA)
bam_path <- list.files(pattern = ".bam$", "/Users/SemiQuant/Downloads/riboDelete/riboSeq_R1_riri", full.names = T)
gtf <- "/Users/SemiQuant/Bioinformatics/Projects/riri/references/NC_000962.gff"

classlabel <- data.frame(condition = c("mutant",  "wildtype"), comparison = c(2, 1))

# https://bioconductor.org/packages/devel/bioc/vignettes/RiboDiPA/inst/doc/RiboDiPA.html#p-site-mapping
## Perform individual P-site mapping procedure
data.psite <- psiteMapping(bam_file_list = bam_path, 
                           gtf_file = gtf, psite.mapping = "auto", cores = 5)

## P-site mapping offset rule generated
data.psite$psite.mapping

# ## Use user specified psite mapping offset rule
# offsets <- cbind(qwidth = c(28, 29, 30, 31, 32), 
#                  psite = c(18, 18, 18, 19, 19))
# data.psite2 <- psiteMapping(bam_path[1:4], bam_path[5], 
#                             psite.mapping = offsets, cores = 2)


## Merge the P-site data into bins with a fixed or an adaptive width
data.binned <- dataBinning(data = data.psite$coverage, bin.width = 0, 
                           zero.omit = FALSE, bin.from.5UTR = TRUE, cores = 2)

## Merge the P-site data on each codon
data.codon <- dataBinning(data = data.psite$coverage, bin.width = 1, 
                          zero.omit = FALSE, bin.from.5UTR = TRUE, cores = 2)

## Merge the P-site data on each exon and perform differential pattern analysis
result.exon <- diffPatternTestExon(psitemap = data.psite, 
                                   classlabel = classlabel, method = c('gtxr', 'qvalue'))

## Perform differential pattern analysis
result.pst <- diffPatternTest(data = data.binned, 
                              classlabel = classlabel, method=c('gtxr', 'qvalue'))


## Plot ribosome per nucleotide tracks of specified genes.
plotTrack(data = data.psite, genes.list = c("YDR050C", "YDR064W"),
          replicates = NULL, exons = FALSE)

## Plot binned ribosome tracks of siginificant genes: YDR086C and YDR210W.
## you can specify the thrshold to redefine the significant level
plotTest(result = result.pst, genes.list = NULL, threshold = 0.05) 