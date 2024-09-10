setwd("/Users/SemiQuant/Bioinformatics/Projects/riri/analysis/testing")
require(Biostrings)
require(GenomeInfoDb)
require(rtracklayer)

download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz")
dna <- readDNAStringSet("GCF_000195955.2_ASM19595v2_genomic.fna.gz")

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))

chrominfo <- getChromInfoFromNCBI("GCF_000195955.2")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences. An alternative would be to rename them to
### chrominfo[ , "SequenceName"] but these names are VERY ugly (e.g.
### "ScRZk8e_1;HRSCAF=1").
names(dna) <- expected_RefSeqAccn

### Export as 2bit.
export.2bit(dna, "GCF_000195955.2.sorted.2bit")


fileConn <- file("BSgenome.Mtuberculosis.NCBI.ASM19595v2-seed")
writeLines(
'Package: BSgenome.Mtuberculosis.NCBI.ASM19595v2
Title: Full genome sequences for Mycobacterium tuberculosis H37Rv (high GC Gram+)
Description: Full genome sequences for Mycobacterium tuberculosis as provided by NCBI (ASM19595v2, 2013-02-01)
Version: 1.0.0
organism: Mycobacterium tuberculosis
common_name: M. tuberculosis
provider: NCBI
source_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/
SrcDataFiles: GCF_000195955.2_ASM19595v2_genomic.fna.gz from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/
organism_biocview: mycobacterium_tuberculosis
genome: ASM19595v2
BSgenomeObjname: Mtuberculosis
seqfiles_suffix: .fna.gz
PkgExamples: genome[["NC_000962.3"]]
seqfile_name: GCF_000195955.2.sorted.2bit
seqs_srcdir: /Users/SemiQuant/Bioinformatics/Projects/riri/analysis/testing/
circ_seqs: character(0)
release_date: 2020/09/02', fileConn)

close(fileConn)

forgeBSgenomeDataPkg("BSgenome.Mtuberculosis.NCBI.ASM19595v2-seed")

install.packages("BSgenome.Mtuberculosis.NCBI.ASM19595v2", 
                 repos = NULL, 
                 type = "source")

require(BSgenome.Mtuberculosis.NCBI.ASM19595v2)
genome <- BSgenome.Mtuberculosis.NCBI.ASM19595v2




