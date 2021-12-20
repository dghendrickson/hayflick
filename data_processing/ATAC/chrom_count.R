library(rtracklayer)
library(GenomicRanges)
source("~/hayflick/data/ATACseq/helper_functions/ATAC_functions.R")


atlas <- import.bed("~/hayflick/annotations/E017_25_imputed12marks_hg38lift_dense.bed")
sample_info <- fread("~/hayflick/data/ATACseq/sample_info.txt")   ###################### PATHS TO BAMS (not included--availible GSE175533)
counts <- countReads(atlas, sample_info$path, sample_info$sample_name) 
saveRDS(counts, "~/hayflick/data/ATACseq/chromHMM_counts_quantile_norm.rdscounts.rds")
