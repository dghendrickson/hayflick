library(DESeq2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
library(ggplot2)
library(ggrepel)
source("~/hayflick/data_processing/ATAC/helper.R")


############
# run FIMO #
############
# create fasta file to run fimo, 500bp center at each peak
dir.create("fimo")
atlas <- import.bed("~/hayflick/data_processing/ATAC/motif_process/atlas.bed") ##import peak atlas as bed file
org <- Hsapiens
atlas <- resize(atlas, width=500, fix="center")
seqs <- get.seqs(org, atlas, no.cores = 1)
writeXStringSet(seqs, "~/hayflick/data_processing/ATAC/motif_process/fimo/atlas.fasta", format="fasta") ###out put as fasta for FIMO to find motifs in peaks

#run FIMO MUST HAVE MEME INSTALLED (USED VERSION 5.4.0)
system("/home/dgh/meme-5.4.0/src/fimo ~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/CIS-BP/Homo_sapiens.meme ~/hayflick/data_processing/ATAC/motif_process/fimo/atlas.fasta")
##################
# process result #
##################
# summarize fimo result
# 1) consider only expressed TFs
system("grep MOTIF ~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/CIS-BP/Homo_sapiens.meme > ~/hayflick/data_processing/ATAC/motif_process/fimo/motif_names.txt")
tfs <- read.table("~/hayflick/data_processing/ATAC/motif_process/fimo/motif_names.txt", stringsAsFactors=F)
tfs$TF <- gsub("_.*|\\(|\\)","",tfs$V3)
rownames(tfs) <- tfs$V2

# 2) process fimo result
library(Biostrings)
library(ggplot2)
library(ggrepel)
library(data.table)
test <- readDNAStringSet("~/hayflick/data_processing/ATAC/motif_process/fimo/atlas.fasta","fasta")
#expressed_TFs <- read.table("fimo/expressed_TFs.txt", stringsAsFactors = F)
#rownames(expressed_TFs) <- expressed_TFs$V2
a <- fread("~/hayflick/data_processing/ATAC/motif_process/fimo_out_han/fimo.tsv", header=T, comment.char= "#")
motifs <- unique(a$motif_id)
colnames(a) <- gsub("-","",colnames(a))
temp <- split(a, f=factor(a$motif_id))

# motif score matrix
motif_p_m <- matrix(1, length(test), length(motifs), dimnames=list(names(test), motifs))
for (i in 1:length(temp)){
  motif <- names(temp)[i]
  pvalue_tmp <- tapply(temp[[i]]$pvalue, temp[[i]]$sequence_name, min) # for each peak, find the min pvalue
  motif_p_m[names(pvalue_tmp), motif] <- pvalue_tmp
  print(i)
}

# make it binary
m <- motif_p_m
colnames(m) <- tfs[colnames(m),"TF"]
m <- m[,order(colnames(m))]
#m <- m[, colnames(m) %in% expressed_TFs$V2]
m <- 1 * (m<1e-4)
save(motif_p_m, m, file="analysis/fimo/motif_hits.rdt")
