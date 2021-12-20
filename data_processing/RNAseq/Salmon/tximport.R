
library(tximport)
library(readr)
library(dplyr)

##TXIMPORT

###meta data
metadata=read.delim("~/hayflick/data_processing/RNAseq/DEseq/RS_column_data.txt",stringsAsFactors = FALSE)



##### IMPORT SALMON output file list----not provided--but this is what we did

setwd("path/to/salmon/output")

# 
file.paths=read.delim("file.paths.tab", stringsAsFactors = FALSE,header = FALSE)
file.paths=file.paths$V1

################ transcript to whole locus conversion file
gene_names <- read.delim(file="tx2gene_hg38_cellRanger3_V2.txt", header=TRUE, stringsAsFactors=FALSE)
t2g <- data.frame(target_id=gene_names$Name, gene_id=gene_names$gene_name)
t2g$gene_id=as.character(t2g$gene_id)
t2g$target_id=as.character(t2g$target_id)
t2g=t2g[,c(1,2)]

##########run tximport
txi_2x <- tximport(e_fil, type = "salmon",txOut = TRUE ,tx2gene = t2g)

########
### PULL OUT TPM table and make data frame for dyplyr filtering later

abundance_raw=txi_2$abundance

### ASSIGN SAMPLE NAMES to columns
colnames(abundance_raw)=s2c$sample
abundance_raw_df=as.data.frame(abundance_raw)

### PULL OUT raw counts for COMBAT AND DESEQ
counts_raw=txi_2$counts


### ASSIGN SAMPLE NAMES to columns
colnames(counts_raw)=s2c$sample
#counts_raw=as.data.frame(counts_raw)


#################   Abundance (TPMs) and counts output from this script can be found in "~/hayflick/data/RNAseq/DEseq/"
#####



