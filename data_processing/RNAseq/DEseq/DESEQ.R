library(DESeq2)
library(apeglm)
library(dplyr)

### CREATe col data object for deseq 
## 
##
coldata=read.delim("~/hayflick/data_processing/RNAseq/DEseq/RS_column_data.txt", stringsAsFactors = FALSE)
counts=read.delim("~/hayflick/data/RNAseq/count_tables/RS_raw_counts.txt", stringsAsFactors = FALSE)

rownames(coldata)=colnames(counts)

###round counts to integers for DESEQ

counts<- round(counts)


############### create DESEQ object; specify model from coldata. This one takes into account both cell line (proportion_timecourse (time) and Batch for late sen samples PDL45,
# PDL52, and PDL53)


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ proportion_timecourse  + batch)


#### create deseq results--all results from specified model in this object

dds=DESeq(dds,test="Wald")

###### run this command to see comparisons made
resultsNames(dds)

######results with shrinkage; 
#####pull out specifically batch corrected RS 

resLFC <- lfcShrink(dds, coef="proportion_timecourse", type="apeglm")

RS_De=(as.data.frame(resLFC))
RS_De$id=rownames(RS_De)

"~/hayflick/data/RNAseq/DESEQ_output/  "

#####################repeat for other conditions i.e. cell density, radiation--  this was script used to generate DESEQ results here:
#####################  "~/hayflick/data/RNAseq/DESEQ_output/" 

