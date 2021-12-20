#################################################
# differential accessibility at each time point #
#################################################
library(data.table)
library(ggplot2)
library(ggrepel)
library(limma)
library(preprocessCore)
source("/home/yuanh/programs/myscripts/NGS_scripts/ATAC_functions.R")

#################
# perform limma #
#################
sample_info <- fread("/home/yuanh/analysis/WI38_ATAC/analysis/atacseq/sample_info.txt")
counts <- readRDS("/home/yuanh/analysis/WI38_ATAC/analysis/atacseq/rdata/raw_counts.rds")

# remove TP3, because TP3 only exist for PDL, but not for hTERT
toremove1 <- which(sample_info$timepoint=="TP3")
toremove2 <- which(colSums(counts) < 1e6)
toremove3 <- which(sample_info$sample_name %in% c("hTERT_TP6_C", "PDL_TP6_C"))
toremove <- Reduce(union, list(toremove1, toremove2, toremove3))
counts <- counts[,-toremove]
coldata <- sample_info[-toremove,c("treatment","timepoint")]
coldata$sample <- paste(coldata$treatment, coldata$timepoint, sep="_")
coldata$treatment <- factor(coldata$treatment, levels = c("hTERT","PDL"))
coldata$timepoint <- factor(coldata$timepoint, levels = c("TP1","TP2","TP4","TP5","TP6","TP7"))

################
# perform voom #
################
design_coldata <- model.matrix(~treatment*timepoint, coldata)

# plot1
png("voom/voom_mean_variance_trend.png", 600, 300)
v <- voom(counts, design_coldata, plot=TRUE, normalize="quantile")
dev.off()
saveRDS(v$E, file="counts.rds")

# plot2
pdf("voom/count_distribution.pdf", 10, 6)
par(mar=c(10, 6, 3, 3))
boxplot(v$E, las=2, outline=F, ylab="log(normalized count)")
dev.off()

# perform limma
fit <- lmFit(v, design_coldata)
fit <- eBayes(fit)
saveRDS(fit, file="interaction_fit.rds")

# output result
output <- list()
for (i in c("TP2", "TP4", "TP5", "TP6", "TP7")) {
    coef_name <- paste0("treatmentPDL:timepoint",i)
    res <- topTable(fit, coef=coef_name, n=nrow(counts))
    output[[i]] <- res[rownames(counts), ]
}
saveRDS(output, file="interaction_res.rds")


##############
# plot DESeq #
##############
library(DESeq2)

dds <- readRDS("../atacseq/rdata/dds.rds")
cnts <- log2(counts(dds, normalize=T)+1)
# plot2
pdf("voom/deseq2_distribution.pdf", 10, 6)
par(mar=c(10, 6, 3, 3))
boxplot(cnts, las=2, outline=F, ylab="log(normalized count)")
dev.off()
