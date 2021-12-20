##################
# motif analysis #
##################
library(Biostrings)
library(data.table)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(dplyr)
library(colorspace)

rna <- read.table("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", sep="\t", header=T)
rna=rna[complete.cases(rna),]
rna_not_down_up=rna %>% filter(!(log2FoldChange < 0 & padj < 0.05 )) #### filter out TFs significantly going down

rownames(rna_not_down_up)=make.unique(rna_not_down_up$gene)

plot_motif_enrich <- function(genes) {
  tmp_genes <- genes[(genes %in% rna$gene) & 
                       (genes %in% atlas[atlas$distToNearest <= 50000]$gene_name) ]
  totest <- (res$TP7$logFC>0) & (res$TP7$adj.P.Val<0.001) & 
    (atlas$distToNearest < 50000) & (atlas$gene_name %in% tmp_genes)
  m_obs <- m[totest, ]
  m_exp <- m[!totest, ]
  # compute FIMO based motif stats
  toplot <- data.frame(obs_1 = apply(m_obs, 2, function(x) sum(x==1)),
                       obs_0 = apply(m_obs, 2, function(x) sum(x==0)),
                       exp_1 = apply(m_exp, 2, function(x) sum(x==1))/nrow(m_exp),
                       exp_0 = apply(m_exp, 2, function(x) sum(x==0))/nrow(m_exp))
  pvals <- apply(toplot, 1, function(x) {
    x <- unlist(x)
    return(binom.test(x[1:2], p=x[3], alternative="greater")$p.value)
  })
  toplot$pvals <- pvals
  toplot <- toplot[order(toplot$pvals),]
  toplot$x <- toplot$obs_1/nrow(m_obs)
  toplot$y <- -log10(toplot$pvals)
  toplot$TF <- rownames(toplot)
  
  return(toplot)
}



#rna$cut.50decile <-Hmisc::cut2(log2(rna$baseMean),g= 20)  # Assign every gene a percentile bin from DESEQ base mean qunatification of normalized counts 
#rna$cut.50decile=factor(rna$cut.50decile)
# load motifs



# add JUN_hocomoco motif


load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/JUN_motif_HOCOMOCO/motif_hits.rdt")
colnames(m) <- "JUN_hm"
m_jun <- m
# add FOXE1 motif
load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/FOXE1_motif/motif_hits.rdt")
m_FOXE1 <- m
load("~/hayflick/data_processing/ATAC/motif_process/fimo/motif_hits.rdt")
m <- m[, colnames(m) %in% (rna_not_down_up[[1]])] # filter for only expressed TFs
m <- cbind(m, m_jun, m_FOXE1)


res <- readRDS("~/hayflick/data/ATACseq/limmaVoom_DE_res.rds")
atlas <- readRDS("~/hayflick/data/ATACseq/peak_atlas.rds")

# RS top 1000 list
genes=rna[complete.cases(rna),]
genes <- genes[genes$padj < 0.01 & genes$log2FoldChange>.5,]
genes <- genes[order(genes$log2FoldChange, decreasing = TRUE),]
#genes <- as.character(genes[1:1000, "gene"])
genes <- as.character(genes[, "gene"])


#genes=sen.features.spec
#genes=con.features.spec

toplot <- plot_motif_enrich(genes)

toplot$padj <- p.adjust(toplot$pvals, method="BH")
toplot$color <- ifelse(toplot$padj<0.001, "red", "grey")
toplot$y <- -log10(toplot$padj)

p <- ggplot(toplot) + geom_point(aes(x, y, color=color)) +
  geom_text_repel(data=toplot[toplot$padj<0.001,],max.overlaps = 30, aes(x, y, label = TF)) +
  scale_color_manual(values=c("red"="red", "grey"="grey")) +  
  geom_hline(yintercept=-log10(0.001), linetype=2) +
  labs(x="percent_hits", y="-log10(padj)", title="motif enrichment (padj < 0.001)") + theme_classic()
pdf("~/hayflick/RS_enrichment_0.001.pdf")
print(p)
dev.off()

################################################### RIDGE REGRESSION

library(glmnet)
library(PRROC)
require(doMC)
library(ggplot2)
registerDoMC(cores = 5)

res <- readRDS("interaction_res.rds")

###############
# load motifs #
###############
#rna <- readRDS("/home/yuanh/analysis/WI38_ATAC/analysis/rnaseq/comparison2/res.rds") # perform RNA-seq myself
# add JUN_hocomoco motif



load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/JUN_motif_HOCOMOCO/motif_hits.rdt")
colnames(m) <- "JUN_hm"
m_jun <- m
# add FOXE1 motif
load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/FOXE1_motif/motif_hits.rdt")
m_FOXE1 <- m
load("~/hayflick/data_processing/ATAC/motif_process/fimo/motif_hits.rdt")
m <- m[, colnames(m) %in% (rna_not_down_up[[1]])] # filter for only expressed TFs
m <- cbind(m, m_jun, m_FOXE1)


select <- abs(res$TP7$logFC)>0 & res$TP7$adj.P.Val<0.001
x <- m[select, ]
y <- 1*(res$TP7$logFC[select]>0)
all_ids <- 1:nrow(x)

set.seed(10)
test_ids <- sample(all_ids, length(all_ids)/5)
train_ids <- all_ids[!all_ids%in%test_ids]

# train elastic model on train
start.time <- proc.time()
cvfit = cv.glmnet(x[train_ids, ], y[train_ids], family="binomial", alpha=0, nfolds=5, parallel=T) # ridge classification
y_pred <- predict(cvfit, newx = x[test_ids, ], s = "lambda.min")
print(proc.time()-start.time)
roc <- roc.curve(y_pred[y[test_ids]==1], y_pred[y[test_ids]==0], curve=T)
#pdf("ridge/test_roc.pdf")
plot(roc)
dev.off()
saveRDS(cvfit, file="ridge/cvfit.rds")

# look at coefficients
cvfit <- readRDS("ridge/cvfit.rds")
set.seed(20)
coef_list <- list()
for (i in 1:10) {
  ids <- sample(all_ids, length(all_ids), replace=T)
  fit <- glmnet(x[ids, ], y[ids], family="binomial", alpha = 0, lambda = cvfit$lambda.min)
  coefs <- coef(fit)
  coef_list[[i]] <- coefs[2:nrow(coefs),1]
  print(i)
}
coefs <- do.call(cbind, coef_list)
toplot <- data.frame(TF=rownames(coefs),
                     coef_mean=rowMeans(coefs),
                     coef_se=apply(coefs, 1, sd)/sqrt(ncol(coefs)))
toplot <- toplot[order(toplot$coef_mean, decreasing=T)[1:20], ]
toplot$TF <- factor(toplot$TF, levels=as.character(toplot$TF))

#pdf("~/hayflick/coef.pdf")
ggplot(toplot) +
  geom_bar( aes(x=TF, y=coef_mean), stat="identity", fill="gray") +
  geom_errorbar(aes(x=TF, ymin=coef_mean-coef_se, ymax=coef_mean+coef_se), width=0.4, size=1.3)  +
  xlab("Transcription Factor") + ylab("Ridge Regression Coefficient") + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


############################ LISA

LISA_new_up=read.delim("~/hayflick/data/RNAseq/LISA/results_RS_ALLsigUp.tsv", stringsAsFactors = FALSE)
unique(LISA_new_up$cell_type)

unique(LISA_new_up$cell_type)

LISA_new_up_fibroblast=LISA_new_up %>% filter(cell_type=="Fibroblast")

unique(LISA_new_up_fibroblast$factor)


LISA_new_down=read.delim("~/hayflick/data/RNAseq/LISA/results_RS_ALLsig_DOWN.tsv", stringsAsFactors = FALSE)

Lisa_AU=read.delim("~/hayflick/data/RNAseq/LISA/astro_up.tsv", stringsAsFactors = FALSE)
LISA_AU_df=Lisa_AU %>% dplyr::group_by(factor) %>% dplyr::summarize(summary_p_value_mean = median(summary_p_value),sd=sd(summary_p_value)) 

Lisa_AD=read.delim("~/hayflick/data/RNAseq/LISA/astro_down.tsv", stringsAsFactors = FALSE)
LISA_AD_df=Lisa_AD %>% dplyr::group_by(factor) %>% dplyr::summarize(summary_p_value_mean = median(summary_p_value),sd=sd(summary_p_value)) 


############
RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = FALSE)
RS_TPMs=read.delim("~/hayflick/data/RNAseq/TPM_tables/RS_batch_tpms.txt", stringsAsFactors = FALSE)
rownames(RS_TPMs)=make.unique(RS_TPMs$gene)
RS_TPMs=RS_TPMs[-1]
sen_TPMS_exp=RS_TPMs[!rowSums(RS_TPMs >5),]
bad_genes=rownames(sen_TPMS_exp)
sen_TPMS_exp=RS_TPMs[!rownames(RS_TPMs) %in% bad_genes,]

#######
#LISA_new=LISA_new %>% filter(ChIP.seq_p_value < 0.00001)

#LISA_new_df_up=LISA_new_up %>% dplyr::group_by(factor) %>% dplyr::summarize(chipP_mean = median(ChIP.seq_p_value),sd=sd( ChIP.seq_p_value)) 

LISA_new_df_up=LISA_new_up %>% dplyr::group_by(factor) %>% dplyr::summarize(summary_p_value_mean = median(summary_p_value),sd=sd(summary_p_value)) 
LISA_new_df_down=LISA_new_down %>% dplyr::group_by(factor) %>% dplyr::summarize(summary_p_value_mean = median(summary_p_value),sd=sd( summary_p_value))


Lisa_RS_merge=merge(LISA_new_df_up,LISA_new_df_down, by="factor")
Lisa_RS_merge=Lisa_RS_merge[Lisa_RS_merge$factor %in% rownames(sen_TPMS_exp),]


Lisa_RS_merge$log10up=-log10(Lisa_RS_merge$summary_p_value_mean.x)
Lisa_RS_merge$log10down=-log10(Lisa_RS_merge$summary_p_value_mean.y)



get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



Lisa_RS_merge$density <- get_density(Lisa_RS_merge$log10up,Lisa_RS_merge$log10down)


Lisa_up_RS_genes=Lisa_RS_merge %>% filter(log10up >5 & log10down < 30) %>% dplyr::select("factor")
Lisa_up_RS_genes=Lisa_up_RS_genes[!grepl("ZNF", Lisa_up_RS_genes$factor),]

Lisa_DON_genes=Lisa_RS_merge %>% filter(log10down >50 & log10up < 10) %>% dplyr::select("factor")
Lisa_DON_genes=Lisa_DON_genes[!grepl("ZNF", Lisa_DON_genes$factor),]

Lisa_both=Lisa_RS_merge %>% filter(log10down >30 & log10up > 10) %>% dplyr::select("factor")
Lisa_both=Lisa_both$factor


Lisa_RS_merge$diff_rank=Lisa_RS_merge$log10up-Lisa_RS_merge$log10down

pdf("/Volumes/GoogleDrive/Shared drives/Hayflick_paper/Figures/fig_chunks/TRD_RS_LISA_UP.down.pdf",height=6, width=6)

ggplot(Lisa_RS_merge) +
  geom_point(aes(x=log10up,y=log10down,color=density),size =2,alpha = .25) +
  scale_color_continuous_sequential("Dark Mint",end = .7, rev = FALSE) +
  #scale_x_continuous(limits=c(0,15)) + 
  scale_y_continuous(limits=c(0,90)) +
  
  geom_text(data=subset(Lisa_RS_merge, Lisa_RS_merge$factor %in% Lisa_up_RS_genes),aes(log10up,log10down,label=factor),size=2.5) +
  geom_text(data=subset(Lisa_RS_merge, Lisa_RS_merge$factor %in% Lisa_DON_genes),aes(log10up,log10down,label=factor),size=2.5) +
  geom_text(data=subset(Lisa_RS_merge, Lisa_RS_merge$factor %in% Lisa_both),aes(log10up,log10down,label=factor),size=2) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
dev.off()



########################

Lisa_A_merge=merge(LISA_AU_df,LISA_AD_df, by="factor")


Lisa_A_merge$log10up=-log10(Lisa_A_merge$summary_p_value_mean.x)
Lisa_A_merge$log10down=-log10(Lisa_A_merge$summary_p_value_mean.y)

Lisa_A_merge$density <- get_density(Lisa_A_merge$log10up,Lisa_A_merge$log10down)


Lisa_A_sen_genes=Lisa_A_merge %>% filter(log10up >5 & log10down < 1) %>% dplyr::select("factor")
Lisa_A_sen_genes=Lisa_A_sen_genes[!grepl("ZNF", Lisa_A_merge$factor),]

Lisa_A_DON_genes=Lisa_A_merge %>% filter(log10down >60 & log10up < 10) %>% dplyr::select("factor")
Lisa_A_DON_genes=Lisa_A_DON_genes[!grepl("ZNF", Lisa_A_DON_genes$factor),]

Lisa_both=Lisa_A_merge %>% filter(log10down >30 & log10up > 10) %>% dplyr::select("factor")
Lisa_both=Lisa_both$factor


RS_dobbla=intersect(Lisa_A_sen_genes,Lisa_up_RS_genes)

Lisa_A_merge$diff_rank=Lisa_A_merge$log10up-Lisa_A_merge$log10down



##############################


Lisa_RS_small=Lisa_RS_merge %>% dplyr::select(c(factor,diff_rank))

Lisa_A_small=Lisa_A_merge %>% dplyr::select(c(factor,diff_rank))

small_rank_merge=merge(Lisa_RS_small,Lisa_A_small,by="factor")


small_rank_merge$rank_rs=rank(small_rank_merge$diff_rank.x,)

small_rank_merge$rank_A=rank(small_rank_merge$diff_rank.y)


small_rank_merge$density=get_density(small_rank_merge$rank_rs,small_rank_merge$rank_A)


same_genes=small_rank_merge %>% filter(rank_rs >500 & rank_A > 200) %>% dplyr::select("factor")
#Lisa_A_sen_genes=Lisa_A_sen_genes[!grepl("ZNF", Lisa_A_merge$factor),]
same_genes=same_genes$factor

#small_rank_merge_pv=  small_rank_merge %>% dplyr::select(c(factor,rank_rs,rank_A))   %>%
# pivot_longer(!c(factor), names_to="group",values_to="rank")




small_rank_merge = within(small_rank_merge, {
  concordance = ifelse(small_rank_merge$factor %in% same_genes , "concordant", "discordant")
})


small_rank_merge = within(small_rank_merge, {
  type = ifelse(grepl("ZNF|ZBT|ZSCAN",small_rank_merge$factor) , "ZNF", "TF")
})


same_genes=same_genes[!grepl("ZNF|ZBT|ZSCAN", same_genes)]

pdf("/Volumes/GoogleDrive/Shared drives/Hayflick_paper/Figures/fig_chunks/A_vs_RS_LISA.rank.pdf",height=6, width=6)

ggplot(small_rank_merge) +
  geom_point(aes(x=rank_rs,y=rank_A,color=factor(concordance),shape=type),size =1.5,alpha = .9) +
  #geom_point(aes(x=PDL50,y=avg50,color=density),size =2,alpha = .4) + 
  # scale_color_continuous_sequential("Dark Mint",end = .7, rev = FALSE) +
  scale_color_discrete_qualitative("Dark 3") +
  geom_text(data=subset(small_rank_merge, small_rank_merge$factor %in% same_genes),
            aes(rank_rs,rank_A,label=factor),size=2,nudge_x = 20, nudge_y = 5) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
dev.off()





writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_6/sessionInfo.txt")








