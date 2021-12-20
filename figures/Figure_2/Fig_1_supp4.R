library(BiocManager)

library(readr)
library(dplyr)

library(ggplot2)
library(cowplot)

library(tidyr)
library(reshape2)
library(fgsea)
library(tximport)
library(pheatmap)
library(yarrr)
library(sva)
library(stringr)
setwd("/home/dgh/hayflick/figures/Figure_1/")

#######################
##FIG sup2  GSEA by individual PDL
##############



TERT=read.delim("~/hayflick/data/RNAseq/TPM_tables/TERT_batch_tpms.txt", stringsAsFactors = FALSE) ##get tetrt info and avg. replicates

TERT_pivot=TERT %>% tidyr::pivot_longer(-gene) 

TERT_pivot$value=TERT_pivot$value +1 #add pseudocount

sample_tracking=str_split_fixed(TERT_pivot$name,"_",3)

TERT_pivot$group=sample_tracking[,2]

TERT_pivot_avg=TERT_pivot %>%
  group_by(group,gene) %>%
  dplyr::summarize(expression_mean = mean(value)) %>% 
tidyr::pivot_wider(id_cols=gene,names_from=group, values_from=expression_mean)

TERT_pivot_avg_fc=log2(TERT_pivot_avg[,3:7]/TERT_pivot_avg$TP1)
TERT_pivot_avg_fc=cbind(TERT_pivot_avg$gene,TERT_pivot_avg_fc)
colnames(TERT_pivot_avg_fc)[1]="gene"

avg.fc=read.delim("~/hayflick/data/RNAseq/TPM_tables/idv.TP_GSEA.txt", stringsAsFactors = FALSE) # read in replicate avg data for rest of samples


avg.fc=merge(TERT_pivot_avg_fc,avg.fc, by="gene")


rownames(avg.fc)=make.unique(avg.fc$gene)

avg.fc=avg.fc[-1]

fullGSEA <- vector("list", length=dim(avg.fc)[2])

rank_list<- list()

for ( i in 1:(dim(avg.fc)[2])){
  rank_list[[i]]= setNames(avg.fc[,i], rownames(avg.fc))
}

names(rank_list)=colnames(avg.fc[1:dim(avg.fc)[2]])

names(fullGSEA) <-colnames(avg.fc)

### Loop for running GSEA ON ALL CONDITIONS IN rank list LIST and GSEA list

pathways2 <- gmtPathways("~/hayflick/annotations/h.all.v7.0.symbols_.gmt") 


for (set in 1:length((rank_list))) {
  print(names(rank_list)[set])
  ranks <- rank_list[[set]]
  
  fgseaRes <- fgseaMultilevel(pathways2, ranks ,minSize=10, maxSize=1500,eps = 0)
  fullGSEA[[set]] <- fgseaRes
}  


pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
#pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
betaGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 6)]) ##col 5 is NES

###this builds the whole matrix, change NAs to Pval of 1 and NES of zero

colnames(pvalGSEA) <- c("pathway", paste(names(fullGSEA), "pval", sep="_"))
colnames(betaGSEA) <- c("pathway", paste(names(fullGSEA), "NES", sep="_"))

pvalGSEA[is.na(pvalGSEA)] <- 1
betaGSEA[is.na(betaGSEA)] <- 0


pvalGSEA[,2:ncol(pvalGSEA)] <- -log10(pvalGSEA[,2:ncol(pvalGSEA)]) * sign(betaGSEA[,2:ncol(betaGSEA)])

###convert to numeric matrices with named rows for easy clustering
pvalGSEA2 <- as.matrix(pvalGSEA[,2:ncol(pvalGSEA)])
rownames(pvalGSEA2) <- pvalGSEA[,1]


betaGSEA2 <- as.matrix(betaGSEA[,2:ncol(betaGSEA)])
rownames(betaGSEA2) <- betaGSEA[,1]


pvalGSEA2=pvalGSEA2[complete.cases(pvalGSEA2),]
betaGSEA2=betaGSEA2[rownames(betaGSEA2) %in% rownames(pvalGSEA2),]


################################## PVALUE CUT OFF FOR PRUNING
pcut <- -log10(.01) #p value threshold

library(reshape2)
tmp <- reshape::melt(apply(pvalGSEA2, 1, function(x){length(which(abs(x)>=pcut))}))
##i verified that the rows are in the same order as the starting matrix, so I'm just going to do direct indexing removal

pvalGSEA2 <- pvalGSEA2[-which(tmp$value==0),]
betaGSEA2 <- betaGSEA2[-which(tmp$value==0),]

rwb<-colorRampPalette(c("steel blue","white","tomato"))
rwb=rwb(100)


pheatmap(pvalGSEA2,
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         gaps_col = c(5,14,19),
         breaks = seq(-4, 4, length.out=101)
)




## ##########################Figure S2 B
###########    retrieve leading edge gened from GSEA output 

temp=data.frame()
for (i in names(fullGSEA[6:14])) {  #just RS
  gs.dat=as.data.frame(fullGSEA[[i]])
  igs.dat=gs.dat %>% filter(padj < .01 )# %>% filter(NES > 0 )
  for (j in seq(1:length(gs.dat$pathway))){
    temp2=as.data.frame(gs.dat[j,8])
    temp2$pathway=gs.dat[j,1]
    temp2=temp2[,c(2,1)]
    temp2$tp=names(fullGSEA[i])
    colnames(temp2)=c("path","gene","PDL")
    temp=rbind(temp,temp2)
  }
  
} 

deduped.data <- unique(temp[,1:2])


EMT=deduped.data %>% filter(path=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")  ## pull leading edge genes from EMT annotation

ids=EMT$gene

small=TPM_FC_MERGE[TPM_FC_MERGE$Row.names %in% ids,] ### TPM fold change matrix from  matrix from Fig 1 script

    rownames(small)=small$Row.names
    small=small[-1]


        pheatmap(small,
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         show_colnames = TRUE,

         #  gaps_col = c(15,42,57),
       
         breaks = seq(-1.5, 1.5, length.out=101)
        )




