library(BiocManager)
library(DESeq2)
library(tximport)
library(readr)
library(dplyr)
library(apeglm)
library(ggplot2)
library(cowplot)
library(sva)
library(tidyr)
library(reshape2)
library(fgsea)
library(tximport)
library(pheatmap)
library(yarrr)
library(sva)
library(stringr)
library()

setwd("/home/dgh/hayflick/figures/Figure_1/")

############################
######################  1C

betaG=read.delim("/home/dgh/hayflick/figures/Figure_1/betaGal.txt", stringsAsFactors = FALSE)
colnames(betaG)

yarrr::pirateplot(formula = percent ~ TP.char + strain,    # DV = height, IV1 = sex, IV2 = headband
                  data = betaG,           
                  theme =4, inf.f.o = 0, avg.line.o = 0,point.o = 1,
                  main = "Pirate Heights",
                  pal = "gray")


############################
######################  1D


##############
####################### 1E

############# Batch correct TPMs for RS
meta_data_RS=read.delim("~/hayflick/data_processing/RNAseq/DEseq/RS_column_data.txt", stringsAsFactors = FALSE)


##### set up batch to correct from S2c
  batch=meta_data_RS$batch

#####set up covartiates i.e. timepoints so COmbat knows what should and shoudlnt look similar
group=meta_data_RS$TP

###get raw counts from RS

counts_raw_sen=read.delim("~/hayflick/data/RNAseq/count_tables/RS_raw_counts.txt", stringsAsFactors = FALSE)
counts_raw_sen=as.matrix(counts_raw_sen)
###
###         combat-seq command and output matrix
#######
adjusted <- ComBat_seq(counts_raw_sen, batch=batch, group=group)

#if out put is list, convert to matrix

### now since we started with counts we need to turn combat adjusted counts back into tpms. so we need lengths from TXimport object from above

gene_length= read.delim("~/hayflick/annotations/RS_gene_length_estimates_by_sample.txt", stringsAsFactors = FALSE)


#####
##    CONVERT TO TPMs
##
tpm=adjusted/gene_length

tpm.mat <- t( t(tpm) * 1e6 / colSums(tpm) )

tpm.mat=as.data.frame(tpm.mat)

##use or load from file (provided)

##LOADTPM tables

### and then convert into log2 fold change for each time point vs. first time point (for that condition)

# FIRST COMBAT RS data

sen_TPMS= read.delim("~/hayflick/data/RNAseq/TPM_tables/RS_batch_tpms.txt", stringsAsFactors = FALSE)
rownames(sen_TPMS)=make.unique(sen_TPMS$gene)
sen_TPMS=sen_TPMS[-1]
sen_TPMS=sen_TPMS+1
sen_TPMS_day1_avg=sen_TPMS %>% dplyr::select(contains("TP1")) %>% rowMeans()
sen_TPMS_l2fc=log2(sen_TPMS/sen_TPMS_day1_avg)


# TERT  data
tert_TPMS= read.delim("~/hayflick/data/RNAseq/TPM_tables/TERT_batch_tpms.txt", stringsAsFactors = FALSE)
rownames(tert_TPMS)=make.unique(tert_TPMS$gene)
tert_TPMS=tert_TPMS[-1]
tert_TPMS=tert_TPMS+1
tert_TPMS_day1_avg=tert_TPMS %>% dplyr::select(contains("TP1")) %>% rowMeans()
tert_TPMS_l2fc=log2(tert_TPMS/tert_TPMS_day1_avg)

# cell density
CD_TPMS= read.delim("~/hayflick/data/RNAseq/TPM_tables/CD_batch_tpms.txt", stringsAsFactors = FALSE)
rownames(CD_TPMS)=make.unique(CD_TPMS$gene)
CD_TPMS=CD_TPMS[-1]
CD_TPMS=CD_TPMS+1
CD_TPMS_day1_avg=CD_TPMS %>% dplyr::select(contains("D1_")) %>% rowMeans()
CD_TPMS_l2fc=log2(CD_TPMS/CD_TPMS_day1_avg)

# radiation
RIS_TPMS= read.delim("~/hayflick/data/RNAseq/TPM_tables/RIS_batch_tpms.txt", stringsAsFactors = FALSE)
rownames(RIS_TPMS)=make.unique(RIS_TPMS$gene)
RIS_TPMS=RIS_TPMS[-1]
RIS_TPMS=RIS_TPMS+1
RIS_TPMS_day1_avg=RIS_TPMS %>% dplyr::select(contains("d0")) %>% rowMeans()
RIS_TPMS_l2fc=log2(RIS_TPMS/RIS_TPMS_day1_avg)


#######merge TPMs
TPM_FC_MERGE=merge(tert_TPMS_l2fc,sen_TPMS_l2fc, by="row.names")

TPM_FC_MERGE=merge(TPM_FC_MERGE,RIS_TPMS_l2fc, by.x="Row.names", by.y="row.names")

TPM_FC_MERGE=merge(TPM_FC_MERGE,CD_TPMS_l2fc, by.x="Row.names", by.y="row.names")

###################### load deseq filter by sig

p=.01

RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = FALSE)
RS_deseq_a=RS_deseq %>% dplyr::filter(padj < p) 
length(RS_deseq_a$gene)

CD_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/CD_deseq.txt", stringsAsFactors = FALSE)
CD_deseq_a=CD_deseq %>% dplyr::filter(padj < p) 
length(CD_deseq_a$gene)

RIS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RIS_deseq.txt", stringsAsFactors = FALSE)
RIS_deseq_a=RIS_deseq %>% dplyr::filter(padj < p) 
length(RIS_deseq_a$gene)

#universe of adjusted p<0.01  gene IDs
sig_ids_a=unique(c(RS_deseq_a$gene,CD_deseq_a$gene,RIS_deseq_a$gene))



###################### merge DEseq output

all_desq=merge(RS_deseq,CD_deseq,by="gene")
all_desq=merge(all_desq, RIS_deseq, by="gene")
rownames(all_desq)=make.unique(all_desq$gene)
all_desq=all_desq[-1]


############################create matrix sig genes


SIG_GENES_TPM_FC_MERGE=TPM_FC_MERGE[TPM_FC_MERGE$Row.names %in% rownames(all_desq),]

SIG_GENES_TPM_FC_MERGE=SIG_GENES_TPM_FC_MERGE[SIG_GENES_TPM_FC_MERGE$Row.names %in% sig_ids_a,]


############################### remmove first TPs used for Fold change calculation from each condition


SIG_GENES_TPM_FC_MERGE = SIG_GENES_TPM_FC_MERGE %>% dplyr::select(!contains("TP1"))

SIG_GENES_TPM_FC_MERGE = SIG_GENES_TPM_FC_MERGE %>% dplyr::select(!contains("d0"))


SIG_GENES_TPM_FC_MERGE = SIG_GENES_TPM_FC_MERGE %>% dplyr::select(!contains("D1_"))

######################### CLUSTER AND HEATMAP


rwb<-colorRampPalette(c("steel blue","white","tomato"))
rwb=rwb(100)

rownames(SIG_GENES_TPM_FC_MERGE)=SIG_GENES_TPM_FC_MERGE$Row.names
SIG_GENES_TPM_FC_MERGE=SIG_GENES_TPM_FC_MERGE[-1]


pheatmap(SIG_GENES_TPM_FC_MERGE,
  
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames = FALSE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         breaks = seq(-3, 3, length.out=101))


############################# 1F GSEA
##
##############################################################

 ##load hallmarks from msigdb

pathways2 <- gmtPathways("~/hayflick/annotations/h.all.v7.0.symbols_.gmt") 

fullGSEA <- vector("list", length=3)
names(fullGSEA) <- c("sen_ranks","con_ranks","IR_ranks")

RIS_deseq_cc=RIS_deseq[complete.cases(RIS_deseq),]
CD_deseq_cc=CD_deseq[complete.cases(CD_deseq),]
RS_deseq_cc=RS_deseq[complete.cases(RS_deseq),]


sen_ranks <- setNames(RS_deseq_cc[,3], RS_deseq_cc[,1])
con_ranks <- setNames(CD_deseq_cc[,3], CD_deseq_cc[,1])
IR_ranks<- setNames(RIS_deseq_cc[,3], RIS_deseq_cc[,1])


rank_list<- vector("list", length=3)
names(rank_list) <- c("sen_ranks","con_ranks","IR_ranks")

rank_list[[1]]=sen_ranks
rank_list[[2]]=con_ranks
rank_list[[3]]=IR_ranks


for (set in 1:length((rank_list))) {
  print(names(rank_list)[set])
  ranks <- rank_list[[set]]
  
  fgseaRes <- fgseaMultilevel(pathways2, ranks ,minSize=10, maxSize=1500,eps = 0)
  fullGSEA[[set]] <- fgseaRes
}

# Create pval and NESobjects objects for pruning

##
pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
#pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
betaGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 6)]) ##col 5 is NES

###this builds the whole matrix, but it's full of NAs
for (l in 2:length(fullGSEA)) {
  pvalGSEA <- merge(pvalGSEA, as.data.frame(fullGSEA[[l]][,c(1, 3)]), by.x=1, by.y=1, all=TRUE)
  #pvalGSEA <- merge(pvalGSEA, as.data.frame(fullGSEA[[l]][,c(1, 3)]), by.x=1, by.y=1, all=TRUE)
  betaGSEA <- merge(betaGSEA, as.data.frame(fullGSEA[[l]][,c(1, 6)]), by.x=1, by.y=1, all=TRUE)
}

colnames(pvalGSEA) <- c("pathway", paste(names(fullGSEA), "pval", sep="_"))
colnames(betaGSEA) <- c("pathway", paste(names(fullGSEA), "NES", sep="_"))


#betaGSEA=betaGSEA[complete.cases(betaGSEA),]

###this is VERY sloppy, but now I pass through again to retrieve the missing values
 for (l in 1:length(fullGSEA)) {
matches <- match(pvalGSEA$pathway, fullGSEA[[l]]$pathway)
 pvalGSEA[,l+1] <- fullGSEA[[l]]$pval[matches]
pvalGSEA[,l+1] <- fullGSEA[[l]]$padj[matches]
  betaGSEA[,l+1] <- fullGSEA[[l]]$NES[matches]
 }

pvalGSEA[,2:ncol(pvalGSEA)] <- -log10(pvalGSEA[,2:ncol(pvalGSEA)]) * sign(betaGSEA[,2:ncol(betaGSEA)])

###convert to numeric matrices with named rows for easy clustering
pvalGSEA2 <- as.matrix(pvalGSEA[,2:ncol(pvalGSEA)])
rownames(pvalGSEA2) <- pvalGSEA[,1]


betaGSEA2 <- as.matrix(betaGSEA[,2:ncol(betaGSEA)])
rownames(betaGSEA2) <- betaGSEA[,1]
#betaGSEA2=betaGSEA2[-1]

pvalGSEA2=pvalGSEA2[complete.cases(pvalGSEA2),]
betaGSEA2=betaGSEA2[rownames(betaGSEA2) %in% rownames(pvalGSEA2),]

################################## PVALUE CUT OFF Figure
pcut <- -log10(.01)


tmp <- melt(apply(pvalGSEA2, 1, function(x){length(which(abs(x)>=pcut))}))
##i verified that the rows are in the same order as the starting matrix, so I'm just going to do direct indexing removal
pvalGSEA2 <- pvalGSEA2[-which(tmp$value==0),]
betaGSEA2 <- betaGSEA2[-which(tmp$value==0),]



pheatmap(pvalGSEA2,
         
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 50,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         breaks = seq(-4, 4, length.out=101))


############################
###################################### FIGURE 1 supp fig 2

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


long_Data = tibble()
files=list.files("~/hayflick/data/RNAseq/TPM_tables/",pattern = "batch",full.names = TRUE)
for (i in seq(1:length(files))) {
  temp=read.delim(files[i],stringsAsFactors = FALSE) 
  temp=cbind(temp$gene,temp[,2:dim(temp)[2]]+1)
  firTP=rowMeans(temp[2:4])
  temp=cbind(temp[,1],temp[2:dim(temp)[2]]/firTP)
  colnames(temp)[1]="gene"      
  colnames(temp)=strReverse(colnames(temp)) ########## kludge fix for condesing reps name inconsistently 
  temp=temp %>% tidyr::pivot_longer(-eneg) ###eneg = gene
  temp_tracking=str_split_fixed(temp$name, "_",2)
  temp$rep=temp_tracking[,2]
  temp_name=str_split_fixed(files[i],"/",9)
  temp_name=temp_name[,9]
  temp_name=str_split_fixed(temp_name,"_",3)
  temp$condition=temp_name[,1]
  
  #temp$rep=factor(temp$rep)
  temp$rep <- factor(temp$rep, levels = unique(temp$rep))
  temp=temp %>%  mutate(numb = as.integer((rep)))
  test= temp %>% filter(eneg=="CD36")
  
  long_Data=rbind(long_Data,temp)
  test2=long_Data %>% filter(eneg=="CD36")
}

dim(TPM_FC_MERGE)


######
gene="NNMT"
gene.use=long_Data %>% filter(eneg==gene) #%>% filter(condition=="RS")
gene.use$numb=factor(gene.use$numb)
gene.use$condition=factor(gene.use$condition, levels = c("TERT","RS","RIS","CD"))
pdfname=paste0("~/hayflick/",gene,"_box.pdf")

paper_cols_builtin_alpha=c("#6EBF91","#F4797E","#FABE79","#70ADD7")
paper_cols=c("#099146","#E8282B","#F39123","#1B76BA")

dodge <- position_dodge(width = 0.9)


######################### BARGRAPH PLOT VERSION

gene=c("SLC2A1","FOXE1","NNMT","CD36","CDKN1A","CDKN2A","TGFB2","GLB1")




for (i in gene){
  gene.use=long_Data %>% filter(eneg==i) %>% filter(condition=="RS"| condition=="TERT") 
  gene.use$numb=factor(gene.use$numb)
  gene.use$condition=factor(gene.use$condition, levels = c("TERT","RS"))
  pdfname=paste0("~/hayflick/",i,"_bar.pdf")
  
  plot=ggplot(gene.use,aes(x=numb , y=value),color=condition) +
    scale_x_discrete() +
    geom_jitter(aes( x = numb), 
                position = position_jitter(width = .05), alpha = 1,size=1) +
    facet_grid(~ condition,space = "free",scales = "free_x") +
    scale_color_manual(values = paper_cols) +
    scale_fill_manual(values = paper_cols) +
    # geom_boxplot(aes(x = numb,color=condition, fill=condition), outlier.colour = NA, position = dodge,alpha=.5) +
    geom_bar(aes(fill=condition,),stat="summary",position = dodge,alpha=.6)+
    theme_bw() + 
    geom_hline(yintercept = 1,linetype="dashed") +
    theme(legend.position = "none")
  
  
  print(pdfname)
  pdf(pdfname,height=6,width=10)
  print(plot)
  dev.off()
  
}






############################ PROTEIN

# data used for  FIGURE 1 supp fig 2 B
pdata=read.delim("~/hayflick/data/proteomics/Proteomics_raw_proportions_RS_and_tert.txt", stringsAsFactors = FALSE)



########################### # FIGURE 1 supp 4


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

###########


#packages <- c("openxlsx", "readxl", "magrittr", "purrr", "ggplot2")
#library(openxlsx)
#library(readxl)
#library(xlsx)


#write.xlsx(fullGSEA, "~/hayflick/ALL_TPs_gsea.xlsx")

#write.table(avg.fc, "~/hayflick/ALL_TPs_input_gsea.txt", sep="\t")


############
pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
#pvalGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 3)]) ##col 3 is padj, col 2 is pval
betaGSEA <- as.data.frame(fullGSEA[[1]][,c(1, 6)]) ##col 5 is NES

###this builds the whole matrix, change NAs to Pval of 1 and NES of zero

###this builds the whole matrix, but it's full of NAs
for (l in 2:length(fullGSEA)) {
  pvalGSEA <- merge(pvalGSEA, as.data.frame(fullGSEA[[l]][,c(1, 3)]), by.x=1, by.y=1, all=TRUE)
  #pvalGSEA <- merge(pvalGSEA, as.data.frame(fullGSEA[[l]][,c(1, 3)]), by.x=1, by.y=1, all=TRUE)
  betaGSEA <- merge(betaGSEA, as.data.frame(fullGSEA[[l]][,c(1, 6)]), by.x=1, by.y=1, all=TRUE)
}

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
       #  gaps_col = c(5,14,19),
         breaks = seq(-4, 4, length.out=101)
)


########################################



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



writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_1/sessionInfo.txt")









