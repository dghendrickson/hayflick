
library("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)
library(colorspace)
library(BiocGenerics)
library(GenomicRanges)
library(Granges)
library(preprocessCore)
library(tidyr)
library(ChIPpeakAnno)
library(regioneR)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)
library(pheatmap)
BiocManager::install("BiocGenerics")
#BiocManager::install("Granges")
#BiocManager::install("GenomicRanges")
#BiocManager::install("org.Hs.eg.db")
#org.Hs.eg.db
#BiocManager::install("ChIPpeakAnno")

#BiocManager::install("preprocessCore")

library(EnsDb.Hsapiens.v75)
library(preprocessCore)

library(AnnotationHub)
#BiocManager::install("AnnotationHub")
#BiocManager::install("org.Hs.egALIAS2EG")

#BiocManager::install("preprocessCore")


#######LOAD PEAKS PEAK ATLAS, PEAK SIGS, and Peak Quants are in same order. 
####I use peakAtlas@elementMetadata$region as the unique identifier to subset and merge when necessary

peakAtlas=readRDS("~/hayflick/data/ATACseq/peak_atlas.rds", refhook = NULL)

peak_sigs= readRDS("~/hayflick/data/ATACseq/limmaVoom_DE_res.rds", refhook = NULL)

peakQuants= readRDS("~/hayflick/data/ATACseq/quantile_norm_peak_counts.rds", refhook = NULL)
peakQuants=as.data.frame(peakQuants)
rownames(peakQuants)=peakAtlas@elementMetadata$region



chrom_quants=readRDS("~/hayflick/data/ATACseq/chromHMM_counts_quantile_norm.rds",refhook = NULL)
chrom_quants=as.data.frame(chrom_quants)


########remove unused columns--Bad data--mislabeled samples

chrom_quants=chrom_quants[-15]
chrom_quants=chrom_quants[-17]
chrom_quants=chrom_quants[-34]
chrom_quants=chrom_quants[-36]

################## Quantile norm for visulaization

chrom_quants_mat=as.matrix(chrom_quants[1:35])
chrom_quants_mat=normalize.quantiles(chrom_quants_mat)

chrom_quants_qn=as.data.frame(chrom_quants_mat)

colnames(chrom_quants_qn)=colnames(chrom_quants)

chrom_quants=chrom_quants_qn

chrom_quants$genes=rownames(chrom_quants)


###load imr90 chromHMM for hg38

imr90chrmm=import("/home/dgh/hayflick/annotations/E017_25_imputed12marks_hg38lift_dense.bed", format = "BED")

chrom_quants$state=imr90chrmm$name
#gdf1 = gather(chrom_quants, "group", "Expression",-genes)

gdf1=pivot_longer(chrom_quants,!c("genes","state"),names_to = "sample", values_to = "reads")

temp=str_split_fixed(gdf1$sample,"_",3)
gdf1$tp=temp[,2]
gdf1$type=temp[,1]


temp=str_split_fixed(gdf1$state,"_",3)
gdf1$cat.num=temp[,1]


### relabel to 4 categories
x=case_when(gdf1$cat.num == 1 ~"promoter",
            gdf1$cat.num == 2 ~"promoter",
            gdf1$cat.num == 3 ~"promoter",
            gdf1$cat.num == 4 ~"promoter",
            gdf1$cat.num == 5 ~"transcribed",
            gdf1$cat.num == 6 ~"transcribed",
            gdf1$cat.num == 7 ~"transcribed",
            gdf1$cat.num == 8 ~"transcribed",
            gdf1$cat.num == 9 ~"enhancer",
            gdf1$cat.num == 10 ~"enhancer",
            gdf1$cat.num == 11 ~"enhancer",
            gdf1$cat.num == 12 ~"enhancer",
            gdf1$cat.num == 13 ~"enhancer",
            gdf1$cat.num == 14~"enhancer",
            gdf1$cat.num == 15 ~"enhancer",
            gdf1$cat.num == 16 ~"enhancer",
            gdf1$cat.num == 17 ~"enhancer",
            gdf1$cat.num == 18 ~"enhancer",
            gdf1$cat.num == 19 ~"misc",
            gdf1$cat.num == 20 ~"misc",
            gdf1$cat.num == 21 ~"misc",
            gdf1$cat.num == 22 ~"promoter",
            gdf1$cat.num == 23 ~"promoter",
            gdf1$cat.num == 24 ~"misc",
            gdf1$cat.num == 25 ~"misc")


gdf1$big_state=x


chrom_quants_sum_ALL=gdf1 %>% group_by(type,big_state) %>% dplyr::summarize(sum_reads = sum(reads))  %>% group_by(type)  %>%  mutate(tp_tot=sum(sum_reads)) %>% 
  group_by(type,big_state) %>% summarise(percent=sum_reads/tp_tot*100) 



####Figure 5A

pdf("/home/dgh/test/tert.pdl_chrom.state.quant.pdf",height=6,width=4)

ggplot(chrom_quants_sum_ALL,aes(x=type, y= percent,fill=big_state )) +
  geom_bar(position = "fill",stat = "identity") + scale_fill_discrete_qualitative(palette="Dark3") +
  # theme(legend.position = "none") +
  theme(axis.line = element_line(size=.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 15)) +
  theme(aspect.ratio = 4/1)

dev.off()


##################################################################### IDIVIDUAL TIME POINTS; Figure 5A right panel

chrom_quants_by_sample_sum=gdf1 %>% group_by(sample,big_state)   %>% dplyr::summarize(sum_reads = sum(reads))   %>% group_by(sample) %>%  mutate(tp_tot=sum(sum_reads)) %>% 
  group_by(sample,big_state) %>% mutate(percent=sum_reads/tp_tot*100) 


temp=str_split_fixed(chrom_quants_by_sample_sum$sample,"_",3)
chrom_quants_by_sample_sum$type=temp[,1]
chrom_quants_by_sample_sum$tp=temp[,2]
chrom_quants_by_sample_sum$re=temp[,3]
#chrom_quants_sum$state[chrom_quants_sum$state == '20_ZNF/Rpts'] <- '20_ZNF_Rpts'


state=c("misc")
state=c("promoter")

state.use=chrom_quants_by_sample_sum %>% filter(big_state==state) #%>% filter(condition=="RS")

 
  state.use$tp=factor(state.use$tp)
  state.use$type=factor( state.use$type, levels = c("hTERT","PDL"))
  pdfname=paste0("~/hayflick/",state,"_bar.pdf")
  
  plot=ggplot(state.use,aes(x=tp , y= percent),color=type) +
    scale_x_discrete() +
    geom_jitter(aes(x = tp), 
                position = position_jitter(width = .05), alpha = 1,size=1) +
    facet_grid(~ type,space = "free",scales = "free_x") +
    scale_color_manual(values = paper_cols) +
    scale_fill_manual(values = paper_cols) +
    # geom_boxplot(aes(x = numb,color=condition, fill=condition), outlier.colour = NA, position = dodge,alpha=.5) +
    geom_bar(aes(fill=type,),stat="summary",position = dodge,alpha=.6)+
    theme_bw() + 
    geom_hline(yintercept = 1,linetype="dashed") +
    theme(legend.position = "none")
  
  
  print(pdfname)
  pdf(pdfname,height=6,width=7)
  print(plot)
  dev.off()
  

################################# PEAKS in STATES HM ALL 25 FOR Figure 5B

#### GET SIG PEAKS FROM LIMMA VOOM & ADD TO PEAK ATALS 
  peak_sigs_tp7=peak_sigs[[5]]
  peak_sigs_tp6=peak_sigs[[4]]
  peak_sigs_tp5=peak_sigs[[3]]
  peak_sigs_tp4=peak_sigs[[2]]
  peak_sigs_tp2=peak_sigs[[1]]
  
  #peakAtlas@elementMetadata$l2_FC=dds_PDL_results$log2FoldChange
  
  peakAtlas@elementMetadata$tp7_p=peak_sigs_tp7$adj.P.Val
  peakAtlas@elementMetadata$tp6_p=peak_sigs_tp6$adj.P.Val
  peakAtlas@elementMetadata$tp5_p=peak_sigs_tp5$adj.P.Val
  peakAtlas@elementMetadata$tp4_p=peak_sigs_tp4$adj.P.Val
  peakAtlas@elementMetadata$tp2_p=peak_sigs_tp2$adj.P.Val
  
  peakAtlas@elementMetadata$tp7_l2=peak_sigs_tp7$logFC
  
  peakAtlas@elementMetadata$atac_peak=rownames(peakAtlas)
  
  forSup=peakAtlas@elementMetadata
  forSup=as.data.frame(forSup)
  
  chrom_ids=data_fram_overla_df_single %>% select(c("region","name"))
  
  forSup_m=merge(forSup,chrom_ids, by.x="region",by.y="region",all=TRUE)
  
  forSup_mm=merge(forSup_m, peakQuants,by.x="region" ,by.y="row.names")
  
  write.table(forSup_m,"~/hayflick/data/ATACseq/peaks_metaData.txt", sep="\t")
  
  write.table(peakQuants,"~/hayflick/data/ATACseq/peak_Quants_qn.txt", sep="\t")
  
  write.table(forSup_mm,"~/hayflick/data/ATACseq/peak_Quants_AND_meta_qn.txt", sep="\t")
  
  
  #### compile Universe of sig peaks 
  p=0.001
  sig=which((peakAtlas$tp7_p < p))
  peakAtlas_sig=peakAtlas[sig]
  tp7sig=peakAtlas_sig$region
  
  
  sig=which((peakAtlas$tp6_p < p))
  peakAtlas_sig=peakAtlas[sig]
  tp6sig=peakAtlas_sig$region
  
  
  sig=which((peakAtlas$tp5_p < p))
  peakAtlas_sig=peakAtlas[sig]
  tp5sig=peakAtlas_sig$region
  
  
  sig=which((peakAtlas$tp4_p < p))
  peakAtlas_sig=peakAtlas[sig]
  tp4sig=peakAtlas_sig$region
  
  
  sig=which((peakAtlas$tp2_p < p))
  peakAtlas_sig=peakAtlas[sig]
  tp2sig=peakAtlas_sig$region
  
  sig_list=c(tp7sig,tp6sig,tp5sig,tp4sig,tp2sig)
  
  sig_list=unique(sig_list)

###################################


  peakAtlas_sig=peakAtlas[peakAtlas$region %in% sig_list]
  
  sig=which(peakAtlas_sig$tp7_l2 > 0)
  peakAtlas_sig_up=peakAtlas_sig[sig]
  
  sig=which(peakAtlas_sig$tp7_l2 < 0)
  peakAtlas_sig_down=peakAtlas_sig[sig]


  ################OVERLAP SIG PEAKS with Chrom states
  library(ChIPpeakAnno)
  
  ol_trd_dmin <- findOverlapsOfPeaks(peakAtlas_sig, imr90chrmm)
  
  #ol_trd_dmin_newSig=ol_trd_dmin
  
  data_fram_overla_df=as.data.frame(ol_trd_dmin$overlappingPeaks[[1]]) ### overlap with ChipPeakANNO, default settings

  data_fram_overla_df_sub=data_fram_overla_df %>% 
    dplyr::select(c("region","annot","gene_name","tp7_p","tp7_l2","name","overlapFeature","shortestDistance")) ## sub on needed columns
  
  
  data_fram_overla_df_inside_include=data_fram_overla_df_sub %>% 
    dplyr::filter(overlapFeature == "inside" | overlapFeature == "includeFeature") ##get unambiguos overalps i.e. peak falls inside ChromHmm state or Chrom state falls inside Peak
  
  lv_m=!data_fram_overla_df_inside_include$region %in% 
    data_fram_overla_df_inside_include$region[duplicated(data_fram_overla_df_inside_include$region)] ## remove peaks that include more than 1 Chrom state
  
  data_fram_overla_df_single=data_fram_overla_df_inside_include[lv_m,] ### whats left!

  peakQuants_single_sig=peakQuants[rownames(peakQuants) %in% data_fram_overla_df_single$region,] #### get Quantile norm counts for unambigouus overlaps of sig peaks
  
  
  ################ FOLD CHANGE CONVERSION
  
  
  peakQuants_single_sig=as.data.frame(peakQuants_single_sig)
  
  
  peakQuants_single_sig_tert=peakQuants_single_sig %>% dplyr::select(contains("hTERT"))
  
  peakQuants_single_sig_tert_d1rm=rowMeans(peakQuants_single_sig_tert[,1:3])
  
  peakQuants_single_sig_tert_l2fc=(peakQuants_single_sig_tert-peakQuants_single_sig_tert_d1rm)
  
  peakQuants_single_sig_PDL=peakQuants_single_sig %>% dplyr::select(contains("PDL"))
  peakQuants_single_sig_PDL_d1rm=rowMeans(peakQuants_single_sig_PDL[,1:3])
  
  peakQuants_single_sig_PDL_l2fc=(peakQuants_single_sig_PDL-peakQuants_single_sig_PDL_d1rm)
  
  peakQuants_single_sig_l2fc=cbind(peakQuants_single_sig_tert_l2fc,peakQuants_single_sig_PDL_l2fc) ########### FOLD CHANGE OBJECT FOR HM
  
  ############# CLEAN IT UP
  
  
  peakQuants_single_sig_l2fc=peakQuants_single_sig_l2fc[complete.cases(peakQuants_single_sig_l2fc),] ## no NAs
  
  peakQuants_single_sig_l2fc=peakQuants_single_sig_l2fc[as.logical(rowSums(peakQuants_single_sig_l2fc != 32)), ] ### remove peaks with replicates that do not contain reads
  
  
  peakQuants_single_sig_l2fc <- peakQuants_single_sig_l2fc[is.finite(rowSums(peakQuants_single_sig_l2fc)),] ## no INF
  
  
  
  
  ##### GET PEAK ANNOTATIONS e.g. chrom state or genomics region based on gene annotation (peakAtlas@elementMetadata$annot)
  
   my_row_annot=as.data.frame(rownames(peakQuants_single_sig_l2fc))
  
  rownames(my_row_annot)=my_row_annot$`rownames(peakQuants_single_sig_l2fc)`
  my_row_annot$anno_state=data_fram_overla_df_single$name
  my_row_annot$anno_type=data_fram_overla_df_single$annot
  my_row_annot=my_row_annot[-1]
  
  
  peaksQ_w_annot=cbind(peakQuants_single_sig_l2fc, my_row_annot$anno_state)
  colnames(peaksQ_w_annot)[33]="state"
  

  
  newcenters=peaksQ_w_annot %>% group_by(state) %>% summarise_all(funs(median))
  newcenters=as.data.frame(newcenters)
  rownames(newcenters)=newcenters$state
  peaksQnewcenters=newcenters[,2:33]
  
  
  
  ###Figure 5B
  
  pheatmap(peaksQnewcenters,
           fontsize_row = 10,
           treeheight_row = 0,
           show_rownames = TRUE,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_colnames = TRUE,
           color = diverging_hcl(100, "Tropic",alpha = .7),
           gaps_col = c(16),
           breaks = seq(-.75, .75, length.out=101))
     

  
  ######################NADS LADS FIGURE 5D

  # TEAD1 and AP1

  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(rtracklayer)
  library(VennDiagram)
  library(seqLogo)

  atlas = peakAtlas
  # load atlas
  NAD <- import.bed("~/hayflick/annotations/hg38_im90_nads.bed")
  LAD <- import.bed("~/hayflick/annotations//GSE49341_LMNB1_Gro.bed.gz")
  
  atlas$NAD <- 0
  atlas$NAD[unique(queryHits(findOverlaps(atlas, NAD)))] <- 1
  atlas$LAD <- 0
  atlas$LAD[unique(queryHits(findOverlaps(atlas, LAD)))] <- 1
  
  # load motifs
  rna  = read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt",stringsAsFactors = FALSE) # load DESEQ results
  # add JUN_hocomoco motif

  ##motif RDT files generated with "~/hayflick/data_processing/ATAC/motif_process/motif_analysis.R"
  
  load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/JUN_motif_HOCOMOCO/motif_hits.rdt")
  colnames(m) <- "JUN_hm"
  m_jun <- m
  # add FOXE1 motif
  load("~/hayflick/data_processing/ATAC/motif_process/motif_databases.12.18/FOXE1_motif/motif_hits.rdt")
  m_FOXE1 <- m
  load("~/hayflick/data_processing/ATAC/motif_process/fimo/motif_hits.rdt")
  m <- m[, colnames(m) %in% rownames(rna[[1]])] # filter for only expressed TFs
  m <- cbind(m, m_jun, m_FOXE1)
  
  # senescent log2fc
  res <- peak_sigs
  
  ########################################################
  # confirm NAD/LAD associated with increased senescence #
  ########################################################
  # NAD
  toplot <- data.frame(log2FC=res$TP7$logFC,
                       NAD=atlas$NAD==1)
  x <- toplot[toplot$NAD==1, "log2FC"]
  y <- toplot[toplot$NAD==0, "log2FC"]
  p <- wilcox.test(x, y, alternative="greater")$p.value
  gplot1 <- ggplot(toplot, aes(x=NAD,y=log2FC))+geom_violin(aes(fill=NAD)) + 
    theme_classic() + ggtitle(sprintf('p < %.2e', p)) + ylab("senescence log2FC")
  
  
  # LAD
  toplot <- data.frame(log2FC=res$TP7$logFC,
                       LAD=atlas$LAD==1)
  x <- toplot[toplot$LAD==1, "log2FC"]
  y <- toplot[toplot$LAD==0, "log2FC"]
  p <- wilcox.test(x, y, alternative="greater")$p.value
  gplot2 <- ggplot(toplot, aes(x=LAD,y=log2FC))+geom_violin(aes(fill=LAD)) + 
    theme_classic() + ggtitle(sprintf('p < %.2e', p))  + ylab("senescence log2FC")
  
  mean(gplot1$data[gplot1$data$NAD,"log2FC"]) - mean(gplot1$data[!gplot1$data$NAD,"log2FC"])
  mean(gplot2$data[gplot2$data$LAD,"log2FC"]) - mean(gplot2$data[!gplot2$data$LAD,"log2FC"])
  
  pdf("analysis/12_14_20_NAD_LAD/NAD_LAD_associated_with_increase_senescence.pdf", 5, 6)
  grid.arrange(gplot1, gplot2, ncol=2)
  dev.off()
  
  ###################################
  # are TEAD1 enriched in NAD/LAD ? #
  ###################################
  # NAD
  x <- sum(atlas$NAD==1 & m[,"TEAD1"]==1) # white ball drawn
  a <- sum(m[,"TEAD1"]==1) # total white balls
  b <- sum(m[,"TEAD1"]==0) # total black balls
  k <- sum(atlas$NAD==1) # total drawn
  p <- phyper(x, a, b, k, lower.tail = F)
  toplot <- table(ifelse(m[,"TEAD1"]==1, "TEAD1+", "TEAD1-"),
                  ifelse(atlas$NAD==1, "NAD+", "NAD-"))
  toplot_norm <- t(t(toplot)/colSums(toplot))
  
  pdf("analysis/12_14_20_NAD_LAD/TEAD1_sites_are_enriched_in_NAD_domains.pdf")
  par(mfrow=c(1,2))
  barplot(toplot, col=c("darkblue","red"), main=sprintf("p-value < %.2e", p))
  legend("topright", fill=c("darkblue","red"), legend=c("TEAD1-","TEAD1+"))
  barplot(toplot_norm, col=c("darkblue","red"), main="normalized barplot")
  dev.off()
  
  # LAD
  x <- sum(atlas$LAD==1 & m[,"TEAD1"]==1) # white ball drawn
  a <- sum(m[,"TEAD1"]==1) # total white balls
  b <- sum(m[,"TEAD1"]==0) # total black balls
  k <- sum(atlas$LAD==1) # total drawn
  p <- phyper(x, a, b, k, lower.tail = F)
  toplot <- table(ifelse(m[,"TEAD1"]==1, "TEAD1+", "TEAD1-"),
                  ifelse(atlas$LAD==1, "LAD+", "LAD-"))
  toplot_norm <- t(t(toplot)/colSums(toplot))
  
  pdf("analysis/12_14_20_NAD_LAD/TEAD1_sites_are_enriched_in_LAD_domains.pdf")
  par(mfrow=c(1,2))
  barplot(toplot, col=c("darkblue","red"), main=sprintf("p-value < %.2e", p))
  legend("topright", fill=c("darkblue","red"), legend=c("TEAD1-","TEAD1+"))
  barplot(toplot_norm, col=c("darkblue","red"), main="normalized barplot")
  dev.off()
  
  ###############
  # interaction #
  ###############
  # NAD: significant
  toplot <- data.frame(log2FC=res$TP7$logFC,
                       NAD=atlas$NAD,
                       TEAD1=m[,"TEAD1"])
  fit <- lm(log2FC ~ NAD*TEAD1, data = toplot)
  summary(fit)
  
  # LAD: not signficant
  toplot <- data.frame(log2FC=res$TP7$logFC,
                       LAD=atlas$LAD,
                       TEAD1=m[,"TEAD1"])
  fit <- lm(log2FC ~ LAD*TEAD1, data = toplot)
  summary(fit)
  
  
  
  ############################################################################## FIGURE 5C
  ###################################
  imr90_nads=import("/home/dgh/test/hg38_im90_nads.bed", format = "BED") #load NADs
  
  
  imr90_lads=import("/home/dgh/test/imr90_hg38_lads.bed", format = "BED") #load LADs
  
  
  
library(TxDb.Hsapiens.UCSC.hg38.knownGene)   ### load UCSC hg38 annotations ENTREZ IDs
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(org.Hs.eg.db)
hs <- org.Hs.eg.db


all.genes <- genes(txdb)        ### gRANGES ALL GENES ENTREZ IDs



human.genome <- getGenomeAndMask("hg38", mask=NA)$genome
human.auto <- filterChromosomes(human.genome, chr.type="autosomal",organism="hg")
human.canon <- filterChromosomes(human.genome, chr.type="canonical",organism="hg")

  
  ######## GET INDUCED GENES
RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = FALSE)
CD_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/CD_deseq.txt", stringsAsFactors = FALSE)
RIS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RIS_deseq.txt", stringsAsFactors = FALSE)


  deseq=CD_deseq  ### set induced genes to test
  deseq=RS_deseq
  deseq=RIS_deseq
  
  my.symbols=deseq$gene
  
  entrez_ids=select(hs,                                 #retrieve gene symbol by ENTREZ ids
                    keys = my.symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  entrez_ids=entrez_ids[complete.cases(entrez_ids),]
  
  
  deseq_Entrez=merge(deseq, entrez_ids,by.x="gene",by.y="SYMBOL") ##merge with DESEQ results
  deseq_Entrez=deseq_Entrez[complete.cases(deseq_Entrez),]
  deseq_Entrez_F=deseq_Entrez[order(deseq_Entrez$log2FoldChange,decreasing = TRUE),] %>% filter(padj < .01)  %>% filter(log2FoldChange > .5) ### filter sig and highest FC

  
  ######################################## THESE ARE GENES TO TEST OVERAL
  top1000=deseq_Entrez_F[,c("ENTREZID")]
  #top1000=deseq_Entrez_F[1:1000,c("ENTREZID")]
  temp_gr=subset(all.genes,  gene_id %in% top1000)               ####subset all genes GRanges by genes indiced (geth their genomic intervals)
  #############

  
  
  
  ############## Permutation test using the above "permuteRegionsMetadata" function
 
  
  genic_test <- permTest(A=temp_gr, B=imr90_lads,
                         ntimes=1000,
                         genome=human.canon,
                         per.chromosome=TRUE,
                         randomize.function=permuteRegionsMetadata,
                         evaluate.function=numOverlaps, verbose = TRUE)
  
  plot(genic_test)
  
  
  sen1000_LAD=genic_test
 con1000_LAD=genic_test
RIS1000_LAD=genic_test
 
  sen1000_gt=genic_test
  
  sen1000_gt$numOverlaps$permuted= sen1000_gt$numOverlaps$permuted/length(temp_gr)
  sen1000_gt$numOverlaps$observed =sen1000_gt$numOverlaps$observed/length(temp_gr)
  
  
    con1000_gt=genic_test
    con1000_gt$numOverlaps$permuted= con1000_gt$numOverlaps$permuted/length(temp_gr)
    con1000_gt$numOverlaps$observed =con1000_gt$numOverlaps$observed/length(temp_gr)
    
    
    
    IR1000_gt=genic_test
    genic_test$numOverlaps$permuted= genic_test$numOverlaps$permuted/length(temp_gr)
    genic_test$numOverlaps$observed =genic_test$numOverlaps$observed/length(temp_gr)
    
    
    
    sen1000_LAD$numOverlaps$permuted= sen1000_LAD$numOverlaps$permuted/length(temp_gr)
    sen1000_LAD$numOverlaps$observed =sen1000_LAD$numOverlaps$observed/length(temp_gr)
    
    

    con1000_LAD$numOverlaps$permuted= con1000_LAD$numOverlaps$permuted/length(temp_gr)
    con1000_LAD$numOverlaps$observed =con1000_LAD$numOverlaps$observed/length(temp_gr)
    
    

    RIS1000_LAD$numOverlaps$permuted= RIS1000_LAD$numOverlaps$permuted/length(temp_gr)
    RIS1000_LAD$numOverlaps$observed =RIS1000_LAD$numOverlaps$observed/length(temp_gr)
    
    
    
    
  #########################
  
  pdf("~/hayflick/GE_sen_nads.pdf", height=5,width=5)
  plot(sen1000_gt,xlim = c(.12,.19), ylim= c(0,70))
  dev.off()
  
  pdf("~/hayflick/GE_con_nads.pdf", height=5,width=5)
plot(con1000_gt,xlim = c(.12,.19), ylim= c(0,70))
  dev.off()
  
  pdf("~/hayflick/GE_IR_nads.pdf", height=5,width=5)
plot(IR1000_gt,xlim = c(.12,.19), ylim= c(0,70))
  dev.off()
  
  ############
  pdf("~/hayflick/GE_sen_lads.pdf", height=5,width=5)
  plot(sen1000_LAD,xlim = c(.35,.50), ylim= c(0,30))
  dev.off()
  
  pdf("~/hayflick/GE_con_lads.pdf", height=5,width=5)
  plot(con1000_LAD,xlim = c(.35,.50), ylim= c(0,30))
  dev.off()
  
  pdf("~/hayflick/GE_IR_lads.pdf", height=5,width=5)
  plot(RIS1000_LAD,xlim = c(.35,.50), ylim= c(0,30))
  dev.off()
  
  
  ############################################################## FIgure 5 supp 5
  
 
  im90_sig=peakAtlas_sig[(elementMetadata(peakAtlas_sig)[, "region"] %in% rownames(peaksQ_w_annot))]
  im90_sig$state=peaksQ_w_annot$state

  
  im90_sig$NAD <- 0
  im90_sig$NAD[unique(queryHits(findOverlaps(im90_sig, NAD)))] <- 1
  
 # im90_sig_up$LAD <- 0
 # im90_sig_up$LAD[unique(queryHits(findOverlaps(im90_sig_up, LAD)))] <- 1
  
  
  peaksQ_w_annot$NAD=im90_sig$NAD
  
 # peaksQ_w_annot$LAD=im90_sig_up$LAD
  
  newcenters=peaksQ_w_annot %>% group_by(state,NAD) %>% summarise_all(list(median))
# newcenters=as.data.frame(newcenters)

  newcenters_nad = newcenters %>% filter(NAD == 1)
  newcenters_nad=newcenters_nad %>% column_to_rownames("state")
  
  newcenters_rest = newcenters %>% filter(NAD == 0)
  newcenters_rest=newcenters_rest %>% column_to_rownames("state")

  newcenters_wide=merge(newcenters_rest,newcenters_nad, by="row.names")

  rownames(newcenters_wide)=newcenters_wide$Row.names

  newcenters_wide=newcenters_wide[-1]
  newcenters_wide=newcenters_wide %>% select(!contains("NAD"))
  
  
  
  
  pdf("~/hayflick/NADS_peaks_states.pdf", height=6,width=8)
  
  pheatmap(newcenters_wide,
           fontsize_row = 10,
           treeheight_row = 0,
           show_rownames = TRUE,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           clustering_distance_rows = "correlation",
           show_colnames = TRUE,
           color = diverging_hcl(100, "Tropic",alpha = .7),
           gaps_col = c(32),
           breaks = seq(-.75, .75, length.out=101)
           
  )
  dev.off()
  #######################
  ########################################################### ### FUNCTION TO PERMUTE SETS OF GENES OF simialr expression profile to top 1000
  permuteRegionsMetadata <- function(A, ...) {
    #   A <- toGRanges(A)
    
    
    res_3_cc=deseq                                   #same as above ID gene symbol to ENTREZ
    my.symbols=res_3_cc$gene
    entrez_ids=select(hs, 
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
    
    entrez_ids=entrez_ids[complete.cases(entrez_ids),]
    
    
    res_3_cc=merge(res_3_cc, entrez_ids,by.x="gene",by.y="SYMBOL")
    
    res_3_cc=res_3_cc[complete.cases(res_3_cc),]
    
    
    res_3_cc$cut.50decile <-Hmisc::cut2(log2(res_3_cc$baseMean),g= 50)  # Assign every gene a percentile bin from DESEQ base mean qunatification of normalized counts 
    res_3_cc$cut.50decile=factor(res_3_cc$cut.50decile)
    
    
    x=res_3_cc[res_3_cc$ENTREZID %in% temp_gr$gene_id,]  %>%   # GET the freq distribution of  expression for the top1000 genes (using the temp_gr granges) object
      group_by(cut.50decile) %>%
      summarise (n = n()) %>%
      mutate(freq = n / sum(n))
    prop.v_eq=x$freq
    names(prop.v_eq)=x$cut.50decile
    
    
    if (length(prop.v_eq) < 50) {
      missing=50-length(prop.v_eq)
      num=50-missing+1
      prop.v_eq[num:50]=0
      names(prop.v_eq)[num:50]=setdiff(levels(x$cut.50decile), x$cut.50decile)                       ######replace empty (missing) bins with 0
    }
    
    test <- res_3_cc %>%
      group_by(cut.50decile) %>%  
      nest()                                                 #### nest all genes by expression bin for sampling based on freq distribution "prop.v_eq"
    
    test$n=round(prop.v_eq*length(temp_gr))  ## set number of genes to sample
    
    test=test %>% mutate(samp = map2(data, n, sample_n)) %>%  #sample and un-nest
      dplyr::select(cut.50decile, samp) %>%
      unnest()
    
    
    all.genes_8_match=subset(all.genes, gene_id %in% test$ENTREZID)  ### return Granages object subset from all genes to test for NAD/LAD overalps
    return(all.genes_8_match)
  }
  
  

###################FIGURE 5 supp fig 4 brower figure
  
  library("karyoploteR")
  
  
  
  data_fram_overla_df_single_up=data_fram_overla_df_single %>% dplyr::filter(tp7_l2  >0)
  data_fram_overla_df_single_up_Q=data_fram_overla_df_single_up %>% dplyr::filter(name=="25_Quies")  
  peakAtlas_Q_up=subset(peakAtlas, region %in%  data_fram_overla_df_single_up_Q$region)  
  
  
  
  peakAtlas_samp=sample(peakAtlas, 15000, replace = FALSE)
  
  kp <- plotKaryotype(plot.type=2, chromosomes = "chr22","hg38")
  
  kpPlotRegions(kp, imr90_nads, col="#AACCFF",r0=.05,r1=.25)
  kpPlotRegions(kp, imr90_lads, col="#ffe1aa",r0=.3,r1=.5)
  
  
  kpPlotRegions(kp, peakAtlas_Q_up, col="#FFAACC",r0=1.05,r1=1.25,avoid.overlapping = FALSE)
  #kpPlotRegions(kp, peakAtlas_ID_up, col="#d3aaff",r0=1.3,r1=1.5,avoid.overlapping = FALSE)
  
  
  kpPlotDensity(kp, peakAtlas_samp, window.size = .2e6, col="#aab8f0",r0=.55,r1=.75)
  
  kpPlotDensity(kp, peakAtlas_Q_up, window.size = .2e6, col="#FFAACC",r0=1.05,r1=1.25)
  
  
  ################################################################
  
  
  writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_5/sessionInfo.txt")
  
  
  
  
  
  
  


