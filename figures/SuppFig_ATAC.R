BiocManager::install("karyoploteR")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library("karyoploteR")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all.genes <- genes(txdb)
head(all.genes)


#################
#################
#
frip_score=read.csv("~/hayflick/data/ATACseq/sample_info.txt", stringsAsFactors = FALSE)

frip_score=frip_score[-15,]
frip_score=frip_score[-35,]

pdf("~/hayflick/FRIP_score.pdf", height=5,width=6)

 ggplot(frip_score,aes(x=timepoint , y= reads_in_peaks),color=treatment) +
  scale_x_discrete() +
  geom_jitter(aes(x = timepoint), 
              position = position_jitter(width = .05), alpha = 1,size=1) +
  facet_grid(~ treatment,space = "free",scales = "free_x") +
  scale_color_manual(values = paper_cols) +
  scale_fill_manual(values = paper_cols) +
  # geom_boxplot(aes(x = numb,color=condition, fill=condition), outlier.colour = NA, position = dodge,alpha=.5) +
  geom_bar(aes(fill=treatment,),stat="summary",position = dodge,alpha=.6)+
  theme_bw() + 
  theme(legend.position = "none")
 
 dev.off()

######################################karyoplot
 
 
 data_fram_overla_df_single_up= data_fram_overla_df_single %>% filter( tp7_l2 > 0 ) ######### sig (padj < .001) peaks overlapping umabigiusly with chrom states from fig5 filtered for induced
 
 data_fram_overla_df_single_up_Q=data_fram_overla_df_single_up %>% dplyr::filter(name=="25_Quies") #### now filtered for the Quiescent state (putative NADS)

 peakAtlas_Q_up=subset(peakAtlas, region %in%  data_fram_overla_df_single_up_Q$region) ##### Granges subset based on these filters for karyoplot
 
 peakAtlas_sig_samp=sample(peakAtlas_sig, 15000, replace = FALSE) ######## sample sig peak atlas--too big
 
 ################## PLOT
 kp <- plotKaryotype(plot.type=2, chromosomes = "chr22","hg38") #chr. 22
 
 kpPlotRegions(kp, imr90_nads, col="#AACCFF",r0=.05,r1=.25) #### NADS
 kpPlotRegions(kp, imr90_lads, col="#ffe1aa",r0=.3,r1=.5)   #### LADs

 kpPlotDensity(kp, peakAtlas_Q_up, window.size = .2e6, col="#FFAACC",r0=.8,r1=1.2)
 kpPlotRegions(kp, peakAtlas_Q_up, col="#FFAACC",r0=.6,r1=.8,avoid.overlapping = FALSE)
  
 kpPlotDensity(kp, peakAtlas_sig_samp, window.size = .2e6, col="#d3aaff",r0=1.05,r1=1.3)
 
 #########################################################
 
 
 imr90chrmm_gr=imr90chrmm
 chrm_names=(unique(imr90chrmm_gr@elementMetadata$name))
 
 chrom_state_v_Nads=list()
 
 for (i in seq(1:length(unique(imr90chrmm_gr@elementMetadata$name)))) {
   
   #for (i in seq(1:5)) {
   
   temp_range=subset(imr90chrmm, name== (unique(imr90chrmm_gr@elementMetadata$name)[i]))
   
   genic_test <- permTest(A=temp_range, B=imr90_nads,ntimes=100,genome=human.canon,
                          randomize.function=resampleRegions,per.chromosome=FALSE, universe=imr90chrmm_gr,
                          evaluate.function=numOverlaps, verbose = TRUE)
   
   x=summary(genic_test)
   chrom_state_v_Nads[[i]]=x
 }
 
 
 names(chrom_state_v_Nads)=chrm_names
 
 chrom_state_v_Nads_df=purrr::reduce((chrom_state_v_Nads), dplyr::bind_rows)
 rownames(chrom_state_v_Nads_df)=chrm_names
 
 
 
 chrom_state_v_Lads=list()
 

 for (i in seq(1:length(unique(imr90chrmm_gr@elementMetadata$name)))) {
   
   
   
   temp_range=subset(imr90chrmm, name== unique(imr90chrmm_gr@elementMetadata$name)[i])
   
   genic_test <- permTest(A=temp_range, B=LAD,ntimes=100,genome=human.canon,
                          randomize.function=resampleRegions,per.chromosome=FALSE, universe=imr90chrmm_gr,
                          evaluate.function=numOverlaps, verbose = TRUE)
   
   x=summary(genic_test)
   chrom_state_v_Lads[[i]]=x
 }
 
 
 names(chrom_state_v_Lads)=chrm_names
 
 chrom_state_v_Lads_df=purrr::reduce((chrom_state_v_Lads), dplyr::bind_rows)
 rownames(chrom_state_v_Lads_df)=chrm_names
 
 
 chrom_state.N.L=cbind(chrom_state_v_Nads_df$zscore,chrom_state_v_Lads_df$zscore)
 
 rownames(chrom_state.N.L)=chrm_names
 
 
 rwb<-colorRampPalette(c("steel blue","white","tomato"))
 rwb=rwb(100)
 
 colnames(chrom_state.N.L)=c("NADs","LADs")
 log2(-1)
 
 
 
 
 
 breaks = seq(-100 ,100, length.out=101)
 
 breaks2 = breaks
 
 breaks2[length(breaks)] <- max(max(chrom_state.N.L),max(breaks))
 breaks2[1] <- min(min(chrom_state.N.L),min(breaks))
 
 
 
 ###########################
 pheatmap(chrom_state.N.L,
          # filename=fNames[1],
          #width = 12, height = 20,
          #cellwidth=18, cellheight=18,
          fontsize_row = 10,
          color=rwb,
          #border_color="black",
          treeheight_row = 0,
          show_rownames = TRUE,
          cluster_rows = FALSE,
          #scale = "row", 
          # clustering_distance_rows = "correlation",
          #kmeans_k = 10,
          #labels_col = FALSE,
          cluster_cols = FALSE,
          show_colnames = TRUE,
          #cutree_cols = 2,
          
          #cutree_rows =10,
          
          #clustering_distance_cols = "correlation",
          
          #color = diverge_hcl(100, "Red-Green",alpha = 0.8),
          #gaps_col = c(16),
          # color=rwb,
          #breaks = seq(-200, 200, length.out=151)
          # breaks = seq(quants[1], quants[2], length.out=101)
          breaks=breaks2
          #color = RdBu(100)
 )
 
 #############################################################
 
 

 ##### GET PEAK ANNOTATIONS e.g. chrom state or genomics region based on gene annotation (peakAtlas@elementMetadata$annot)
 
 my_row_annot=as.data.frame(rownames(peakQuants_single_sig_l2fc))
 
 rownames(my_row_annot)=my_row_annot$`rownames(peakQuants_single_sig_l2fc)`
 my_row_annot$anno_state=data_fram_overla_df_single$name
 my_row_annot$anno_type=data_fram_overla_df_single$annot
 my_row_annot=my_row_annot[-1]
 
 
 peaksQ_w_annot=cbind(peakQuants_single_sig_l2fc, my_row_annot$anno_type)
 colnames(peaksQ_w_annot)[33]="state"
 
 
 
 newcenters=peaksQ_w_annot %>% group_by(state) %>% summarise_all(funs(median))
 newcenters=as.data.frame(newcenters)
 rownames(newcenters)=newcenters$state
 peaksQnewcenters=newcenters[,2:33]
 
 
 
 
 
 pheatmap(peaksQnewcenters,
          fontsize_row = 10,
          treeheight_row = 0,
          show_rownames = TRUE,
          cluster_cols = FALSE,
          cluster_rows = TRUE,
          clustering_distance_rows = "correlation",
          show_colnames = TRUE,
          color = diverging_hcl(100, "Tropic",alpha = .7),
          gaps_col = c(16),
          breaks = seq(-.75, .75, length.out=101)
          
 )
 
 
