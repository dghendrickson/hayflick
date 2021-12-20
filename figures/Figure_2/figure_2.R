library(Seurat)
library(SeuratDisk)
library(pheatmap)
library(plyr)
library(dplyr)
library(colorspace)
library(reshape2)
library(devtools)
library(hdf5r)
library(BiocParallel) 
library(ggplot2) 
library(stringr) 
library(tibble)
library(data.table)


###load SC object

###download from GEO from 
S3_hay <- LoadH5Seurat("~/big_files_git/sceasy_S3hay.h5ad.h5seurat")

# set an unadulterated version
S3_hay.r=S3_hay

DimPlot(S3_hay.r)

### normalize with SCtransform
S3_hay <- SCTransform(object = S3_hay, verbose = TRUE)


### run PCA
S3_hay <- RunPCA(object = S3_hay, verbose = FALSE,assay = "SCT")


### run UMAP
S3_hay <- RunUMAP(object = S3_hay, verbose = TRUE,
                  n.components=2, dims = 1:50,
                  ,metric = "euclidean"
                  ,n.neighbors =50
                  ,min.dist = .3 ,reduction = "pca")

### flip UMAP coordinates to match previous implementation 

umap.cord=as.data.frame(S3_hay[["umap"]]@cell.embeddings)
umap.cord$UMAP_2=umap.cord$UMAP_2*-1
umap.cord=as.matrix(umap.cord)

S3_hay[["umap"]]@cell.embeddings=umap.cord

### check
DimPlot(S3_hay)

######CELL CYCLE scoring (seurat vignette) https://satijalab.org/seurat/articles/cell_cycle_vignette.html
cc.genes <- readLines(con = "~/hayflick/annotations/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]


S3_hay <- CellCycleScoring(S3_hay, s.features = s.genes, g2m.features = g2m.genes)

DimPlot(S3_hay, group.by = "Phase")

### set cell cycle colors

cc_cols=c('#1b9e77','#d95f02','#7570b3')



############# Clustering (used in pseudotime later)

S3_hay <- FindNeighbors(object = S3_hay, dims = 1:50, verbose = FALSE,assay = "SCT")


S3_hay <- FindClusters(object = S3_hay, verbose = FALSE,resolution = 1.3,algorithm = 2) 

DimPlot(S3_hay,label = TRUE, group.by = "seurat_clusters")

############################ Seurat Object just PDL exps 
tp_target=c("PDL_25","PDL_29","PDL_33","PDL_37","PDL_46","PDL_50")

#cl_target=c("0","1","2","3","4","5","6","7","8","9","10" ,"11","12")

TP.WT_S3_hay=S3_hay
Idents(TP.WT_S3_hay)=TP.WT_S3_hay$PDL
TP.WT_S3_hay=subset(TP.WT_S3_hay,idents=tp_target)




TP.WT_S3_hay <- SCTransform(object = TP.WT_S3_hay, verbose = TRUE)


### run PCA
TP.WT_S3_hay <- RunPCA(object = TP.WT_S3_hay, verbose = FALSE,assay = "SCT")


### run UMAP
TP.WT_S3_hay <- RunUMAP(object = TP.WT_S3_hay, verbose = TRUE,
                  n.components=2, dims = 1:50,
                  ,metric = "euclidean"
                  ,n.neighbors =50
                  ,min.dist = .3 ,reduction = "pca")


Idents(TP.WT_S3_hay)=TP.WT_S3_hay$PDL

##check order

unique(TP.WT_S3_hay$PDL)
TP.WT_S3_hay$PDL=as.character(TP.WT_S3_hay$PDL)

#####################################################################Figure 2A

DimPlot(TP.WT_S3_hay, group.by = "Phase" ,split.by = "PDL", cols =cc_cols)

############################################### 
unique(S3_hay$PDL)
S3_hay$PDL=as.character(S3_hay$PDL)

####### subset to just focus on 1st last tert to show no change over course of the experiment 

tp_target=c("PDL_25","PDL_29","PDL_33","PDL_37","PDL_46","PDL_50","htert_2","htert_7")

TP_S3_hay=S3_hay
Idents(TP_S3_hay)=TP_S3_hay$PDL
TP_S3_hay_fig=subset(TP_S3_hay,idents=tp_target)



################## make cell cycle proportions

tp=prop.table(x = table((TP_S3_hay_fig$Phase), (TP_S3_hay_fig$PDL)), margin = 2)     
tp=melt(tp)


tp = within(tp, {
  exp = ifelse(grepl("PDL", Var2) , "WT", "hTERT")
})
tp$exp=factor(tp$exp)
unique(levels(tp$exp))
tp$exp = factor(tp$exp,levels(tp$exp)[c(2,1)])
unique(levels(tp$exp))


#####################################################################Figure 2B

tp %>% 
  ggplot(aes(x = Var2,  y = value, fill = Var1))  +  
  geom_col(position = "fill",stat = "identity") +
  facet_grid(~exp, scales = "free_x", space = "free_x", switch = "x") + 
  scale_fill_manual(values = cc_cols) +
  theme(legend.position = "none") +
  
  theme(axis.line = element_line(size=.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 20))


###################################################### FIGURES 2C and 2D

DimPlot(TP_S3_hay_fig, pt.size = 0.1,group.by = "Phase", cols = cc_cols ) 
 


DimPlot(TP_S3_hay_fig, pt.size = 0.1,group.by = "batch2",order = c("H","G","F","C","B","A") ) + 
  scale_color_discrete_diverging(palette = "Blue-Red 2")


#####################################################################################
####################################################

################ 
############FIND sig Genes DEseq with increasing PDL in different phases of cell cycle using WALD test
##

####

##empty list for Deseq results
results=list()


Idents(TP.WT_S3_hay)=TP.WT_S3_hay$Phase
for (i in levels(Idents(TP.WT_S3_hay))) {                                 ##loop through cell cycle phase
  S_small=subset(TP.WT_S3_hay,  idents = c(i))
  c_count <- as(as.matrix(S_small@assays$SCT@counts), 'Matrix')              ##get counts for DEseq--I used SCT corrected--could aslo use raw in "RNA" slot
  c_count=as.data.frame(c_count)
  sumcheck=as.data.frame(rowSums(c_count))                                   ## data frame for filtering on 8000 most expressed genes-- more genes gets into very lowly expressed and starts to really slow down DEseq
  colnames(sumcheck)="sum"
  sumcheck$id=rownames(sumcheck)
  sumcheck=sumcheck[order(sumcheck$sum,decreasing = TRUE),]
  sumcheck=rownames(sumcheck[1:8000,])                                     ##filter top 8k genes
  sc_lookup=as.data.frame(names(c_count))
  colnames(sc_lookup)="barcode"
  tracking <- as.data.frame(str_split_fixed(sc_lookup$barcode, "-", 2))
  sc_lookup$PDL=tracking$V2
  sc_lookup$PDL=plyr::revalue(sc_lookup$PDL, c("1"="1","2"="2","3"="2",       ##collapse tech reps from each PDL
                                         "4"="3","5"="3","6"="4",
                                         "7"="4","8"="5","9"="5",
                                         "10"="6","11"="6"))
  sc_lookup$PDL=as.numeric(sc_lookup$PDL)
  rownames(sc_lookup)=sc_lookup$barcode
  
  count=(sc_lookup %>% dplyr::group_by(PDL) %>% dplyr::summarise(barcode=n()))    ## remove any phase by PDL partition with less than 15 member cells
  for(j in seq(count$PDL)) {
    if(count[j,2] < 15) 
      
      sc_lookup=sc_lookup %>% filter(!PDL == count$PDL[j])
    
  }
  
  c_count=c_count[rownames(c_count) %in% sumcheck,names(c_count) %in% sc_lookup$barcode] ## creat Phase by PDL subset matrix for DEseq
  
  ## DESEQ!
  dds <- DESeqDataSetFromMatrix(countData = c_count,
                                colData = sc_lookup,
                                design =  ~ PDL)
  
  
  dds_res=DESeq(dds,test = "Wald", parallel=TRUE, BPPARAM=MulticoreParam(10))    ## use bioCparallel to speed up
  
  results[[i]]=dds_res
  

###############################################################FIGURE 2E


###### CREATE LOG2 FC vs PDL25 CountsPerMillion (cpm) for all cells agreggated 

library(Seurat)
  
data_c <- as(as.matrix(TP.WT_S3_hay@assays$SCT@counts), 'Matrix')
data_c=as.data.frame(data_c)



PDL_meta=data.frame(cbind(as.character(TP.WT_S3_hay$barcode),as.character(TP.WT_S3_hay$PDL),as.character(TP.WT_S3_hay$PDL)),stringsAsFactors = FALSE)


TRD_pdl_squish=list()
for (i in (unique(PDL_meta$X2))) {
  
  ref="PDL_25"
  
  r1_list=PDL_meta %>% filter(X2 == ref)
  
  real_list=PDL_meta %>% filter(X2 == i)
  
  
  r1_sums=rowSums(data_c[,r1_list$X1])#+1
  r1_cpm <-(r1_sums/sum(r1_sums))*1000000 
  real_sums=rowSums(data_c[,real_list$X1])#+1
  real_cpm=(real_sums/sum(real_sums))*1000000 
  TRD_pdl_squish[[i]]=log2((real_cpm)/(r1_cpm))
}


PDL_cpms = do.call(cbind.data.frame, TRD_pdl_squish)
PDL_cpms=PDL_cpms[complete.cases(PDL_cpms),]
PDL_cpms <- PDL_cpms[!is.infinite(rowSums(PDL_cpms)),]



###### CREAT LOG2 FC vs PDL25 CountsPerMillion (cpm) for agregagetd cells BY CELL CYCLE PHASE and PDL

TP.WT_S3_hay <- FindNeighbors(object = TP.WT_S3_hay, dims = 1:50, verbose = FALSE,assay = "SCT")
TP.WT_S3_hay <- FindClusters(object = TP.WT_S3_hay, verbose = FALSE,resolution = .8,algorithm = 2)  

DimPlot(TP.WT_S3_hay,label = TRUE)

data_c <- as(as.matrix(TP.WT_S3_hay@assays$SCT@counts), 'Matrix')
data_c=as.data.frame(data_c)


DimPlot(TP.WT_S3_hay) ##UMAP cells used for 2E

DimPlot(TP.WT_S3_hay,group.by = "Phase") ##UMAP cells used for 2E colored by cell cycle Phase


cluster.PDL_meta=data.frame(cbind(as.character(TP.WT_S3_hay$barcode),as.character(TP.WT_S3_hay$Phase),as.character(TP.WT_S3_hay$PDL)),stringsAsFactors = FALSE)

cluster.PDL_meta=cluster.PDL_meta %>% filter(!X3=="PDL_50") ####remove PDL 50-- scored as "G1" but purpose of analysis is to focus on Sen signature in prescenescent cells

cluster_squish=list()
pdl_squish=list()
for (i in (unique(cluster.PDL_meta$X2))) {
  pdl_squish=list()
  c1_list=cluster.PDL_meta %>% filter(X2 == i)
  c1_list=c1_list[order(c1_list$X3),]
  for (j in (unique(c1_list$X3))) {
    
    p1_list=c1_list %>% filter(X3 == j)
    
    
   # if(dim(p1_list)[1] >30) {p1_sums=rowSums(data_c[,p1_list$X1])+1
    
    #} else {p1_sums=0}
    
   p1_sums=rowSums(data_c[,p1_list$X1])
    

    
    
    p1_cpm <-(p1_sums/sum(p1_sums))*1000000 
    
    if(length(pdl_squish) <1) {pdl_squish[[j]]=log2(p1_cpm)
    }else {pdl_squish[[j]] = unlist(lapply(list(log2(p1_cpm)), "-" , unlist(pdl_squish[1], use.names=FALSE)))}
    
    
  }
  
  cluster_squish[[i]]=pdl_squish
  rm(pdl_squish)
}

##########make big data table from list

fc_pdl_cluster=data.frame(matrix(, nrow=19807, ncol=0))
for (i in seq(length(cluster_squish))) {
  PDLC_cpms = do.call(rbind, cluster_squish[[i]])
  PDLC_cpms=as.data.frame(t(PDLC_cpms))
  #PDLC_cpms=PDLC_cpms[colSums(!is.na(PDLC_cpms)) > 0]
  
  colnames(PDLC_cpms) <- paste("C",names(cluster_squish[i]), colnames(PDLC_cpms), sep = "_")
  fc_pdl_cluster=cbind(fc_pdl_cluster,PDLC_cpms)
}

##remove NAs +/- INF
fc_pdl_cluster=fc_pdl_cluster[colSums(!is.na(fc_pdl_cluster)) > 0]

fc_pdl_cluster[sapply(fc_pdl_cluster, is.infinite)] <- NA

fc_pdl_cluster=fc_pdl_cluster[complete.cases(fc_pdl_cluster),]


colnames(fc_pdl_cluster)


############### Get Sig genes (by PDL) from each Phase and subset agreggated log2 CPM matrix
  c0_betas_sig_S=as.data.frame(results(results[["S"]])) %>% rownames_to_column('id') %>% filter(padj < .0001)
  c0_betas_sig_G1=as.data.frame(results(results[["G1"]])) %>% rownames_to_column('id') %>% filter(padj < .0001) 
  c0_betas_sig_G2M=as.data.frame(results(results[["G2M"]])) %>% rownames_to_column('id') %>% filter(padj < .0001) 
  
  c0_betas_sig=as.data.frame(unique(c(c0_betas_sig_S$id,c0_betas_sig_G1$id,c0_betas_sig_G2M$id)))
  colnames(c0_betas_sig)="gene"

  
  ccpdl_0=fc_pdl_cluster[rownames(fc_pdl_cluster) %in% c0_betas_sig$gene,]  #names(fc_pdl_cluster) %like% paste0("C_",i),]
  
  
  ccpdl_0[,names(ccpdl_0) %like% "*PDL_25"]=0 #set ref PDL to 0
  
  ccpdl_0_PDLtc=merge(ccpdl_0,PDL_cpms,by="row.names") ## merge with matrix all cells aggregated by PDL
  
  rownames(ccpdl_0_PDLtc)=ccpdl_0_PDLtc$Row.names
  
  ccpdl_0_PDLtc=ccpdl_0_PDLtc[-1]
  

  wt_cluster_avg_exp_mc =ccpdl_0_PDLtc
  
  rwb<-colorRampPalette(c("steel blue","white","tomato"))

  rwb=rwb(100)
  
  
pdf("~/hayflick/newSC_single.pdf", height=6, width=8)
 pheatmap::pheatmap(wt_cluster_avg_exp_mc,
           
           fontsize_row = 10,
           color=rwb,
           treeheight_row = 0,
           show_rownames = FALSE,
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cluster_cols = FALSE,
           show_colnames = TRUE,
           breaks = seq(-1.5, 1.5, length.out=101))
  dev.off()
  

  
  
  
  
  
  ########################################## FIGURE 2 SUP 2
  
  #############################################FIGURE 2 supp 4
  
  library(viridis)
  
  #######panel A
  
  RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = TRUE)
  
  DimPlot(TP_S3_hay_fig, pt.size = 0.1,group.by = "Phase", cols = cc_cols )  #from fig2 scripts
  
  sen.features=RS_deseq %>% filter(baseMean > 50) %>% filter(padj < .01) #filter out least expressed
  sen.features=sen.features[order(sen.features$log2FoldChange, decreasing = TRUE),] 
  sen.features=sen.features$gene 
  
  pro.features=RS_deseq %>% filter(padj < .01 & log2FoldChange < 0  & baseMean > 10)
  pro.features=pro.features[order(pro.features$log2FoldChange, decreasing = FALSE),]
  pro.featuresS=pro.features[1:500,1]
  
  
  
  ###### specific con score
  
  CD_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/CD_deseq.txt", stringsAsFactors = TRUE)
  
  
  con.features=CD_deseq %>% filter(baseMean > 50) %>% filter(padj < .01) #filter out least expressed
  con.features=con.features[order(con.features$log2FoldChange, decreasing = TRUE),] 
  con.features=con.features$gene
  
  pro.features=CD_deseq %>% filter(padj < .01 & log2FoldChange < 0  & baseMean > 10)
  pro.features=pro.features[order(pro.features$log2FoldChange, decreasing = FALSE),]
  pro.featuresC=pro.features[1:500,1]
  
  ####
  sen.features.spec=setdiff(sen.features,con.features)
  con.features.spec=setdiff(con.features, sen.features)
  
  TP_S3_hay_fig <- CellCycleScoring(TP_S3_hay_fig, s.features = sen.features.spec, g2m.features = pro.featuresS,)
  TP_S3_hay_fig$sen.score=TP_S3_hay_fig$S.Score
  
  TP_S3_hay_fig <- CellCycleScoring(TP_S3_hay_fig, s.features = con.features.spec, g2m.features = pro.featuresC,)
  TP_S3_hay_fig$con.score=TP_S3_hay_fig$S.Score
  
  
  
  plot1=FeaturePlot(TP_S3_hay_fig, "sen.score", max.cutoff = "q99.5",pt.size = .5) + scale_color_viridis(option = "magma")
  plot2=VlnPlot(TP_S3_hay_fig,"sen.score", group.by = "PDL" )
  
  plot3=FeaturePlot(TP_S3_hay_fig, "con.score",max.cutoff = "q99.5",pt.size = .5) + scale_color_viridis(option = "magma")
  plot4=VlnPlot(TP_S3_hay_fig,"con.score", group.by = "PDL" )
  
  pdf("~/hayflick/SEN.SCOREumap.pdf", height=6,width=7)
  print(plot1)
  dev.off()
  
  pdf("~/hayflick/SEN.SCORE.vln.pdf", height=6,width=7)
  print(plot2)
  dev.off()
  
  pdf("~/hayflick/con.SCORE.umap.pdf", height=6,width=7)
  print(plot3)
  dev.off()
  
  pdf("~/hayflick/con.SCORE.vln.pdf", height=6,width=7)
  print(plot4)
  dev.off()
  
  
  
  #############################################FIGURE 2 supp 3
  
  
  DimPlot(TP.WT_S3_hay, group.by = "Phase")
  
  Idents(TP.WT_S3_hay)=TP.WT_S3_hay$Phase
  phase=c("S","G2M")
  
  CC_hay=subset(TP.WT_S3_hay,idents=phase)
  CC_hay=SCTransform(CC_hay, verbose = TRUE,return.only.var.genes = FALSE,variable.features.n = 3000)
  CC_hay <- RunPCA(CC_hay,assay = "SCT")
  CC_hay <- RunUMAP(CC_hay, dims = 1:50, umap.method = "uwot",min.dist = .3,n.neighbors = 50)
  
  pdf("~/hayflick/CC_umap.pdf", height=6,width=6.5)
  
  DimPlot(CC_hay, group.by = "PDL",pt.size = 1) +  scale_color_discrete_diverging(palette = "Blue-Red 2")
  
  dev.off()
  
  
  
  writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_2/sessionInfo.txt")
  
  
  
  