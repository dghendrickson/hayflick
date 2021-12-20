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
library(monocle3)


#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))



#BiocManager::install(c("spdep"))
#install_github("r-spatial/sf")

#nlibrary(devtools)
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)

S3_hay <- LoadH5Seurat("~/big_files_git/sceasy_S3hay.h5ad.h5seurat")



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



S3_hay <- FindNeighbors(object = S3_hay, dims = 1:50, verbose = FALSE,assay = "SCT")


S3_hay <- FindClusters(object = S3_hay, verbose = FALSE,resolution = 1.3,algorithm = 2) 

DimPlot(S3_hay,label = TRUE, group.by = "seurat_clusters")


temp=data.frame(cbind(as.character(S3_hay$barcode),as.character(S3_hay$seurat_clusters)),stringsAsFactors = FALSE)
rownames(temp)=temp$X1
S3_hay=AddMetaData(S3_hay, temp$X2,"PS_seurat_clusters")




#target <- c("3", "9", "8","1","0","2","4","6")
target <- c("13", "21", "15","8","10","1","4","3","7","16","6") # isolate main WT trajectory 
Idents(S3_hay)=S3_hay$PS_seurat_clusters
Cs_WT.S3_hay=S3_hay
Cs_WT.S3_hay=subset(Cs_WT.S3_hay,idents=target)

Idents(Cs_WT.S3_hay)=Cs_WT.S3_hay$exp
target="WT"
Cs_WT.S3_hay=subset(Cs_WT.S3_hay,idents=target)

Cs_WT.S3_hay=SCTransform(Cs_WT.S3_hay, verbose = TRUE,return.only.var.genes = FALSE,variable.features.n = 3000)



Cs_WT.S3_hay <- RunPCA(Cs_WT.S3_hay,assay = "SCT")

Cs_WT.S3_hay <- RunUMAP(Cs_WT.S3_hay, dims = 1:50, umap.method = "uwot",min.dist = .3,n.neighbors =50,metric ="euclidean")

DimPlot(Cs_WT.S3_hay)
#########
cc.genes <- readLines(con = "~/hayflick/annotations/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]


Cs_WT.S3_hay <- CellCycleScoring(Cs_WT.S3_hay, s.features = s.genes, g2m.features = g2m.genes)

DimPlot(Cs_WT.S3_hay,group.by = "exp")


data <- as(as.matrix(Cs_WT.S3_hay@assays$SCT@counts), 'sparseMatrix')
pd <- (Cs_WT.S3_hay@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(expression_data = data,
                         cell_metadata = pd,
                         gene_metadata = fData)



############# get Genes

RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = FALSE)

RS_deseq=RS_deseq[order(RS_deseq$log2FoldChange, decreasing = TRUE),]

sen_genes_up= RS_deseq %>% filter(padj < 0.01)  %>% filter(log2FoldChange > 3)

sen_genes_down= RS_deseq %>% filter(padj < 0.01)  %>% filter(log2FoldChange < -3)


s.genes.4m_B=c(sen_genes_up$gene,sen_genes_down$gene)


sen_genes_up=sen_genes_up$gene

sen_genes_up=sen_genes_up[sen_genes_up %in% rowData(cds)@rownames]



s.genes.4m_B=read.delim("~/s.genes.4m_B.txt",stringsAsFactors = FALSE)
s.genes.4m_B=s.genes.4m_B$s.genes.4m_B

sen_genes_up=s.genes.4m_B[s.genes.4m_B %in% rowData(cds)@rownames]




#cds.5 <- preprocess_cds(cds, num_dim = 100,method = "PCA",residual_model_formula_str = "~S.Score" + "G2M.Score" ,use_genes = sen_genes_up)
 
   #cds.5 <- align_cds(cds.5,num_dim=100,residual_model_formula_str =  "~sphase_score + mphase_score")
   
 #  cds.5 <- reduce_dimension(cds.5,reduction_method = "UMAP",preprocess_method = "PCA",umap.min_dist = .05,umap.metric = "cosine",umap.n_neighbors = 50,)
 
 #plot_cells(cds.5,color_cells_by =  "PDL",cell_size = 1,label_groups_by_cluster = TRUE)

 
 #cds.5=cluster_cells(cds.5,resolution = .005)
 
 #cds.5=learn_graph(cds.5,,close_loop = FALSE,learn_graph_control = list(minimal_branch_len = 25,
  #                                                                      geodesic_distance_ratio = .5, euclidean_distance_ratio=1)) 
 #cds.5=order_cells(cds.5)
 
# plot_cells(cds.5,color_cells_by =  "pseudotime",cell_size = 1,show_trajectory_graph = TRUE)# +scale_color_viridis_c(begin = 0,end = 1)
 


 
 #colData(cds.5)$pseudotime <- pseudotime(cds.5)
 #meta=as.data.frame(colData(cds.5))
 
 #ggplot(meta, aes(x=pseudotime)) +
   #geom_histogram(bins =50)


 
 
 cds <- preprocess_cds(cds, num_dim = 100,method = "PCA" ,use_genes = sen_genes_up)
 cds <- reduce_dimension(cds,reduction_method = "UMAP",preprocess_method = "PCA",umap.min_dist = .05,umap.metric = "cosine",umap.n_neighbors = 50)
  plot_cells(cds,color_cells_by =  "PDL",cell_size = 1,label_groups_by_cluster = TRUE)
 
  plot_cells(cds,color_cells_by =  "PDL",cell_size = .75,show_trajectory_graph = FALSE) +
    scale_color_discrete_diverging(palette = "Blue-Red 2",rev = TRUE) + 
    theme(legend.position = "none") 
  
  cds=cluster_cells(cds,resolution = .005)
  
  
  cds=learn_graph(cds,close_loop = FALSE,learn_graph_control = list(minimal_branch_len = 35,
                                                                    geodesic_distance_ratio = .5, euclidean_distance_ratio=1)) 
  cds=order_cells(cds)
  
  
  ########### FIGURE 7A
  plot_cells(cds,color_cells_by =  "pseudotime",cell_size = 1,show_trajectory_graph = TRUE)# +scale_color_viridis_c(begin = 0,end = 1)
  
  
  
  
  
  #### FIGURE 7b
  
  
  
  AFD_genes <- c("CENPK","PAPPA","SNAI2")

  AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes]
  
  #The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
  
plot_genes_in_pseudotime(AFD_lineage_cds,
                              
                              color_cells_by="PDL",
                              min_expr=.5)
  
  
  
  
  
  
 ######################## 
  
  colData(cds)$pseudotime <- pseudotime(cds)
  meta=as.data.frame(colData(cds))
  
  ggplot(meta, aes(x=pseudotime)) +
    geom_histogram(bins =50)
  
  #######################
  
  subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
 
  
  subset_pr_test_res=subset_pr_test_res[order(subset_pr_test_res$morans_test_statistic,decreasing = TRUE),]
  
  
  pr_deg_ids_s= row.names(subset_pr_test_res[1:5000,])
  
  AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pr_deg_ids_s]
  
  ######
  colData(AFD_lineage_cds)$pseudotime <- pseudotime(AFD_lineage_cds )
  
  
  
  model_tbl = fit_models(AFD_lineage_cds , "~ splines::ns(pseudotime, df=3)")
  
  
  model_expectation <- model_predictions(model_tbl, new_data = colData(AFD_lineage_cds))
  
  colnames(model_expectation) <- colnames(AFD_lineage_cds)
  
########################  
  
  cds_subset=AFD_lineage_cds
  f_id <- NA
  Cell <- NA
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  
  min_expr <- 0
  
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", 
                     by.y = "row.names")

  cds_exprs <- base::merge(cds_exprs, cds_colData, by.x = "Cell", 
                           by.y = "row.names")
  cds_exprs$adjusted_expression <- cds_exprs$expression
  

  cds_exprs$f_id <- as.character(cds_exprs$f_id)

  
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                             function(x) {
                               data.frame(expectation = model_expectation[x$f_id, 
                                                                          x$Cell])
                             })
  
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  

  #####################################
  

  library(Hmisc)
  bins=60
  cds_exprs$ps_bin=cut_interval(cds_exprs$pseudotime, n=bins,labels=as.character(seq(1:bins)))

  df1= cds_exprs %>%
    select(c("f_id","expectation","ps_bin")) %>% 
    dplyr::group_by(ps_bin,f_id)
  
  ps_med <- dplyr::summarize(df1, med_exp_pred = median(expectation))  
  
  library(reshape2)
  
  ps_med_cast=ps_med %>% pivot_wider(names_from = ps_bin,values_from = med_exp_pred) %>%
    column_to_rownames( var = "f_id")
  
#  ps_med_cast = as.matrix(ps_med_cast)
  
  M_PC1_mono_WT=ps_med_cast

  ###center only  preserves magnitude
  # M_PC1_mono_WT = Matrix::t(scale(Matrix::t(M_PC1_mono_WT), center = TRUE,scale=FALSE))
  
  ## rescale 0 to 1
  ##   By default, this scales the given range of s onto 0 to 1, but either or both of those can be adjusted. For example, if you wanted it scaled from 0 to 10, 
  
  M_PC1_mono_WT= t(apply(M_PC1_mono_WT, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  
    M_PC1_mono_WT = M_PC1_mono_WT[is.na(row.names(M_PC1_mono_WT)) == FALSE, ]
  M_PC1_mono_WT[is.nan(M_PC1_mono_WT)] = 0
  heatmap_matrix <- M_PC1_mono_WT
 # heatmap_matrix=as.data.frame(heatmap_matrix)
  write.table(heatmap_matrix,"~/hayflick/PS_5kgenes_heatmap_matrix.txt",sep="\t" )
  ################################################
  
  
  heatmap_matrix=read.delim("~/hayflick/data/scRNAseq/PS_5kgenes_heatmap_matrix.txt", stringsAsFactors = FALSE)
  
  ###################################################################################
  ###########  NEW PS K MEDOID
  library(cluster)
  library(pamr)
  #BiocManager::install("hopach")
  library(hopach)
  heatmap_matrix=as.matrix(heatmap_matrix)
  
  distMat3 <- as.matrix( distancematrix(heatmap_matrix, d = "cosangle", na.rm = TRUE) )
  
  
  start_time <- Sys.time()
  
  
  PS_clusters_COS_k25 <- pam(distMat3, k = 25,   diss = TRUE, keep.diss = TRUE)
  
  heatmap_matrix=as.data.frame(heatmap_matrix)
  
  PS_df_clutsers=cbind(heatmap_matrix, PS_clusters_COS_k25$cluster)
  
  #  cluster_df_clutsers=cbind(cluster_df, cluster_df_CO_20$cluster)
  
  
  # cluster_df_clutsers=cbind(cluster_df, cluster_df_CO_euc$cluster)
  
  colnames(PS_df_clutsers)[dim(PS_df_clutsers)[2]]="cluster"
  
  table(PS_df_clutsers$cluster)
  

  og_ps=read.delim("~/hayflick/data/scRNAseq/Clusterlabels_PseudoTime.txt", stringsAsFactors = FALSE)  #Cluster solution from paper slightly different pam vs pamr versions cluster run
  og_ps=og_ps[-dim(og_ps)[2]]
  og_ps=og_ps %>% column_to_rownames("gene")
  cols <- c(colnames(og_ps)[!(colnames(og_ps) %in% c("cluster"))])
  
  
  gdf2 =og_ps %>% tidyr::pivot_longer(-cluster, names_to = "group", values_to = "expression")
  final_df2=gdf2 %>% dplyr::group_by(cluster, group) %>% dplyr::summarize(expression_mean = median(expression)) %>% tidyr::spread(., group, expression_mean)
  PS_df_clutsers_meds=as.data.frame(final_df2)
  rownames(PS_df_clutsers_meds) =PS_df_clutsers_meds$cluster
  PS_df_clutsers_meds=PS_df_clutsers_meds[-1]
  
  PS_df_clutsers_meds <- PS_df_clutsers_meds[,cols]
  
  
  rwb<-colorRampPalette(c("steel blue","white","tomato"))
  rwb=rwb(100)
  paletteLength=100
  

  PS_df_clutsers_meds=PS_df_clutsers_meds[-61]
  colnames(PS_df_clutsers_meds)=seq(1:60)
  
  maxes_order=colnames(PS_df_clutsers_meds)[max.col(PS_df_clutsers_meds,ties.method="first")]
  maxes_order=as.numeric(maxes_order)
  names(maxes_order)=seq(1:25)
  maxes_order=order(maxes_order)
  
  
  nc= sequential_hcl(100,rev = FALSE,alpha = 1,palette = "Inferno",power=1)
  nc=nc[20:length(nc)]
  
  
  
  nc= sequential_hcl(100,rev = TRUE,alpha = .9,palette = "Purple-Orange",power=1)
  nc=nc[1:(length(nc)-20)]
  
  pdf("/Volumes/GoogleDrive/Shared drives/Hayflick_paper/Figures/fig_chunks/PS_5kg_25k_cosangle_scale1_Purple.pdf")
  ph_cl=  pheatmap(PS_df_clutsers_meds[maxes_order,],
                   # filename=fNames[1],
                   #width = 12, height = 20,
                   #cellwidth=18, cellheight=18,
                   fontsize_row = 10,
                   # color=rwb,
                   border_color="black",
                   treeheight_row = 100,
                   #cluster_rows =  tf_tr,
                   cluster_rows = FALSE,
                   #  scale = "row",   
                   # clustering_distance_rows = "correlation",
                   #clustering_distance_rows = "euclidean",
                   cluster_cols = FALSE,
                   #clustering_distance_cols = "correlation",
                   #color = diverging_hcl(100, "Blue-Red"),
                   #gaps_col = c(18,48,66,69,99),
                   color=nc,
                   #  color=sequential_hcl(100, "Blue_Yellow",rev = TRUE,alpha = 1),
                   breaks = seq(0, 1, length.out = length(nc)+1)
                   # breaks = seq(quants[1], quants[2], length.out=101)
                   #breaks=myBreaks
                   #color = RdBu(100))
                   
  )
  
  
  
  
###############################################################  

  
  
  
  
  ####################LISA
  
  library(stringr)
  library(mat)
  
  RS_tpms=read.delim("~/hayflick/data/RNAseq/TPM_tables/RS_batch_tpms.txt", stringsAsFactors = FALSE)
  rownames(RS_tpms)=make.unique(RS_tpms$gene)
  RS_tpms=RS_tpms[-1]
  sen_TPMS=RS_tpms
  
  lisa_files=read.delim("~/hayflick/data/scRNAseq/LISA_results/PS_k25_samples.tab", stringsAsFactors = FALSE,header = FALSE)
  
  base="~/hayflick/data/scRNAseq/LISA_results/"
  fName <- vector()
  
  for (i in 1:length(lisa_files$V1)) { 
    fName[i] <- paste(base, lisa_files[i,],sep="")
    
  }
  cName=strsplit(fName,"/", 8)
  
  cluster_lab <- vector()
  for (i in 1:length(cName)) { 
    cluster_lab[i] <- cName[[i]][8]
    
  }
  #clus_temp=str_split_fixed(cluster_lab,"_",4)
  clus_temp=str_split_fixed(cluster_lab,"_",5)
  clus_lab=clus_temp[,4]
  
  clus_temp=str_split_fixed(clus_lab,"c",2)
  clus_lab=clus_temp[,2]
  
  
  lisa_results=data.frame()
  for (i in 1:length(fName)) {
    stderr=read.csv(file=fName[i], stringsAsFactors = FALSE)
    expthresh=as.matrix(sen_TPMS)
    expthresh=expthresh[which(rowMedians(expthresh) > 5) , ]
    stderr=stderr[stderr$Transcription.Factor %in% rownames(expthresh),1]
    

    stderr=data.frame(stderr,stringsAsFactors = FALSE)
    
   
    stderr$rank=((as.numeric(rownames(stderr))))
  
    stderr=stderr[order(stderr$stderr),]
  
    lisa_results[1:dim(stderr)[1],i]=(stderr$rank)
    rownames(lisa_results)=stderr$stderr

  }
  colnames(lisa_results)=clus_lab
  
  
  
  top_tfs=vector()
  for (i in 1:dim(lisa_results)[2]) { 
    mt <- lisa_results[order(lisa_results[,i]), ]
    mt=mt[1:10,]
    top_tfs=c(top_tfs,rownames(mt))
    
  }
  top_tfs=unique(top_tfs)
  
  
  length(top_tfs)
  
  
  lisa_results_small=lisa_results[rownames(lisa_results) %in% top_tfs ,]
  lisa_results_small=t(lisa_results_small) 

  lisa_results_small=as.data.frame(lisa_results_small)
  
  
  floor_u=100
  floor_d=1
  lisa_results_small[lisa_results_small>floor_u]=floor_u
  lisa_results_small[lisa_results_small<floor_d]=floor_d
  
  
  lisa_results_small$num=as.numeric(rownames(lisa_results_small))
  lisa_results_small=lisa_results_small[order(lisa_results_small$num,decreasing = FALSE),]
  x=dim(lisa_results_small)[2]
  lisa_results_small=lisa_results_small[-dim(lisa_results_small)[2]]
  
  lisa_results_small_scale=scale(lisa_results_small,center = TRUE,scale = TRUE)
  lisa_results_small_scale=lisa_results_small_scale[,which(colMins(lisa_results_small_scale) < -1.5) ] ### specific filter
 
  library(dplyr)
  
 # tf_sp <- hclust(as.dist(1-cor(t(lisa_results_small_scale), method="kendall")), method="complete") 
 # tf_sp_c <- hclust(as.dist(1-cor((lisa_results_small_scale), method="kendall")), method="complete") 
  
  tf_sp_c <- hclust(as.dist(1-cor((lisa_results_small_scale), method="spearman")), method="complete") 
  
  
  nc= sequential_hcl(100,rev = FALSE,alpha = 1,palette = "Reds 3",power=4)
  nc=nc[15:length(nc)]
  
  #lisa_results_small_scale=lisa_results_small_scale[c(3,10,14,20,21,25,13,5,11,16,18,23,19,8,24,1,4,6,17,2,7,15,22,9,12),]
  
  #lisa_results_small_scale=lisa_results_small_scale[maxes_order,]
  
  #reorg=c("MITF","NRF1","KDM2B","ZNF282","RARA","ZNF202","ERCC6","TRIM28","RUNX1","MAZ","SP1","MXI1","NFKB1","FLI1","CCNT2","REST","RELA","KMT2A","EGR1","CREB1")
  
  pdf("/Volumes/GoogleDrive/Shared drives/Hayflick_paper/Figures/fig_chunks/reorg_PS_autoclust_TFs.pdf", width=10,height = 6)
  ph_tf=pheatmap((lisa_results_small_scale[maxes_order,]),#angle_col = 45,
                 # filename=fNames[1],
                 # width = 12, height = 25,
                 #  cellwidth=11, cellheight=35,
                 fontsize_row = 10,               
                 #color=rwb,
                 # border_color="black", 
                 treeheight_col =100,
                 #treeheight_row = 0,
                 #scale = "column",
                 cluster_rows = FALSE,
                 #cluster_rows = tf_sp,
                 # clustering_distance_rows = "correlation",
                 # clustering_distance_rows = "euclidean",
                 # cluster_cols = tf_sp_c,
                 # clustering_distance_cols = "euclidean",
                 cluster_cols = tf_sp_c,
                 #clustering_distance_cols = "euclidean",
                 #color = diverging_hcl(100, "Blue-Red"),
                 #gaps_col = c(18,48,66,69,99),
                 #color=r2,
                 #cutree_rows =5,
                 #cutree_cols = 9,   
                 # color = sequential_hcl(100,rev = FALSE,alpha = .9,h1=10,c1=65,l1=20,l2=97,cmax=150,power=3,fixup = TRUE),
                 #color = sequential_hcl(100,rev = FALSE,alpha = .9,palette = "Reds 3",power=4),
                 color=nc,
                 # color=magma(40, direction = -1),
                 breaks = seq(-1.5,1.5, length.out=length(nc)+1)
                 #breaks = seq(0,75, length.out=101)
                 #  breaks = seq(quants[1], quants[2], length.out=101)
                 #breaks=myBreaks
                 #color = RdBu(100)
  )
  dev.off()
  tf_tr=ph_tf$tree_row  
  
  
  
  #################### GO TERMS
  ()
  install.packages("gprofiler2")
  library(gprofiler2)
  
  #lapply(pathways[12254:12290], write, "test.txt", append=TRUE, ncolumns=1000)
  #lapply(names(pathways[12254:12290]), write, "names.txt", append=TRUE, ncolumns=1000)
  
  c_GMT= upload_GMT_file(gmtfile = "~/hayflick/annotations/custom_GMT.gmt") 
  
  custom_bg=RS_deseq$gene
  
  
  
  
  PS_df_clutsers$id=rownames(PS_df_clutsers)
  gp_result_list_PS=list()
  gp_custom_result_list_PS=list()
  for (cl in seq(1:length(unique(PS_df_clutsers$cluster)))) {
    
    iso_clust=PS_df_clutsers %>% dplyr::filter(PS_df_clutsers$cluster == cl)  %>% dplyr::select(id) 
    

    gostres <- gost(query = iso_clust$id,
                    organism = "hsapiens", ordered_query = FALSE, custom_bg = custom_bg, domain_scope = "custom",evcodes = TRUE,significant = FALSE,sources=c("GO","MIRNA"))
    
    
    
    x=as.data.frame(gostres$result)
    gp_result_list_PS[[cl]]=x[,3]
    names(gp_result_list_PS[[cl]])=x[,11]
    
    go_custom<- gost(query = iso_clust$id,
                     organism = c_GMT, ordered_query = FALSE, custom_bg = custom_bg, domain_scope = "custom",evcodes = TRUE,significant = FALSE)
    
    
    
    
    y=as.data.frame(go_custom$result)
    gp_custom_result_list_PS[[cl]]=y
    
    gp_custom_result_list_PS[[cl]]=y[,3]
    names(gp_custom_result_list_PS[[cl]])=y[,9]
    
    
    
  }
  
  
  
  writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_7/sessionInfo.txt")
  
  
  
  
  
  
                               