

VPdata=read.delim("~/hayflick/data/RNAseq/TPM_tables/VP_tpms_cb.txt",stringsAsFactors = FALSE)



library(tidyverse)
library(pheatmap)
library(ggplot2)
library(reshape2)

#setwd("~/YAP_RNAseq/YAPi_compound_validation_RNAseq/")
VPdata$gene=make.unique(VPdata$gene)

TPM_matrix = VPdata %>% column_to_rownames("gene")
TPM_matrix=TPM_matrix+1
notrtment = TPM_matrix %>% select( contains("X2hrs_0uM")) %>% rowMeans()

TPM_matrix_fc=TPM_matrix/notrtment

TPM_matrix_fc_l2=log2(TPM_matrix_fc)

#################genes myofivs



lit_gene=c("ACTA2","VIM","DES","PALLD","TIMP1","FN1","SPARC","RUNX1","PPIC","POSTN","FSTL1","FBN1","CTHRC1","CALD1",
           "COL5A2","COL5A1","COL3A1","COL1A2","COL4A6","COL4A5","COL4A4","COL4A2","COL4A1","COL1A1",
           "WNT5A","THBS1","TGFBI","TGFB3","TGFB2","TGFB1")


lit_gene=c("ACTA2","VIM","DES","PALLD","TIMP1")



lit_gene=c("FN1","SPARC","RUNX1","PPIC","POSTN","FSTL1","FBN1","CTHRC1","CALD1")

VPdata_l2fc_small=TPM_matrix_fc_l2[rownames(TPM_matrix_fc_l2) %in% lit_gene,]

VP.yap.merge=merge(sen_TPMS_l2fc,VPdata_l2fc_small,by="row.names")

VP.yap.merge=VP.yap.merge %>% select(!contains("PDL20")) %>% select(!contains("0hrs")) %>% select(!contains("2hrs_0"))


rownames(VP.yap.merge)=VP.yap.merge$Row.names
VP.yap.merge=VP.yap.merge[-1]



rwb<-colorRampPalette(c("steel blue","white","tomato"))
rwb=rwb(100)

pdf("~/hayflick/VP_NF_HM.pdf",height = 6,width = 6)
pheatmap(VP.yap.merge,
         
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames =TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         breaks = seq(-2, 2, length.out=101))
dev.off()
############################################# YAP

yap_zhang_jbc=read.delim("~/hayflick/annotations/published_yap_targets/yap_zhang_jbc_2009.txt", stringsAsFactors = FALSE)
yap_zhang_jbc=yap_zhang_jbc$TEAD.dependent.YAP.targets.from.Zhang.et.al..Journal.of.Biological.Chemistry.2009.284.13355.13362

yap_dupont=read.delim("~/hayflick/annotations/published_yap_targets/yap_dupont_et_al_2018.txt", stringsAsFactors = FALSE)
yap_dupont=yap_dupont$YAP.signature.from.Dupont.et.al..Nature.2018.474.179.183


yap_cordendossi_msigDB=read.delim("~/hayflick/annotations/published_yap_targets/yap_cordendossi_msigDB_2011.txt", stringsAsFactors = FALSE)
yap_cordendossi_msigDB=yap_cordendossi_msigDB$YAP.signature.from.Cordenonsi.et.al..Cell.2011.147.759.772..CORDENONSI_YAP_CONSERVED_SIGNATURE.in.MSigDB

yap_wang=read.delim("~/hayflick/annotations/published_yap_targets/yap_s_wang_2018.txt", stringsAsFactors = FALSE)
yap_wang=yap_wang$YAP.signature.from.Wang.et.al..Cell.Reports.2018.25.1304.1317

yap_zhang_msigdb=read.delim("~/hayflick/annotations/published_yap_targets/yap_zhang_et_al_msigDB_yap_up.txt", stringsAsFactors = FALSE)
yap_zhang_msigdb=yap_zhang_msigdb$YAP.signature.from.Zhang.et.al..Cancer.Research.2008.68.2789.2794..YAP1_UP.in.MSigDB

yap_all=c(yap_zhang_jbc,yap_dupont,yap_cordendossi_msigDB,yap_wang,yap_zhang_msigdb)


oc=table(yap_all)

oc=as.data.frame(oc)
oc_trd=oc %>% filter(Freq >1)


VPdata_l2fc_small=TPM_matrix_fc_l2[rownames(TPM_matrix_fc_l2) %in% oc_trd$yap_all,]



rwb<-colorRampPalette(c("steel blue","white","tomato"))
rwb=rwb(100)


VP.yap.merge=merge(sen_TPMS_l2fc_small,VPdata_l2fc_small,by="row.names")

VP.yap.merge=VP.yap.merge %>% select(!contains("PDL20")) %>% select(!contains("0hrs")) %>% select(!contains("2hrs_0"))


rownames(VP.yap.merge)=VP.yap.merge$Row.names
VP.yap.merge=VP.yap.merge[-1]

pdf("~/hayflick/yap_oc_VP.pdf",height=6,width=6)
pheatmap(VP.yap.merge,
         
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         breaks = seq(-2, 2, length.out=101))
dev.off()










VPdata_l2fc_small=TPM_matrix_fc_l2[rownames(TPM_matrix_fc_l2) %in% ids,]



VP.yap.merge=merge(sen_TPMS_l2fc,VPdata_l2fc_small,by="row.names")

VP.yap.merge=VP.yap.merge %>% select(!contains("PDL20")) %>% select(!contains("0hrs")) %>% select(!contains("2hrs_0"))


rownames(VP.yap.merge)=VP.yap.merge$Row.names
VP.yap.merge=VP.yap.merge[-1]

pdf("~/hayflick/EMT_VP.pdf",height=6,width=6)
pheatmap(VP.yap.merge,
         
         fontsize_row = 10,
         color=rwb,
         treeheight_row = 0,
         show_rownames = FALSE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = FALSE,
         show_colnames = TRUE,
         breaks = seq(-2, 2, length.out=101))
dev.off()





########################
all_protein_data_small=all_protein_data[rownames(all_protein_data) %in% oc_trd$yap_all,]





############


nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=1)
breaks = seq(-2 ,2, length.out=101)
breaks2 = breaks

breaks2[length(breaks)] <- max(max(all_protein_data[,dim(all_protein_data)[2]]),max(breaks))
breaks2[1] <- min(min(all_protein_data[,dim(all_protein_data)[2]]),min(breaks))

all_protein_data_small=all_protein_data_small %>% select(!contains("TERT"))


pheatmap(all_protein_data_small,
         fontsize_row = 10,
         color=nc,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         
         clustering_distance_rows = "euclidean",
         
         cluster_cols = FALSE,
         show_colnames = TRUE,
         
         breaks = breaks2)



















#######################

all_protein_data_small=all_protein_data[rownames(all_protein_data) %in% lit_gene,]




nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=1)
breaks = seq(-2 ,2, length.out=101)
breaks2 = breaks

breaks2[length(breaks)] <- max(max(all_protein_data[,dim(all_protein_data)[2]]),max(breaks))
breaks2[1] <- min(min(all_protein_data[,dim(all_protein_data)[2]]),min(breaks))

all_protein_data_small=all_protein_data_small %>% select(!contains("TERT"))

pdf("~/hayflick/protein_myf_lit_HM.pdf", height=5,width=5)
pheatmap(all_protein_data_small,
         fontsize_row = 10,
         color=nc,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = FALSE,
         
         clustering_distance_rows = "euclidean",
         
         cluster_cols = FALSE,
         show_colnames = TRUE,
         
         breaks = breaks2)
dev.off()


all_protein_data_small=all_protein_data[rownames(all_protein_data) %in% oc_trd$yap_all,]



nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=1)
breaks = seq(-2 ,2, length.out=101)
breaks2 = breaks

breaks2[length(breaks)] <- max(max(all_protein_data[,dim(all_protein_data)[2]]),max(breaks))
breaks2[1] <- min(min(all_protein_data[,dim(all_protein_data)[2]]),min(breaks))

all_protein_data_small=all_protein_data_small %>% select(!contains("TERT"))

pdf("~/hayflick/protein_myf_YAP_HM.pdf", height=5,width=5)
pheatmap(all_protein_data_small,
         fontsize_row = 10,
         color=nc,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = FALSE,
         
         clustering_distance_rows = "euclidean",
         
         cluster_cols = FALSE,
         show_colnames = TRUE,
         
         breaks = breaks2)
dev.off()

