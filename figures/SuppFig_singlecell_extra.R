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



#######panel B 

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

