library(stringr)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(colorspace)
library(tibble)
library(lubridate)
library(broom)
library(plyr)

library("devtools")
#install_github("jdstorey/qvalue")
library("qvalue")
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(purrr)
library(qvalue)

library(fgsea)
###########load data

pdl_protein_data_l2fc=read.delim("~/hayflick/data/proteomics/PDL_average_replicates_l2fc.txt", stringsAsFactors = FALSE)

tert_protein_data_l2fc=read.csv("~/hayflick/data/proteomics/tert_average_replicates_l2fc.csv", stringsAsFactors = FALSE)
tert_protein_data_l2fc=tert_protein_data_l2fc %>% dplyr::select(contains("TERT" ))

all_protein_data=cbind(pdl_protein_data_l2fc,tert_protein_data_l2fc)

all_protein_data=all_protein_data[complete.cases(all_protein_data),]




############# Figure 3 supp fig 1 


nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=1)
breaks = seq(-3 ,3, length.out=101)
breaks2 = breaks

breaks2[length(breaks)] <- max(max(all_protein_data[,dim(all_protein_data)[2]]),max(breaks))
breaks2[1] <- min(min(all_protein_data[,dim(all_protein_data)[2]]),min(breaks))



pheatmap(all_protein_data[,2:dim(all_protein_data)[2]],
         fontsize_row = 10,
         color=nc,
         treeheight_row = 0,
         show_rownames = FALSE,
         cluster_rows = TRUE,
      
         clustering_distance_rows = "euclidean",
   
         cluster_cols = FALSE,
         show_colnames = TRUE,
      
         breaks = breaks2)
      

################################################### GSEA (FIGURE 3A)

pathways2 <- gmtPathways("~/hayflick/annotations/h.all.v7.0.symbols_.gmt")

rownames(all_protein_data)=make.unique(all_protein_data$Gene)


all_protein_data=all_protein_data[!(is.na(all_protein_data$Gene) | all_protein_data$Gene==""), ] #remove rows w/o gene names


all_protein_data=all_protein_data[-1]

fullGSEA <- vector("list", length=dim(all_protein_data)[2])

rank_list<- list()

for ( i in 1:(dim(all_protein_data)[2])){
  rank_list[[i]]= setNames(all_protein_data[,i], rownames(all_protein_data))
}

names(rank_list)=colnames(all_protein_data[1:dim(all_protein_data)[2]])

names(fullGSEA) <-colnames(all_protein_data)

### Loop for running GSEA ON ALL CONDITIONS IN rank list LIST and GSEA list

for (set in 1:length((rank_list))) {
  print(names(rank_list)[set])
  ranks <- rank_list[[set]]
  
  fgseaRes <- fgseaMultilevel(pathways2, ranks ,minSize=10, maxSize=1500,eps = 0)
  fullGSEA[[set]] <- fgseaRes
}
##### UNADULTERATED RAW GSEA OUTPUT--all annoatations



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
# for (l in 1:length(fullGSEA)) {
#matches <- match(pvalGSEA$pathway, fullGSEA[[l]]$pathway)
# pvalGSEA[,l+1] <- fullGSEA[[l]]$pval[matches]
#pvalGSEA[,l+1] <- fullGSEA[[l]]$padj[matches]
#  betaGSEA[,l+1] <- fullGSEA[[l]]$NES[matches]
# }


pvalGSEA[,2:ncol(pvalGSEA)] <- -log10(pvalGSEA[,2:ncol(pvalGSEA)]) * sign(betaGSEA[,2:ncol(betaGSEA)])


pvalGSEA2 <- as.matrix(pvalGSEA[,2:ncol(pvalGSEA)])

rownames(pvalGSEA2) <- pvalGSEA[,1]


betaGSEA2 <- as.matrix(betaGSEA[,2:ncol(betaGSEA)])
rownames(betaGSEA2) <- betaGSEA[,1]
#betaGSEA2=betaGSEA2[-1]

pvalGSEA2=pvalGSEA2[complete.cases(pvalGSEA2),]
betaGSEA2=betaGSEA2[rownames(betaGSEA2) %in% rownames(pvalGSEA2),]



pcut <- -log10(.01)



library(reshape2)
tmp <- reshape2::melt(apply(pvalGSEA2, 1, function(x){length(which(abs(x)>=pcut))}))
##i verified that the rows are in the same order as the starting matrix, so I'm just going to do direct indexing removal
pvalGSEA2 <- pvalGSEA2[-which(tmp$value==0),]
betaGSEA2 <- betaGSEA2[-which(tmp$value==0),]



nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=2)
breaks = seq(-4 ,4, length.out=101)
#breaks2 = breaks

#breaks2[length(breaks)] <- max(max(pvalGSEA2),max(breaks))
#breaks2[1] <- min(min(pvalGSEA2),min(breaks))



pheatmap(pvalGSEA2,
         fontsize_row = 10,
         color=nc,
         treeheight_row = 0,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         
         clustering_distance_rows = "euclidean",
         
         cluster_cols = FALSE,
         show_colnames = TRUE,
         
         breaks = breaks)



#################################### Pull in TPMs

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


columns= c("eneg","name","value","group","condition","numb") 

long_Data = tibble(matrix(nrow = 0, ncol = length(columns))) 

# assign column names
colnames(long_Data) = columns

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

#dim(TPM_FC_MERGE)


######
gene="NNMT"
gene.use=long_Data %>% filter(eneg==gene) #%>% filter(condition=="RS")
gene.use$numb=factor(gene.use$numb)
gene.use$condition=factor(gene.use$condition, levels = c("TERT","RS","RIS","CD"))
pdfname=paste0("~/hayflick/",gene,"_box.pdf")

paper_cols_builtin_alpha=c("#6EBF91","#F4797E","#FABE79","#70ADD7")
paper_cols=c("#099146","#E8282B","#F39123","#1B76BA")

dodge <- position_dodge(width = 0.9)


######################### BOX PLOT VERSION

plot= ggplot(gene.use,aes(x=numb , y=value),color=condition) +
  scale_x_discrete() +
  geom_jitter(aes( x = numb), 
              position = position_jitter(width = .05), alpha = 1,size=1) +
  facet_wrap(~ condition, nrow = 1,scales = "free_x",) +
  scale_color_manual(values = paper_cols) +
  scale_fill_manual(values = paper_cols) +
  geom_boxplot(aes(x = numb,color=condition, fill=condition), outlier.colour = NA, position = dodge,alpha=.5) +

  theme_bw() +
  theme(legend.position = "none")

pdf(pdfname,height=5,width=7)
plot
dev.off()


######################### BARGRAPH PLOT VERSION

gene=c("SLC2A1","FOXE1","NNMT","CD36")

gene=c("FOXE1")


for (i in gene){
gene.use=long_Data %>% filter(eneg==i) #%>% filter(condition=="RS")
gene.use$numb=factor(gene.use$numb)
gene.use$condition=factor(gene.use$condition, levels = c("TERT","RS","RIS","CD"))
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


##############################  LOAD IN METABALOMICS AND SIG TESTING

#read in QE1 and QE2 RAW data

Qe1.2_raw=read.delim("~/hayflick/data/metabalomics/raw_Qe1_QE2.txt", stringsAsFactors = FALSE)
rownames(Qe1.2_raw)=Qe1.2_raw$labelcompound_QE2_parent_medRt

Qe1.2_raw=Qe1.2_raw[-2]



################
#######################  SPLIT tert and PDL timecourse
Qe1.2_tert =Qe1.2_raw %>%  dplyr::select(contains("hT"))
Qe1.2_tert=Qe1.2_tert+1 #add pseudocount

Qe1.2_pdl= Qe1.2_raw%>%  dplyr::select(contains("WI"))
Qe1.2_pdl=Qe1.2_pdl+1 #add pseudocount

##################### NORMALIZE METAB DATA TO TOTAL PROTEIN CONTENT

##median BCA values for tert samples HARD CODED IN
htert_pc=c(rep(x = 1.01589452,4),rep(x = 1.073200718,4),rep(x = 1.000611485,4),rep(x = 0.848035814,4),rep(x = 1.104643641,3))

##median BCA values for PDL samples HARD CODED IN
wi_pc=c(rep(x = 1.219290158,4),rep(x = 1.126617193,4),rep(x = 1.122333678,4),rep(x = 0.918550028,4),rep(x = 0.819630173,4),rep(x = 0.915329854,2))

##norm values bt protein BCA
Qe1.2_tert=cbind(Qe1.2_raw[1],sweep(Qe1.2_tert, MARGIN=2, htert_pc, `*`)) #for tert
colnames(Qe1.2_tert)[1]="CID"

Qe1.2_pdl=cbind(Qe1.2_raw[1],sweep(Qe1.2_pdl, MARGIN=2, wi_pc, `*`))#for pdl
colnames(Qe1.2_pdl)[1]="CID"


##### significance testing of slope over time (linear in log space-- so exponential change given biggest changes clearly at PDL50) We did in linear space, r.quared better in log

metab_pivot=Qe1.2_pdl %>% pivot_longer(-CID)
temp=str_split_fixed(metab_pivot$name,"_",2)
metab_pivot$tp=temp[,1]
metab_pivot$value=log2(metab_pivot$value)
temp=str_split_fixed(metab_pivot$CID,"_",2)
metab_pivot$ID=temp[,1]
metab_pivot$exp="pdl"

long_metab=metab_pivot

metab_pivot=Qe1.2_tert %>% pivot_longer(-CID)
temp=str_split_fixed(metab_pivot$name,"_",2)
metab_pivot$tp=temp[,1]
metab_pivot$value=log2(metab_pivot$value)
temp=str_split_fixed(metab_pivot$CID,"_",2)
metab_pivot$ID=temp[,1]
metab_pivot$exp="tert"

long_metab=rbind(long_metab,metab_pivot)

metab="L-Cysteine"
metab="Quinoline"
metab= "1-Methylnicotinamide"
metab= "Oxoglutaric acid"
metab="trans-4-Hydroxy-L-proline"
panel1=c("Uridine diphosphate glucose",
         "Uridine diphosphate-N-acetylglucosamine","FRUCTOSE 6-PHOSPHATE",
         "GLYCERALDEHDYE-3-PHOSPHATE","Glycerol 3-phosphate",
         "L-Lactic acid",
         "ALPHA-D-GLUCOSE 1-PHOSPHATE","2-Phosphoglyceric acid", "Pyruvic acid",
         "RIBOSE 5-PHOSPHATE",
         "CARNITINE","DEOXYCARNITINE","L-Acetylcarnitine",
         "LysoPC(16:0)","LysoPC(18:1(9Z))","LysoPE(16:0/0:0)","LysoPE(18:0/0:0)",
         "CITRIC ACID","ISOCITRIC ACID","Oxoglutaric acid","L-Malic acid","Fumaric acid")


metab_pivot_small=long_metab %>% filter(ID == metab)

ggplot(metab_pivot_small,aes(x=tp , y=value)) +
  scale_x_discrete() +
  geom_jitter(aes( x = tp,color=exp), 
              position = position_jitter(width = .05), alpha = 1,size=1) +
  #facet_wrap(~ condition, nrow = 1,scales = "free_x",) +
  #scale_color_manual(values = paper_cols) +
  #scale_fill_manual(values = paper_cols) +
  geom_boxplot(aes(x = tp,color=exp), outlier.colour = NA, position = dodge,alpha=.5) +
  
  theme_bw() +
  theme(legend.position = "none")
##############################################




##########################
###remove PDL 28--no htert to norm to

Qe1.2_pdl_fc=Qe1.2_pdl %>% dplyr::select(!contains("PDL28"))

#Qe1.2__gdf1_tert

#Qe1.2_tert

Qe1.2_pdl_fc=Qe1.2_pdl_fc[order(match(rownames(Qe1.2_pdl_fc) ,Qe1.2_tert$labelcompound_QE2_parent_medRt)),]

pdls=c("PDL25","PDL33","PDL37","PDL46","PDL50")
build_list=list()
for (i in pdls){

  temp_df=Qe1.2_pdl_fc %>% dplyr::select(contains(i))
 
   temp_tert_ref= as.data.frame(Qe1.2_tert) %>% dplyr::select(contains(i))
  temp_df=sweep(temp_df, MARGIN=1, temp_tert_ref[,1], `/`) #for tert
 build_list[[i]]=temp_df
}

metab_raw_fc_tert_norm=do.call(cbind.data.frame, build_list)


##clean up
metab_raw_fc_tert_norm=metab_raw_fc_tert_norm[complete.cases(metab_raw_fc_tert_norm),]
metab_raw_fc_tert_norm <- metab_raw_fc_tert_norm[!is.infinite(rowSums(metab_raw_fc_tert_norm)),]

metab_raw_fc_tert_norm_2=metab_raw_fc_tert_norm
metab_raw_fc_tert_norm_2$CID=rownames(metab_raw_fc_tert_norm)


regressions.pdl.tertcor <- metab_raw_fc_tert_norm_2 %>%
  gather(v, value, -CID,) %>%
  group_by(CID) %>% 
  nest() %>%
  mutate(data = map(data, ~ .x %>% mutate(pdl_n = c(rep(.16,4),rep(.43,4),rep(.56,4),rep(.86,4),rep(1,2))))) %>% # add numeric value for PDL or tert TP to percent of timecourse between PDL20 to PDL 50
  mutate(
    lin.fit = map(data, ~lm(value ~ pdl_n, data = .)),
    lin.tidied = map(lin.fit, tidy),
    lin.glanced = map(lin.fit, glance))  %>%
  mutate(
    log.fit = map(data, ~lm(log(value) ~ pdl_n, data = .)),
    log.tidied = map(log.fit, tidy),
    log.glanced = map(log.fit, glance)) #%>%


##### object with slopes and pvalues for each metabolite
regressions.pdl.tertcor= regressions.pdl.tertcor %>%
  unnest(lin.tidied,names_sep = ".")  %>%
  filter(lin.tidied.term == "pdl_n")

#### get fdr values

pvalues=regressions.pdl.tertcor$lin.tidied.p.value
qobj.tert <- qvalue(p = pvalues)
regressions.pdl.tertcor$fdr=qobj.tert$lfdr
regressions.pdl.tertcor$qval=qobj.tert$qvalues

regressions.pdl.tertcor=regressions.pdl.tertcor %>% dplyr::select(c("CID","lin.tidied.estimate","lin.tidied.std.error","lin.tidied.statistic","lin.tidied.p.value","fdr","qval"))


sig.metab.pdl.COR=regressions.pdl.tertcor %>% filter(qval < 0.05)  %>% pull(CID)


############################################################ FOLD CHANGES FOR PLOTTING
#metab_raw_fc_tert_norm

fc_pdl25_rowM=metab_raw_fc_tert_norm %>% select(contains("PDL25"))  %>% rowMeans()
FC_metab=sweep(metab_raw_fc_tert_norm, MARGIN=1, fc_pdl25_rowM, `/`) 
FC_metab=log2(FC_metab)

FC_metab$id=rownames(FC_metab)



FC_metab_gdf1 = gather(FC_metab, "group", "Expression",-id)

FC_metab_gdf1$tgroup = (str_split_fixed(FC_metab_gdf1$group, "[.]", 2)[, c(1)])

FC_metab_medl2fc=FC_metab_gdf1 %>% dplyr::group_by(id, tgroup) %>% 
  dplyr::summarize(expression_mean = median(Expression)) %>% 
  tidyr::spread(., tgroup, expression_mean)

temp=str_split_fixed(FC_metab_medl2fc$id,"_",2)
FC_metab_medl2fc$common_id=temp[,1]

#FC_metab_plot_sig=FC_metab_medl2fc[FC_metab_medl2fc$id %in% sig.metab.pdl.COR,]

panel1=c("Uridine diphosphate glucose",
         "Uridine diphosphate-N-acetylglucosamine","FRUCTOSE 6-PHOSPHATE",
         "GLYCERALDEHDYE-3-PHOSPHATE","Glycerol 3-phosphate",
         "L-Lactic acid",
         "ALPHA-D-GLUCOSE 1-PHOSPHATE","2-Phosphoglyceric acid", "Pyruvic acid",
         "RIBOSE 5-PHOSPHATE",
         "CARNITINE","DEOXYCARNITINE","L-Acetylcarnitine",
         "LysoPC(16:0)","LysoPC(18:1(9Z))","LysoPE(16:0/0:0)","LysoPE(18:0/0:0)",
         "CITRIC ACID","ISOCITRIC ACID","Oxoglutaric acid","L-Malic acid","Fumaric acid")

##Figure 4
panel1=c("1-Methylnicotinamide",
         "L-Methionine",
         "S-Adenosylhomocysteine",
         "S-Adenosylmethionine",
         "NAD",
         "NAD+",
         "NICOTINAMIDE")



FC_metab_medl2fc_plot = FC_metab_medl2fc[FC_metab_medl2fc$common_id %in% panel1,]
FC_metab_medl2fc_plot= FC_metab_medl2fc_plot %>% column_to_rownames("id")
FC_metab_medl2fc_plot=FC_metab_medl2fc_plot[-6]

nc= diverge_hcl(100,rev = FALSE,alpha = .9,palette = "Green-Orange",power=1)
breaks = seq(-1.5 ,1.5, length.out=101)

plot=pheatmap(FC_metab_medl2fc_plot,
              fontsize_row = 10,
              color=nc,
              treeheight_row = 0,
              show_rownames = TRUE,
              cluster_rows = FALSE,
              clustering_distance_rows = "correlation",
              
              cluster_cols = FALSE,
              show_colnames = TRUE,
              
              breaks = breaks)
pdf("~/hayflick/TRDFINALmetabs_NNMT.pdf", height=5,width=6)
print(plot)
dev.off()




##figure 3
panel1=c("Uridine diphosphate glucose",
         "Uridine diphosphate-N-acetylglucosamine","FRUCTOSE 6-PHOSPHATE",
         "GLYCERALDEHDYE-3-PHOSPHATE","Glycerol 3-phosphate",
         "L-Lactic acid",
         "ALPHA-D-GLUCOSE 1-PHOSPHATE","2-Phosphoglyceric acid", "Pyruvic acid",
         "RIBOSE 5-PHOSPHATE",
         "CARNITINE","DEOXYCARNITINE","L-Acetylcarnitine",
         "LysoPC(16:0)","LysoPC(18:1(9Z))","LysoPE(16:0/0:0)","LysoPE(18:0/0:0)",
         "CITRIC ACID","ISOCITRIC ACID","Oxoglutaric acid","L-Malic acid","Fumaric acid")

##Figure 4
panel1=c("1-Methylnicotinamide",
         "L-Methionine",
         "S-Adenosylhomocysteine",
         "S-Adenosylmethionine",
         "NAD",
         "NAD+",
         "NICOTINAMIDE")

##Kennedy supp

panel1=c("CDP-Ethanolamine",
         "O-Phosphoethanolamine",
         "Glycerophosphocholine")


panel1=c("L-Proline",
         "trans-4-Hydroxy-L-proline",
         "4-Hydroxyproline",
         "L-Lysine","L-Alanine","L-Arginine",
         "L-Asparagine","L-Aspartic acid","L-Cysteine",
         "L-Glutamic acid","L-Glutamine","L-Histidine", "L-Isoleucine",
        "L-Leucine","L-Lysine","L-Histidine",
         "L-Methionine","L-Phenylalanine","L-Serine",
         "L-Threonine","L-Tryptophan","L-Tyrosine","L-Valine","L-Cysteine")


panel1="Oxoglutaric acid"


###############################






###############################




writeLines(capture.output(sessionInfo()), "~/hayflick/figures/Figure_3/sessionInfo.txt")




