BiocManager::install("ATACseqQC")


BiocManager::install(c( "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg38"))

library(ATACseqQC)
library(txs)
library(GenomicAlignments)



set()
  bams=list.files(path="/home/yuanh/analysis/WI38_ATAC/process/bam",pattern="filter.rmdup.bam",full.names = TRUE)
bams=bams[!grepl("bai",bams)]
bams=bams[grepl("omni",bams)]
bams=bams[-39]
bams=bams[-36]
bams=bams[-18]
bams=bams[-15]





######################## TSSE

params = ScanBamParam(mapqFilter = 1, 
                      flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
                                                                                         "mapq", "isize"), which = GRanges("chr1", IRanges(1, 248956422)))
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)



plotty=list()
for (i in seq(1:length(bams))){
  #gal <- readBamFile(bams[i], tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  #gal1 <- shiftGAlignmentsList(gal)
  gal1 <- readGAlignments(bams[i],param = params)
  tsse <- TSSEscore(gal1, txs)
  mper[i,6]=tsse$TSSEscore  
  
 # pdfname=paste0(mper[i,2],"_",mper[i,3],"_",mper[i,4],".pdf")
plotty[[i]]=
  print(ggplot(as.data.frame(tsse),aes(x=100*(-9:10-.5),y=values)) +
   geom_smooth(method="gam",se = FALSE) +
   theme_bw() +
   ylim(c(0,70)) +
   xlab(label = "distance from TSS") +
   ylab(label = "aggregate TSS score"))
}  

library(gridExtra)

pdf("plots_TSSe.pdf", onefile = TRUE)
#for (i in seq(length(plotty))) {
  do.call("grid.arrange", plotty[[i]])  
#}
dev.off()

  pdf(pdfname, height=4,width=4)
    dev.off()


mper$strain="hTERT"
mper[19:39,7]="WT"
mper=mper[-18,]
mper=mper[-15,]
mper=mper[-39,]




######################################### FRAG LENGTHS


library(magrittr)
library(dplyr)
library(ggplot2)

for (i in seq(1:length(bams))){

atacReads <- readGAlignmentPairs(bams[i], param = ScanBamParam(mapqFilter = 1, 
                                                                flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
                                                                                                                           "mapq", "isize"), which = GRanges("chr1" ,IRanges(1, 10025520))))
temp=bams[i]
temp=str_split_fixed(temp,"/",8)
temp=temp[,8]
temp=str_split_fixed(temp,"_",4)
pdfname=paste0("~/hayflick/",temp[,1],"_",temp[,2],"_fragLen.pdf")
 
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)

fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                 Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

pdf(pdfname)
print(fragLenPlot + theme_bw() + xlim(25, 250))
dev.off()
}




################# MITO PERCENT

library(Rsamtools)
library(ggplot2)
library(magrittr)

sortedBAM=(bams[40])
mper=as.data.frame(bams)
temp=str_split_fixed(mper$bams,"_",6)
temp2=str_split_fixed(temp[,2],"/",4)
mper$line=temp2[,4]
mper$tp=temp[,3]
mper$rep=temp[,5]
mper$mper=0
for (i in seq(1:length(bams))){
  idxstat=idxstatsBam(bams[i])
sum=sum(idxstat$mapped)
mnum=idxstat[25,3]/sum *100
mper[i,5]=mnum


}

################## PLOTS


mper$strain="hTERT"
mper[17:35,7]="WT"

library(yarrr)


colnames(mper)

pdf("/home/dgh/multiqc/hayflick_atac/mito_perecent.pdf", height=6,width=7)
yarrr::pirateplot(formula = mper ~tp + strain ,    # DV = height, IV1 = sex, IV2 = headband
                  data = mper, 
                  ylim=c(0,5),         
                  theme =4, inf.f.o = 0, avg.line.o = 0,point.o = 1,
                  main = "Pirate Heights",
                  pal = "gray")

dev.off()

pdf("/home/dgh/multiqc/hayflick_atac/tp_perecent.pdf", height=6,width=7)
yarrr::pirateplot(formula = V6 ~tp + strain ,    # DV = height, IV1 = sex, IV2 = headband
                  data = mper, 
                  ylim=c(0,110),         
                  theme =4, inf.f.o = 0, avg.line.o = 0,point.o = 1,
                  main = "Pirate Heights",
                  pal = "gray")

dev.off()




#################### FRIP SCore

frip=read.csv("~/hayflick/data/ATACseq/FRIP_results.txt", stringsAsFactors = FALSE)

frip=frip[-39,]
frip=frip[-36,]
frip=frip[-18,]
frip=frip[-15,]

plot=ggplot(frip,aes(x=timepoint , y=reads_in_peaks),color=treatment) +
  scale_x_discrete() +
  geom_jitter(aes( x = timepoint), 
              position = position_jitter(width = .05), alpha = 1,size=1) +
  facet_grid(~ treatment,space = "free",scales = "free_x") +
  scale_color_manual(values = paper_cols) +
  scale_fill_manual(values = paper_cols) +
  # geom_boxplot(aes(x = numb,color=condition, fill=condition), outlier.colour = NA, position = dodge,alpha=.5) +
  geom_bar(aes(fill=treatment,),stat="summary",position = dodge,alpha=.6)+
  theme_bw() + 
  theme(legend.position = "none")


frip_pdl=frip %>% filter(treatment== "PDL")

frip_pdl=frip_pdl %>% dplyr::mutate(pdl_n = c(rep(.14,3),rep(.29,3),rep(.43,3),rep(.57,3),rep(.71,3),rep(.85,2),rep(1,2)))




t=lm(reads_in_peaks ~ pdl_n, data = frip_pdl)
t.t=tidy(t)



frip_tert=frip %>% filter(treatment== "hTERT")
frip_tert=frip_tert %>% dplyr::mutate(pdl_n = c(rep(.14,3),rep(.43,3),rep(.57,3),rep(.71,3),rep(.85,2),rep(1,2)))
s=lm(reads_in_peaks ~ pdl_n, data = frip_tert)
s.t=tidy(t)

