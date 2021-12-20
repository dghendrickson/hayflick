library(Hmisc)
library(ggpubr)

################ recalculate bulk RNAseq vs. PDL 25 (single cell starts with PDL 25) so that we can compare PDL 50v25 for bulk and SC

RS=read.delim("~/hayflick/data/RNAseq/TPM_tables/RS_batch_tpms.txt", stringsAsFactors = FALSE) ##get tetrt info and avg. replicates

colnames(RS)

RS_pivot=RS %>% tidyr::pivot_longer(-gene) 

RS_pivot$value=RS_pivot$value +1 #add pseudocount

sample_tracking=str_split_fixed(RS_pivot$name,"_",3) #### break up sample name to isolate identifier common to replicate

RS_pivot$group=sample_tracking[,1] #### use identifier common to group to average

RS_pivot_avg=RS_pivot %>%                      #### average and recast wide
  group_by(group,gene) %>%
  dplyr::summarize(expression_mean = mean(value)) %>% 
  tidyr::pivot_wider(id_cols=gene,names_from=group, values_from=expression_mean)


RS_pivot_avg_fc=log2(RS_pivot_avg[,4:9]/RS_pivot_avg$PDL25 )  #### compute log2 FC vs PDL 25 for this specific analyis
RS_pivot_avg_fc=cbind(RS_pivot_avg$gene,RS_pivot_avg_fc)
colnames(RS_pivot_avg_fc)[1]="gene"

RS_PDL50v25=as.data.frame(RS_pivot_avg_fc$PDL50, row.names = RS_pivot_avg_fc$gene)


################################Pull in SC from PDL vs PDL 25 Single cell matrix used in FIg2 script


###### CREATE LOG2 FC vs PDL25 CountsPerMillion (cpm) for all cells agreggated 

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


SC_pdl50v25=as.data.frame(PDL_cpms$PDL_50, row.names = rownames(PDL_cpms))

##############combined object

combined_bulk_SC=merge(RS_PDL50v25, SC_pdl50v25, by="row.names")

##bring in average ex from deseq results fig1
RS_deseq=read.delim("~/hayflick/data/RNAseq/DESEQ_output/RS_deseq.txt", stringsAsFactors = FALSE)

baseMean=RS_deseq[,c(1:2)]
baseMean$log10=log10(baseMean$baseMean +1) #pseudocount



combined_bulk_SC=merge(combined_bulk_SC, baseMean, by.x="Row.names",by.y="gene")
colnames(combined_bulk_SC)=c("gene","bulkFC50","SC50","bm","log10bm")

combined_bulk_SC$QuartilesBM= cut(combined_bulk_SC$log10bm, breaks=c(quantile(combined_bulk_SC$log10bm, probs = seq(0, 1, by = 0.25))), 
      labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)


ggplot(combined_bulk_SC) +
  geom_point(aes(x=SC50,y=bulkFC50,color=log10bm),size =1,alpha = .2) + 
  scale_color_continuous_sequential("Reds 3",begin = 0,end = 1, rev = TRUE) +
  #scale_colour_gradient(low = "white", high = "red") +
   geom_smooth(aes(x=SC50,y=bulkFC50,group = QuartilesBM),method = "lm") +
   stat_cor(aes(x=SC50,y=bulkFC50,group = QuartilesBM), method = "pearson") +
  #geom_smooth( aes(x=big_merge$b.Confluence,y=big_merge$b.repSen)) +
  scale_x_continuous(limits=c(-6,6)) + 
  scale_y_continuous(limits=c(-6,6)) +
  theme_bw(base_size = 14)  +
  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x=element_text(colour="black",size = 14),
        axis.text.y=element_text(colour="black",size = 14),
        axis.title = element_text(colour="black",size = 20)) #+
# theme(legend.position = "none")












