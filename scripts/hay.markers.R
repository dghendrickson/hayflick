

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






############################ PROTEI

pdata=read.delim("~/hayflick/data/proteomics/Proteomics_raw_proportions_RS_and_tert.txt", stringsAsFactors = FALSE)
pdata=pdata[complete.cases(pdata),]

pdata_pdl=pdata %>% select(contains("PDL"))

tp1=rowMeans(pdata_pdl[,1:3])

pdata_pdl=sweep(pdata_pdl,1,tp1,"/")


pdata_tert=pdata %>% select(contains("TERT"))


