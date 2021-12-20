library(org.Hs.eg.db)
library(rtracklayer)

read_narrowPeak <- function(narrowPeak_file) {
  extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  gr_narrowPeak <- import(narrowPeak_file, format = "BED",
                          extraCols = extraCols_narrowPeak)
  return(gr_narrowPeak)
}

###########################
# create annotation table #
###########################
makeGTF1 <- function(annotation_gtf, org=org.Hs.eg.db){
  annotationTable<-fread(annotation_gtf)
  temp<-gsub("gene_id \"","",annotationTable$V9)
  annotationTable$gene_id<-gsub("\".*","",temp)
  annotationTable$symbol<-NA
  temp<-gsub(".*transcript_id \"","",annotationTable$V9)
  annotationTable$transcript_id<-gsub("\";","",temp)
  map <- select(org, keys=annotationTable$gene_id, columns=c("REFSEQ","SYMBOL"), keytype="REFSEQ")
  map <- map[!duplicated(map[,1]),]
  rownames(map) <- map[,1]
  annotationTable$symbol <- map[annotationTable$gene_id,2] 
  a <- annotationTable[,list(seqname=V1,annot=V3,start=V4,end=V5,strand=V7,symbol=symbol,transcript_id=transcript_id)]
  a$width<-abs(a$end-a$start)+1
  a<-a[a$seqname%in%c(paste0("chr", 1:22), "chrX", "chrY","chrM"),]
  return(a)
}

###########################
# create annotation table #
###########################
annotateSites <- function(sitesGR, annotationTab, exon_ovlp_ratio_threshold=0.5) {  
  
  # make sure peaks have names
  if (is.null(names(sitesGR))) {
    stop("peak names must be specified!")
  }
  
  # change it to a transcript table, where all rows of the same transcript becomes the same, with the same start and end site
  #transcripts <- annotationTab[, list(start = min(start), end = max(end), symbol = symbol), by = list(transcript_id, seqname, strand)]
  transcripts <- unique(annotationTab[, list(start = min(start), end = max(end), symbol = symbol), by = list(transcript_id, seqname, strand)])
  
  #################################
  # create a transcript GR object #
  #################################
  transcriptsGR <- GRanges(seqnames = transcripts$seqname,
                           ranges = IRanges(start = transcripts$start,end = transcripts$end),
                           strand = transcripts$strand,
                           transcript_id = transcripts$transcript_id,
                           symbol = transcripts$symbol
  )
  names(transcriptsGR) <- transcriptsGR$transcript_id
  
  ##################
  # Promoter sites #
  ##################
  tss <- rbind(transcripts[strand == "+", list(transcript_id, seqname, strand, pos = start, symbol)],
               transcripts[strand == "-", list(transcript_id, seqname, strand, pos = end, symbol)])
  tssGR <- GRanges(seqnames = tss$seqname,
                   ranges = IRanges(start = tss$pos,end = tss$pos),
                   strand = tss$strand,
                   transcript_id = tss$transcript_id,
                   symbol = tss$symbol)
  names(tssGR) <- tssGR$transcript_id
  strand(tssGR) <- "*" # ignore the strand of tssGR
  
  #one index for each row of sitesGR, shows where is its closest tss
  nearest_id <- nearest(x = sitesGR, subject = tssGR)
  
  output_table <- data.table(peak_id = names(sitesGR),
                             transcript_id = tssGR[nearest_id]$transcript_id,
                             symbol = tssGR[nearest_id]$symbol,
                             annot="promoter",
                             dist = distance(x = sitesGR,y = tssGR[nearest_id]))
  output_table <- output_table[dist <= 2000]
  
  ###################
  # Intergenic sites.
  ###################
  # sites that are not in matchTranscript2SiteTab(promoter), and not in transcripts
  strand(transcriptsGR) <- "*" # ignore strand of transcriptsGR from this point.
  overlapWithTranscr <- subsetByOverlaps(sitesGR[!names(sitesGR) %in% output_table$peak_id], transcriptsGR) # anything that ovlps with transcriptGR
  intergPks <- sitesGR[!names(sitesGR) %in% union(output_table$peak_id, names(overlapWithTranscr))] # anything that not promoter or overlapWithTranscr
  
  # compute distance
  nearest_id <- nearest(x = intergPks, subject = transcriptsGR)
  
  # add intergenic sites to the table
  intergenicSites <- data.table(peak_id = names(intergPks),
                                transcript_id = transcriptsGR[nearest_id]$transcript_id,
                                symbol = transcriptsGR[nearest_id]$symbol,
                                annot = "intergenic",
                                dist = distance(x = intergPks, y = transcriptsGR[nearest_id]))
  
  output_table <- rbind(output_table, intergenicSites)
  
  ##############
  # Exon sites #
  ##############
  exonsGR <- GRanges(seqnames = annotationTab[annot == "exon"]$seqname,
                     ranges = IRanges(start = annotationTab[annot == "exon"]$start,
                                      end = annotationTab[annot == "exon"]$end),
                     strand = "*",
                     transcript_id = annotationTab[annot == "exon"]$transcript_id,
                     symbol = annotationTab[annot == "exon"]$symbol
  )
  names(exonsGR) <- 1:length(exonsGR)
  ovlWithExon <- findOverlaps(overlapWithTranscr,exonsGR)
  ovlWithExonTab <- data.table(peak_id = names(overlapWithTranscr[queryHits(ovlWithExon)]),
                               exon_id = names(exonsGR[subjectHits(ovlWithExon)]))
  intersectRanges <- pintersect(overlapWithTranscr[ovlWithExonTab$peak_id],
                                exonsGR[ovlWithExonTab$exon_id])
  ovlWithExonTab$intersectWidth <- width(intersectRanges)
  ovlWithExonTab$peak_width <- width(overlapWithTranscr[ovlWithExonTab$peak_id])
  ovlWithExonTab$intersectRatio <- ovlWithExonTab$intersectWidth/ovlWithExonTab$peak_width
  
  # filtered exon pks table
  exonPksTab <- ovlWithExonTab[, list(exon_id = exon_id[which.max(intersectRatio)],intersectRatio = max(intersectRatio)), by = peak_id][intersectRatio >= exon_ovlp_ratio_threshold]
  exonPksTab$transcript_id <- exonsGR[exonPksTab$exon_id]$transcript_id
  exonPksTab$symbol <- exonsGR[exonPksTab$exon_id]$symbol
  
  # exon pk annotations
  exonSites <- data.table(peak_id = exonPksTab$peak_id,
                          transcript_id = exonPksTab$transcript_id,
                          symbol = exonPksTab$symbol,
                          annot = "exon",
                          dist=0)
  
  output_table <- rbind(output_table, exonSites)
  
  ################
  # Intron sites #
  ################
  intronPks <- sitesGR[setdiff(names(sitesGR), output_table$peak_id)]
  nearest_id <- nearest(x = intronPks, subject = transcriptsGR)
  intronSites <-  data.table(peak_id = names(intronPks),
                             transcript_id = transcriptsGR[nearest_id]$transcript_id, 
                             symbol = transcriptsGR[nearest_id]$symbol, 
                             annot = "intron",
                             dist = 0)
  
  output_table <- rbind(output_table, intronSites)
  
  #######################
  # reorder ouput table #
  #######################
  setkey(output_table, peak_id)
  output_table <- output_table[names(sitesGR)] # reorder match
  
  return(output_table)
}

add_annot1 <- function(atlas, annotation_gtf="/home/yuanh/programs/genomes/hg38/hg38_refGene.gtf", org=org.Hs.eg.db) {
  # use REFSEQ annotation here.
  # only consider annotated genes.
  annot <- makeGTF1(annotation_gtf, org)
  annot <- annot[!is.na(symbol),]
  annot <- annot[symbol!=""]
  
  match <- annotateSites(atlas, annot) # match is not the same order as atlas
  atlas$annot <- match$annot
  atlas$transcript_id <- match$transcript_id
  atlas$gene_name <- match$symbol
  atlas$distToNearest <- match$dist
  
  return(atlas)
}

###############
# merge peaks #
###############
collapse.pairwise.celltype.peaks <- function (peaks1, peaks2, overlap.ratio.cutoff=0.75) {
  ## Function to find overlapping ratios
  find.overlap.ratio <- function (gr1, gr2) {
    ## Find Overlaps
    overlaps <- as.matrix (findOverlaps (gr1, gr2))
    if (nrow(overlaps)==0){
      output = NULL
    } else {
      ## Build ranges
      ranges <- cbind (start (gr1)[overlaps[,1]], end (gr1)[overlaps[,1]],
                       start (gr2)[overlaps[,2]], end (gr2)[overlaps[,2]])
      ranges <- t(apply (ranges, 1, sort))
      ## Min widths
      widths <- pmin (width (gr1)[overlaps[,1]], width (gr2)[overlaps[,2]])
      ## Overlap ratios
      overlap.ratio <- (ranges[,3] - ranges[,2])/widths
      #for querys mapping to the same subject, find the rows in overlaps that represent the query with best mapping
      best.gr1 <- tapply (1:nrow (overlaps), overlaps[,1], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
      #for subjects mapping to the same query, find the rows in overlaps that represent the subject with best mapping
      best.gr2 <- tapply (1:nrow (overlaps), overlaps[,2], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
      common <- intersect (best.gr1, best.gr2)      
      output <- list (overlap.ratio = overlap.ratio,
                      ranges = ranges,
                      overlaps = overlaps,
                      best.common = common)
    } 
    return (output)
  }
  ## Reset metadata
  values (peaks1) <- values (peaks2) <- NULL
  ## Overlapping peaks which exceed overlap.ratio.cutoff
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common[or$overlap.ratio[or$best.common] > overlap.ratio.cutoff]
  union <- GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                    IRanges (or$ranges[common,1], or$ranges[common,4])) # modified from Manu's code. He's taking the intersection, I am taking the union.
  ## Overlapping peaks with ratio < overlap.ratio.cutoff
  peaks1 <- peaks1[countOverlaps (peaks1, union) == 0]
  peaks2 <- peaks2[countOverlaps (peaks2, union) == 0]    
  or <- find.overlap.ratio (peaks1, peaks2)
  if (is.null(or)){
    common <- or$best.common
    union <- c(union,GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                              IRanges (or$ranges[common,1], or$ranges[common,4])))    
  }
  ## Non overlapping peaks
  union <- c(union, peaks1[countOverlaps (peaks1, union) == 0])
  union <- c(union, peaks2[countOverlaps (peaks2, union) == 0])
  return (union)
}
mergeCTpeaks <- function(pkLst) {
  combined.peaks <- collapse.pairwise.celltype.peaks (pkLst[[1]], pkLst[[2]])
  if (length(pkLst)>2){
    for (i in 3:length(pkLst)){
      combined.peaks <- collapse.pairwise.celltype.peaks(combined.peaks,pkLst[[i]])
      print(i)
    }
  }
  ## annotate each sample where peak was called.
  for (sample in names(pkLst)) {    
    overlaps <- countOverlaps (combined.peaks, pkLst[[sample]])
    mcols (combined.peaks)[,paste0(sample, ".peak")] <- 0
    mcols (combined.peaks)[,paste0(sample, ".peak")][overlaps > 0] <- 1
  }
  names(combined.peaks) <- 1:length(combined.peaks)
  combined.peaks$pattern <- apply(as.matrix(mcols(combined.peaks)), 1, function(row) {
    argsToPaste0 <- as.list(row)
    do.call(paste0, argsToPaste0)
  })
  return(combined.peaks)
}

###############
# count reads #
###############
countReads <- function(atlas, bam_file_list, sample_names) {
  # bam_file_list must be rds files containing shifted reads.
  output <- list()
  for (i in 1:length(bam_file_list)) {
    tmp <- readRDS(bam_file_list[[i]])
    output[[i]] <- countOverlaps(atlas, tmp)
    print(i)
  }
  output <- do.call(cbind, output)
  colnames(output) <- sample_names
  return(output)
}

######################
# plotting functions #
######################
# plot pie chart
piechart<-function(v,outname){
  freq<-table(v)
  pdf(outname)
  pie(freq, 
      labels=paste0(names(freq),":",freq), 
      main=paste0("pie chart"))
  dev.off()  
}

###########
# KS test #
###########
motif_KS <- function (log2FC, outdir, compare, m) {
  m_list <- lapply(1:ncol(m), function(i) {
    bound <- log2FC [m[,i]==1]
    return(bound)
  })
  names(m_list) <- colnames(m)
  
  # compute % motif hits and KS test statistic
  result <- matrix(0, length(m_list), 4, dimnames=list(names(m_list), c("direction","z","effective_size","%hits")))
  # use log2FC on the atlas as background
  for (i in 1:length(m_list)) {
    bound <- m_list[[i]]
    n1 <- as.numeric(length(bound))
    n2 <- as.numeric(length(log2FC))
    right <- ks.test(bound, log2FC, alternative = "less")$statistic
    left <- ks.test(bound, log2FC, alternative = "greater")$statistic
    
    result[i, 1] <- which.max(c(right, left))
    result[i, 2] <- max(c(right, left))
    result[i, 3] <- result[i, 2] / sqrt(n1*n2/(n1+n2))
    result[i, 4] <- round( length(bound) / nrow(m) * 100, 1)
    
    if (i%%100==0) print(i) # print for every 100 TFs
  }
  
  # plot
  toplot <- data.frame(result)
  toplot$x <- result[,4]
  toplot$y <- ifelse(result[,1]==1, 1, -1) * result[,2]
  toplot$TF <- rownames(toplot)
  toplot$color <- "gray"
  toplot$color[order(toplot$y)[1:10]] <- "blue"
  toplot$color[order(toplot$y,decreasing=T)[1:10]] <- "red"
  toplot$color <- factor(toplot$color, levels = c("blue","gray","red"))
  fwrite(toplot, file=paste0(outdir, "/fimo_KS_",compare,".csv"))
  
  pdf(paste0(outdir, "/fimo_KS_",compare,".pdf"))
  p <- ggplot(toplot, aes(x, y)) + geom_point(aes(x,y,color=color)) +
    xlim(c(0, 50)) + ylim(c(-0.4, 0.4)) + xlab("% motif hits") + ylab("KS statistic") + ggtitle(compare) + theme_classic() +
    #geom_point(data = toplot[toplot$color!="gray", ], aes(x, y, color=color)) +
    scale_colour_manual(values = c("blue","gray", "red")) +
    geom_text_repel(data = toplot[toplot$color!="gray", ], aes(x, y, label = TF)) +
    geom_hline(yintercept=0, linetype = 2)
  print(p)
  dev.off()
  
  return(toplot) 
}



