##################
# help functions #
##################
get.seqs <- function (org, regions, no.cores = 1) {
  get.seq.chr <- function(chr) {
    seq <- org[[chr]]
    if (class(seq) == "MaskedDNAString") 
      seq <- unmasked(seq)
    chr.regions <- regions[which(as.character(seqnames(regions)) == chr)]
    region.seq <- DNAStringSet(Views(seq, start = start(chr.regions), 
                                     end = end(chr.regions)))
    rev <- as.logical(strand(chr.regions) == "-")
    region.seq[rev] <- reverseComplement(region.seq[rev])
    names(region.seq) <- values(chr.regions)$.ZZ
    gc()
    return(region.seq)
  }
  values(regions)$.ZZ <- sprintf("Unit%d", 1:length(regions))
  seqs <- unique(as.character(seqnames(regions)))
  all.seqs <- mclapply(seqs, get.seq.chr, mc.cores = no.cores)
  all.seqs <- do.call(c, all.seqs)
  inds <- sort(as.numeric(gsub("Unit", "", names(all.seqs))), 
               index.return = TRUE)$ix
  all.seqs <- all.seqs[inds]
  all.seqs <- as.character(all.seqs)
  gc()
  all.seqs <- DNAStringSet(all.seqs)
  gc()
  return(all.seqs)
}

library(ggplot2)
draw_diamond <- function(rna.data, atac.data, genes=NULL, top.n=NULL, max.gene.dist=50000, diamond.dist=0.1) {
  # inspired by Yuri's code for drawing
  # @rna.data requires column: gene, log2FC, padj
  # @atac.data requires column: gene, gene.dist, log2FC, padj
  # @either genes or top.n need to be specified
  
  # data preparation
  rna.data$gene <- as.character(rna.data$gene)
  atac.data$gene <- as.character(atac.data$gene)

  rna.data <- rna.data[rna.data$gene %in% atac.data$gene,] # only consider genes with atac peak
  rownames(rna.data) <- rna.data$gene
  atac.data <- atac.data[atac.data$gene.dist <= max.gene.dist, ]
  
  # plot1 is RNA-seq for selected genes
  if (!is.null(genes)) {
    toplot1 <- rna.data[genes, c("gene","log2FC","padj")]
    toplot1 <- toplot1[order(toplot1$log2FC), ]
    toplot1$rank <- 1:nrow(toplot1)
  }
  if (!is.null(top.n)) {
    toplot1 <- rna.data[, c("gene","log2FC","padj")]
    toplot1 <- toplot1[order(toplot1$log2FC), ]
    toplot1 <- toplot1[c(1:top.n, (nrow(toplot1)-top.n+1):nrow(toplot1)), ]
    toplot1$rank <- 1:nrow(toplot1)
  }
  
  # plot2 is ATAC
  toplot2 <- atac.data[atac.data$gene %in% toplot1$gene, c("gene", "log2FC", "padj")]
  toplot2$x <- toplot1[toplot2$gene,"rank"]
  toplot2$y <- 0
  toplot2 <- split(toplot2, toplot2$gene)
  toplot2 <- lapply(toplot2, function(tmp) {
    tmp <- tmp[order(tmp$log2FC), ]
    y <- rep(toplot1[unique(tmp$gene), "log2FC"], nrow(tmp))
    y <- y + seq(diamond.dist, diamond.dist * nrow(tmp), by = diamond.dist)
    tmp$y <- y
    return(tmp)
  })
  toplot2 <- do.call(rbind, toplot2)
  
  # set maximum and min value of atac log2FC to be 75% quantile
  max_val <- quantile(toplot2$log2FC[toplot2$log2FC>0], 0.75)
  min_val <- quantile(toplot2$log2FC[toplot2$log2FC<0], 0.25)
  toplot2$log2FC[toplot2$log2FC>max_val] <- max_val
  toplot2$log2FC[toplot2$log2FC<min_val] <- min_val


  # plot all
  y_max <- max(abs(toplot1$log2FC))+1
  p <- ggplot(toplot1,  aes(rank, log2FC, label = gene)) + 
    geom_text(aes(angle=90), hjust = 1) +
    geom_point(data=toplot2, aes(x, y, color=log2FC)) +
    geom_point(data=toplot2, aes(x, y), shape=1, size=2, stroke=0.2) +
    labs(colour = "ATAC_log2FC") +
    scale_colour_gradient2(low="royalblue3", mid="white", high="firebrick3", midpoint = 0) +
    theme_classic() +  ylim(-y_max,y_max) + xlab("") + ylab("RNA_log2FC") +
    geom_hline(yintercept=0, linetype = 2)
  return(p)
  
}


