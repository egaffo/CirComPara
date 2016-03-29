#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated table of raw read counts for different samples as given by DESeqDataSetFromHTSeqCount()"),
  make_option(c("-s", "--samples"), action="store", type="character", 
              help="A comma separated table, describing samples and conditions"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The normalized read count matrix file (CSV format)"),
  make_option(c("-d", "--plotDirectory"), action="store", type="character", default = ".",
              help="The directory where to store plots")
)

parser <- OptionParser(usage="%prog -i rawReads.csv -o normalizedReads.csv -s sampleTable.csv", option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

if(length(arguments$input)>0){
  
  sampleTable <- read.csv(file = arguments$sample)
  rownames(sampleTable) <- as.character(sampleTable[, 1])
  countData   <- read.csv(file = arguments$input, row.names = 1)
  
  suppressPackageStartupMessages(library(DESeq2))
  
  ddsHTSeq <- DESeqDataSetFromMatrix(countData = countData,
                                colData = sampleTable,
                                design = ~ condition)
  dds <- DESeq(ddsHTSeq)
  
  norm.counts <- as.data.frame(counts(dds, normalized = TRUE))
  write.csv(norm.counts, file = arguments$output)
  
  ## PCA plot
  rld <- rlog(dds, blind = T)
  pca <- prcomp(t(assay(rld)))
  
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggthemes))
  suppressPackageStartupMessages(library(RSvgDevice))
  
  intgroup <- c("condition")
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  ## compute PC variance explained % and plot
  percentVar.df <- data.frame(PCs = paste0("PC", 1:length(percentVar)), 
                              Variation_percentage = percentVar*100)
  pc.var.explained.plot <- ggplot(data = percentVar.df, aes(x = PCs, y = Variation_percentage)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = round(Variation_percentage, 2)), vjust = -.5) + 
    ylab("Variation explained %") + theme_bw()
  
  print("Saving principal components variation explained plot ...")
  devSVG(file = file.path(arguments$plotDirectory, "pc_variation_plot.svg"), 
         width = 10, height = 10, 
         bg = "white", fg = "black", onefile = TRUE, xmlHeader = TRUE)
  print(pc.var.explained.plot)
  dev.off()
  
  ## compute and plot PCs
  intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  } else {
    colData(rld)[[intgroup]]
  }
  d <- data.frame(pca$x, group = group, intgroup.df, name = colnames(rld))
  pca.plot <- ggplot(data = d, aes(PC1, PC2, color=condition, shape=condition, label=name)) + 
    geom_point(size = 3) + coord_fixed() + 
    geom_text(position="jitter", hjust=0.2, vjust=-0.2, size=6, show_guide = F) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+ 
    scale_color_brewer(palette = "Dark2") + scale_x_continuous(expand = c(0.3,0)) +
    scale_y_continuous(expand = c(0.3,0)) + theme_bw() + 
    theme(legend.position="bottom")
  
  print("Saving principal components analysis plot ...")
  devSVG(file = file.path(arguments$plotDirectory, "pca_plot.svg"), 
         width = 10, height = 10, 
         bg = "white", fg = "black", 
         onefile = TRUE, xmlHeader = TRUE)
  print(pca.plot)
  dev.off()
  
  writeLines(capture.output(sessionInfo()), 
             file.path(dirname(arguments$output), "normalize_DESeq2_counts.sessionInfo"))
  
}else{
  cat("Please give a valid input table file\n")
  print_help(parser)
  stop()
}
