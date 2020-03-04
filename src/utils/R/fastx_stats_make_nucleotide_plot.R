#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(reshape2))

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", default="fastx_quality_stats.tsv",
              help="The FASTX-Toolkit fastx_quality_stats output table"),
  make_option(c("-o", "--output"), action="store", type="character", default="fastx_nucleotides.svg",
              help="Plot output file name")
)
parser <- OptionParser(usage = "%prog [options]", 
                       option_list=option_list,
                       description="Make nucleotides per position barplot")
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
fastx_stats <- opt$input

if(file.access(fastx_stats) == -1){
  stop(sprintf("Specified file ( %s ) does not exist", fastx_stats))
} else {
  
  output <- opt$output
  
  image_width <- 12
  image_height <- 7
  
  title <- fastx_stats
  
  plot_template <- theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12), 
                         axis.title.x=element_text(size=16), axis.title.y=element_text(size=16, angle=90))
  
  fastx <- read.table(file=fastx_stats, header=TRUE, sep="\t")
  fastx$cycle <- ordered(fastx$cycle)
  
  base_count_headers <- c("A_count", "C_count", "T_count", "G_count", "N_count")
  base_counts <- fastx[, c("cycle", base_count_headers)]
  m.base_counts <- melt(base_counts, id.vars = "cycle", value.name = "amount", variable.name = "nucleotide")
  
  nucleotide_plot <- ggplot(m.base_counts, aes(x=cycle, y=amount, fill=nucleotide)) + 
    geom_bar(position="fill", stat="identity") + scale_y_continuous(labels=percent_format()) + 
    xlab("Posiotion in read sequences (nt)") + theme_bw() + plot_template + ggtitle(title) + 
    scale_x_discrete(breaks=seq(min(as.integer(fastx$cycle)), max(as.integer(fastx$cycle)), 5))
  ggsave(filename=output, plot=nucleotide_plot, width=image_width, height=image_height)
}