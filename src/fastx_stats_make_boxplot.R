#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))   # needed for formatting y-axis labels to non-scientific type

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", default="fastx_quality_stats.tsv",
              help="The FASTX-Toolkit fastx_quality_stats output table"),
  make_option(c("-o", "--output"), action="store", type="character", default="fastx_quality_stats.svg",
              help="Plot output file name")
  )
parser <- OptionParser(usage = "%prog [options]", 
                       option_list=option_list,
                       description="Make boxplot for FASTQ statistics on quality per base")
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
  
  quals_plot <- ggplot(fastx, aes(x=cycle, lower=ALL_Q1, upper=ALL_Q3, middle=ALL_med, ymin=ALL_lW, ymax=ALL_rW, 
                                  group=cycle)) + 
    geom_boxplot(stat="identity", fill=I("grey"), alpha=0.5) +
    geom_line(stat="smooth",aes(y=ALL_med, group=1), se=F, colour=I("red"), lwd=2, alpha=0.5, method="loess") +
    xlab("Posiotion in read sequences (nt)") + ylab("Phred quality score") + 
    theme_bw() + plot_template + ggtitle(title) + 
    scale_y_continuous(breaks=min(fastx$ALL_lW):max(fastx$ALL_rW)) + 
    scale_x_discrete(breaks=seq(min(as.integer(fastx$cycle)), max(as.integer(fastx$cycle)), 5))
  ggsave(filename=output, plot=quals_plot, width=image_width, height=image_height)
  
}