#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

cumdiff <- function(x){
  y <- as.vector(c(0))
  for(i in seq_along(x)){
    y[i] <- (x[i] - x[i+1])
  }
  end <- length(x)
  y[end] <- x[end]
  y
}

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", default="fastx_quality_stats.tsv",
              help="The FASTX-Toolkit fastx_quality_stats output table"),
  make_option(c("-o", "--output"), action="store", type="character", default="read_lengths.svg",
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

  title <- fastx_stats
  
  image_width <- 12
  image_height <- 7

  plot_template <- theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12), 
                         axis.title.x=element_text(size=16), axis.title.y=element_text(size=16, angle=90))
  
  fastx <- read.table(file=fastx_stats, header=TRUE, sep="\t")
  fastx$cycle <- ordered(fastx$cycle)
  fastx$length_hist <- cumdiff(fastx$ALL_count)
  
  lengths_plot <- ggplot(aes(x=cycle, y=length_hist), data=fastx, main=title) +  
			xlab("Read lengths") + ylab("Read count") + 
			geom_bar(stat="identity", aes(fill=I("blue"))) + 
                  theme_bw() + plot_template + scale_y_continuous(labels=comma, expand=c(0.01,0)) + 
                  scale_x_discrete(breaks=seq(min(as.integer(fastx$cycle)), max(as.integer(fastx$cycle)), 5))
  ggsave(filename=output, plot=lengths_plot, width=image_width, height=image_height)  
}

