#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character", 
              help="A comma separated list of labels, one for each find_circ.py output given as input"),
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated list of find_circ.py outputs. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output (CSV format)")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1/circ_candidates.bed,S2/circ_candidates.bed,S3/circ_candidates.bed -o findcirc_compared.csv ", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

if(length(arguments$input)>0){
  
  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  colnames <- c("chrom", "start", "end", "name", "n_reads", "strand", "n_uniq", "best_qual_A", "best_qual_B", 
                "spliced_at_begin", "spliced_at_end", "tissues", "tiss_counts", "edits", "anchor_overlap", 
                "breakpoints", "unknown")
  
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file=x, sep="\t", header=F, col.names = colnames); a}
  )
  colnames(combined.df)[1] <- "sampleID"
  
  if(length(arguments$output)>0){
    write.csv(combined.df, file = arguments$output, row.names = F)
  }else{
    write.csv(combined.df, file = stdout(), row.names = F)
  }
}else{
  cat("Please give a valid input table file\n")
  print_help(parser)
  stop()
}