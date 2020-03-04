#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character", 
              help="A comma separated list of labels, one for each CircRNACount output given as input"),
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated list of CircRNACount outputs. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output file (CSV format)")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1/CircRNACount,S2/CircRNACount,S3/CircRNACount -o dcc_compared.csv ", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)
opt <- arguments$options
minreads <- opt$minreads

if(length(arguments$input)>0){
  
  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file = x, sep="\t", header = T, stringsAsFactors = F); a}
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
