#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character", 
              help="A comma separated list of labels, one for each filteredJunctions.bed output given as input"),
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated list of filteredJunctions.bed outputs. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output file (CSV format)")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1/filteredJunctions.bed,S2/filteredJunctions.bed,S3/filteredJunctions.bed -o cfinder_compared.csv ", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)
opt <- arguments$options
minreads <- opt$minreads

if(length(arguments$input)>0){
  
  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file = x, sep="\t", header = F, stringsAsFactors = F); a}
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
