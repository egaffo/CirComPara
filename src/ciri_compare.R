#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character", 
              help="A comma separated list of labels, one for each CIRI output given as input"),
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated list of CIRI outputs. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output  (CSV format)")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1_ciri.out,S2_ciri.out,S3_ciri.out -o ciri_compared.csv ", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

if(length(arguments$input)>0){
  
  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  #colnames <- c("circRNA_ID", "chr", "circRNA_start", "circRNA_end", "#junction_reads", "SM_MS_SMS", 
  #              "#non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "junction_reads_ID")
  
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file=x, sep="\t", header=T, comment.char="", stringsAsFactors = F); a}
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
