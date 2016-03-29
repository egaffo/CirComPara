#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character", 
              help="A comma separated list of labels, one for each CIRCexplorer_circ.txt file given as input"),
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated list of CIRCexplorer_circ.txt files. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output file (CSV format)")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1/CIRCexplorer_circ.txt,S2/CIRCexplorer_circ.txt,S3/CIRCexplorer_circ.txt -o CIRCexplorer_compared.csv ", 
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)
opt <- arguments$options

if(length(arguments$input)>0){
  
  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  colnms <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
              "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType",
              "geneName", "isoformName", "exonIndex_intronIndex", "flankIntron")
  
  colcls <- c("factor", "integer", "integer", "factor", "numeric", "factor", "integer", "integer",
              "character", "integer", "character", "character", "integer", "factor",
              "factor", "factor", "character", "character")
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file=x, sep="\t", header=F, stringsAsFactors = F, col.names=colnms, colClasses=colcls); a}
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
  