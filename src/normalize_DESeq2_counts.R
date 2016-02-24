#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated table of raw read counts for different samples as given by DESeqDataSetFromHTSeqCount()"),
  make_option(c("-s", "--samples"), action="store", type="character", 
              help="A comma separated table, describing samples and conditions"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The normalized read count matrix file (CSV format)")
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
  
}else{
  cat("Please give a valid input table file\n")
  print_help(parser)
  stop()
}
