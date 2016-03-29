#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character", 
              help="A comma separated table, describing samples, read count file and conditions"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="The output read count matrix file (CSV format)")
)

parser <- OptionParser(usage="%prog -i sampleTable.csv -o countMatrix.csv ", option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

if(length(arguments$input)>0){
    

    sampleTable <- read.csv(file=arguments$input)

    suppressPackageStartupMessages(library(DESeq2))
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~condition)

    raw.counts <- as.data.frame(counts(ddsHTSeq))

    write.csv(raw.counts, file = arguments$output)
    
    writeLines(capture.output(sessionInfo()), 
               file.path(dirname(arguments$output), "make_count_matrix.sessionInfo"))
    
}else{
    cat("Please give a valid input table file\n")
    print_help(parser)
    stop()
}
