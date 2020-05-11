#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                help = "Findcirc result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "integer",
                default = 40,
                help = "Minimum quality of backsplice end mapping"),
    make_option(c("-f", "--filter_tags"), action = "store", type = "character",
                default = 'CIRCULAR',
                help = "Comma separated list of flags that have to be present in rows of the Find_circ BED file"),
    make_option(c("-o", "--output"), action = "store", type = "character",
                default = "filteredFindcirc.bed", help = "The output file")
)

parser <- OptionParser(usage = "%prog -i circ_candidates.bed -q 40 -f CIRCULAR -o filteredFindcirc.bed",
                       option_list = option_list,
                       description = "Filter Findcirc results by quality")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
output <- arguments$output
minqual <- arguments$minqual
flags <- unlist(strsplit(arguments$filter_tags, ','))
if(length(flags) < 1){
    flags <- '.'
}

circ_candidates.bed <-
    fread(cmd = paste0(c(paste0("grep -v '^#' ", input),
                         flags),
                       collapse = " | grep -w "),
          showProgress = F)

if(nrow(circ_candidates.bed) > 0){
    circ_candidates.bed <- circ_candidates.bed[V9 >= minqual & V10 >= minqual]
}

write.table(x = circ_candidates.bed,
            file = output, row.names = F, quote = F, sep = "\t", col.names = F)
