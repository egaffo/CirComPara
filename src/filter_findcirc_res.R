#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                help = "Findcirc result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "integer",
                default = 40,
                help = "Minimum quality of backsplice end mapping"),
    make_option(c("-o", "--output"), action = "store", type = "character",
                default = "filteredFindcirc.bed", help = "The output file")
)

parser <- OptionParser(usage = "%prog -i circ_candidates.bed -q 40 -o filteredFindcirc.bed",
                       option_list = option_list,
                       description = "Filter Findcirc results by quality")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
output <- arguments$output
minqual <- arguments$minqual

circ_candidates.bed <-
    fread(input, showProgress = F, skip = 1)[
        grepl('CIRCULAR', V18)][
            grepl('UNAMBIGUOUS_BP', V18)][
                grepl('ANCHOR_UNIQUE', V18)][
                    V9 >= minqual & V10 >= minqual]

write.table(x = circ_candidates.bed,
            file = output, row.names = F, quote = F, sep = "\t", col.names = F)
