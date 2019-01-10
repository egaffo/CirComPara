#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "sample.sngl.bed",
                help = "Segemehl result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "integer",
                default = 25,
                help = "Minimum quality of backsplice end mapping"),
    make_option(c("-o", "--output"), action = "store", type = "character",
                default = "splicesites.bed",
                help = "The output file")
)

parser <- OptionParser(usage = "%prog -i sample.sngl.bed -q 25 -o splicesites.bed",
                       option_list = option_list,
                       description = "Filter segemehl results by quality and compute backsplices' counts")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
# input <- "/blackhole/enrico/circular/circompara_testing/circompara/test_circompara/analysis/samples/sample_A/processings/circRNAs/segemehl/sample_A.unmappedSE.fq.sngl.bed"
output <- arguments$output
minqual <- arguments$minqual

sege_circ <- fread(input, header = F, skip = 1)

if(nrow(sege_circ) > 0){
    sege_circ <-
        sege_circ[grepl(";B|C;",
                        V4)][, .(n = .N,
                                 median_qual = median(V5)),
                             by = .(chr = V1, left = V2, right = V3,
                                    strand = V6)][, .(chr, left, right, n, median_qual,
                                                      strand)][order(left, right, n, median_qual,
                                                                     strand)][median_qual >= minqual]
}

write.table(x = sege_circ,
            file = output,
            row.names = F, quote = F, sep = "\t", col.names = T)
