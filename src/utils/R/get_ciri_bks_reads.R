#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "ciri.out",
                help = "CIRI2 output file (default: ciri.out)"),
    make_option(c("-b", "--reads_bed"), action = "store", type = "character",
                default = "circular.reads.bed.gz",
                help = "Output filename of BED file with read IDs for each backsplice (default: circular.reads.bed.gz)"),
    make_option(c("-l", "--read_list"), action = "store", type = "character",
                default = "bks.reads",
                help = "Output filename of all read ID list (default: bks.reads)")
)

parser <- OptionParser(usage = "%prog -i ciri.out -b circular.reads.bed.gz -l bks.reads",
                       option_list = option_list,
                       description = "")

arguments <- parse_args(parser, positional_arguments = F)

input <- arguments$input
reads.bed <- arguments$reads_bed
read.list <- arguments$read_list

# input <- "/blackhole/enrico/circular/circompara_testing/circompara/test_circompara/analysis/samples/sample_A/processings/circRNAs/ciri_out/sample_A_ciri.out"

## N.B: CIRI output is 1-based as GTFs.
## We need to decrease start position to comply with BED format
ciri.out <- fread(input = input, select = c(2,3,4,11,12), showProgress = F)[, circRNA_start := circRNA_start - 1]

## method jaap_DT2 from
## https://stackoverflow.com/questions/13773770/split-comma-separated-strings-in-a-column-into-separate-rows
bks.reads <-
    ciri.out[, strsplit(as.character(junction_reads_ID), ",", fixed = T),
         by = .(chr, circRNA_start, circRNA_end, strand,
                junction_reads_ID)][, .(chr, circRNA_start, circRNA_end,
                                        read_id = V1, score = 0, strand)]

## write gzipped file for circular reads
reads_output.gz <- gzfile(reads.bed, "w")
write.table(bks.reads, file = reads_output.gz,
            sep = "\t", col.names = F, row.names = F, quote = F)
close(reads_output.gz)

all.reads <-
    bks.reads[, .N, by = read_id][order(-N), .(N, read_id)]

write.table(all.reads, file = read.list,
            sep = "\t", col.names = F, row.names = F, quote = F)
