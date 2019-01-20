#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "sample.sngl.bed",
                help = "Segemehl result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "integer",
                default = 25,
                help = "Minimum mapping quality of backsplice-reads"),
    make_option(c("-o", "--count_output"), action = "store", type = "character",
                default = "splicesites.bed",
                help = "The filtered circrnas in BED format. The 7th column is just a copy of the score field and reports the read count."),
    make_option(c("-r", "--reads_output"), action = "store", type = "character",
                default = "sample.circular.reads.bed.gz",
                help = "The output file")
)

parser <- OptionParser(usage = "%prog -i sample.sngl.bed -q 25 -o splicesites.bed",
                       option_list = option_list,
                       description = "Filter segemehl results by quality and compute backsplices' counts")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
count_output <- arguments$count_output
reads_output <- arguments$reads_output
minqual <- arguments$minqual

# Start and end position indicate the genomic range of the predicted intron.
# The name has the format (read-group;type;read-name;mate-status), the bed
# score is the alignment score of the respective alignment. The type is either
# 'R' (in case of a regular, collinear split), 'C' (circular split) or 'B' (backsplice)

sege_circ <- fread(cmd = paste0('grep ";B\\|C;" ', input),
                   header = F, skip = 1)[V5 >= minqual]

splicesites.bed <- data.table()
reads_output.dt <- data.table()

if(nrow(sege_circ) > 0){

    sege_circ[, c("read.group", "type", "read.name", "mate.status"):=(tstrsplit(V4, ";"))]
    sege_circ <- sege_circ[, .(multi.mapping = .N, map.qual = median(V5)),
                           by = .(chr = V1, left = V2,
                                  right = V3, read.name,
                                  strand = V6)]

    splicesites.bed <-
        sege_circ[, .(n = .N,
                      median_qual = median(map.qual)),
                  by = .(chr, left, right,
                         strand)][, .(chr, left, right, n, median_qual, strand,
                                      score = n)][order(chr, left,
                                                        right)]

    reads_output.dt <- sege_circ[, .(chr, left, right, read.name, map.qual, strand)]
}

## write backsplice counts
write.table(x = splicesites.bed,
            file = count_output,
            row.names = F, quote = F, sep = "\t", col.names = T)

## write gzipped file for circular reads
reads_output.gz <- gzfile(reads_output, "w")
write.table(x = reads_output.dt,
            file = reads_output.gz,
            row.names = F, quote = F, sep = "\t", col.names = F)
close(reads_output.gz)
