#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "circular.reads.bed.gz.txt",
                help = "File with the list of circular.reads.bed.gz files to merge, each sample and method (default: circular.reads.bed.gz.txt)"),
    make_option(c("-q", "--min_methods"), action = "store", type = "integer",
                default = 2,
                help = "Minimum number of methods (default: 2)"),
    make_option(c("-o", "--output_prefix"), action = "store", type = "character",
                default = "bks.counts.",
                help = "A prefix for output file names (default: bks.counts.intersect.csv, bks.counts.union.csv, bks.counts.union.intersected.csv)"),
    # make_option(c("-f", "--format"), action = "store", type = "character",
    #             default = "BED",
    #             help = "The file format of output files {GTF,BED}"),
    make_option(c("-s", "--stranded"), action = "store_true", default = F,
                help = "Set this option if circRNA strand has to be considered")
)

parser <- OptionParser(usage = "%prog -i circular.reads.bed.gz.txt -q 2 -o bks.counts.",
                       option_list = option_list,
                       description = "Compute backsplices' read counts")

arguments <- parse_args(parser, positional_arguments = F)

input <- arguments$input
output_prefix <- arguments$output_prefix
min_methods <- arguments$min_methods
# file.ext <- tolower(arguments$format)
file.ext <- "csv"
stranded <- arguments$stranded

circular.reads.bed.gz.txt <- readLines(input)

bks.read.method <-
    unique(rbindlist(sapply(circular.reads.bed.gz.txt, fread, header = F,
                            col.names = c("chr", "start", "end", "read_id", "score", "strand"),
                            simplify = F),
                     idcol = "sample_id")[, `:=`(sample_id = sub(".circular.reads.bed.gz", "",
                                                                 basename(sample_id)),
                                                 circ_method = sub(".+circRNAs/([^/]+)/.*", "\\1",
                                                                   sample_id))][])
## fix DCC coordinates
bks.read.method[circ_method == "dcc", start := start -1L]

## save collected reads
filename <- paste0(output_prefix, "collected_reads.", file.ext)
write.csv(x = bks.read.method[, .(sample_id, chr, start, end, read_id, strand, circ_method)],
          file = filename,
          row.names = F)

##
if(stranded){
    bks.reads <-
        bks.read.method[, .(n_methods = length(unique(circ_method)),
                            methods = paste0(sort(unique(circ_method)), collapse = "|")),
                        by = .(sample_id, chr, start, end, strand, read_id)]
}else{
    bks.reads <-
        bks.read.method[, .(n_methods = length(unique(circ_method)),
                            methods = paste0(sort(unique(circ_method)), collapse = "|")),
                        by = .(sample_id, chr, start, end, read_id)]
}

## STRATEGY A)
## for each backsplice, count only reads detected by >= min_methods
## This should improve reliability of read counts
## n_methods_partials: comma-separated list of read.count@n_methods.
##                     Example
##                     33@5,56@2,8@3,1@2 means that for 98 reads in total:
##                     - 33 reads were commonly detected by 5 methods,
##                     - 56 reads were commonly detected by 2 methods,
##                     - 8 reads were commonly detected by 3 methods,
##                     - 1 reads were commonly detected by 2 methods (in a different combination with respect to the other 56 reads)
## methods_partials: comma-separated list of read.count@method_names.
##                   As for n_methods_partials, giving the method names combinations
if(stranded){
    bks.read.counts.intersect <-
        bks.reads[n_methods >= min_methods, .N,
                  by = .(sample_id, chr, start, end, strand, methods,
                         n_methods)][, .(read.count = sum(N),
                                         n_methods_partials = paste0(paste0(N, "@", n_methods),
                                                                     collapse = ","),
                                         methods_partials = paste0(paste0(N, "@", methods),
                                                                   collapse = ",")),
                                     by = .(sample_id, chr, start, end, strand)]
}else{
    bks.read.counts.intersect <-
        bks.reads[n_methods >= min_methods, .N,
                  by = .(sample_id, chr, start, end, methods,
                         n_methods)][, .(read.count = sum(N),
                                         n_methods_partials = paste0(paste0(N, "@", n_methods),
                                                                     collapse = ","),
                                         methods_partials = paste0(paste0(N, "@", methods),
                                                                   collapse = ",")),
                                     by = .(sample_id, chr, start, end)]
}

filename <- paste0(output_prefix, "intersect.", file.ext)
write.csv(x = bks.read.counts.intersect,
          file = filename,
          row.names = F)

## STRATEGY B)
## count all reads for each backsplice and keep it if >= min_methods detected
## any read for it, i.e: keep the backsplice and its read count if different
## methods detected it also if by different read ids
## F.i, if min_methods = 2:
## - keep bks Y with count "100 reads" in total: 50 reads are detected only
## by DCC and 50 only by CIRI, no reads detected by both methods
## - discard bks Z with 100 reads detected only by DCC
## This should improve sensibility of detection
if(stranded){
    bks.read.counts.union <-
        bks.read.method[, .(circ_methods = list(unique(circ_method))),
                        by = .(sample_id, chr, start, end, strand,
                               read_id)][, .(read.count = .N,
                                             circ_methods = list(unique(unlist(circ_methods)))),
                                         by = .(sample_id, chr, start, end,
                                                strand)][, .(n_methods = length(unlist(circ_methods)),
                                                             circ_methods = paste0(sort(unlist(circ_methods)), collapse = "|")),
                                                         by = .(sample_id, chr, start, end,
                                                                strand, read.count)] #[n_methods >= min_methods] ## do not filter already
}else{
    bks.read.counts.union <-
        bks.read.method[, .(circ_methods = list(unique(circ_method))),
                        by = .(sample_id, chr, start, end,
                               read_id)][, .(read.count = .N,
                                             circ_methods = list(unique(unlist(circ_methods)))),
                                         by = .(sample_id, chr, start, end)][, .(n_methods = length(unlist(circ_methods)),
                                                                                 circ_methods = paste0(sort(unlist(circ_methods)), collapse = "|")),
                                                                             by = .(sample_id, chr, start, end,
                                                                                    read.count)] #[n_methods >= min_methods] ## do not filter already
}

filename <- paste0(output_prefix, "union.", file.ext)
write.csv(x = bks.read.counts.union,
          file = filename,
          row.names = F)

## STRATEGY C)
## for each backsplice, if one read was detected by >= min_methods then count
## also reads detected by < min_methods
## F.i, if min_methods = 2:
## - keep bks Y with count "100 reads" in total: 50 reads are detected only
## by DCC and 49 only by CIRI, plus 1 read which is detected by both methods;
## - discard bks Z with count "100 reads" in total: 50 reads are detected only
## by DCC and 50 only by CIRI, no reads detected by both methods
if(stranded){
    bks <- unique(bks.reads[n_methods >= min_methods, .(chr, start, end, strand)])
    bks.read.counts.union.intersected <-
        bks.reads[bks, on = c("chr", "start", "end", "strand")][, .(read.count = .N),
                                                                by = .(sample_id, chr, start, end, strand)]
}else{
    bks <- unique(bks.reads[n_methods >= min_methods, .(chr, start, end)])
    bks.read.counts.union.intersected <-
        bks.reads[bks, on = c("chr", "start", "end")][, .(read.count = .N),
                                                      by = .(sample_id, chr, start, end)]

}

filename <- paste0(output_prefix, "union.intersected.", file.ext)
write.csv(x = bks.read.counts.union.intersected,
          file = filename,
          row.names = F)

# ## save in a matrix-like format
# bks.read.counts <- bks.read.counts.union
# bks.read.counts.wide <-
#     dcast(data = bks.read.counts,
#           formula = chr + start + end + strand ~ sample_id,
#           value.var = "read.count",
#           fill = 0)
#
# filename <- paste0(output_prefix, "")
# write.csv(x = ,
#           file = filename,
#           row.names = F)
