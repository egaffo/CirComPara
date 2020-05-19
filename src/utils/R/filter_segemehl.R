#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "sample.sngl.bed",
                help = "Segemehl result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "character",
                default = "median_10",
                help = paste("Filter function and minimum mapping quality of backsplice-reads. ",
                       "Function supported: ",
                       "- median: for each backsplice compute the read mapping median quality and discard backsplices with median quality <= N",
                       "- mean: as for median, but compute the mean instead",
                       "- any: first, discard any alignment with mapping quality <= N. Then, report the median as backsplice quality",
                       sep = "\n")),
    make_option(c("-o", "--count_output"), action = "store", type = "character",
                default = "splicesites.bed",
                help = "The filtered circrnas in BED format. The 7th column is just a copy of the score field and reports the read count."),
    make_option(c("-r", "--reads_output"), action = "store", type = "character",
                default = "sample.circular.reads.bed.gz",
                help = "The output file"),
    make_option(c("-m", "--keep_mates"), action = "store_true", #type = "logical",
                default = F,
                help = "Whether to include /1 or /2 in read names to distinguish readmates")
)

parser <- OptionParser(usage = "%prog -i sample.sngl.bed -q median_10 -o splicesites.bed -r sample.circular.reads.bed.gz",
                       option_list = option_list,
                       description = "Filter segemehl results by alignments' quality and compute backsplices' read count")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
count_output <- arguments$count_output
reads_output <- arguments$reads_output
min.qual <- arguments$minqual

split.qual.par <- strsplit(min.qual, "_", fixed = T)[[1]]
if(length(split.qual.par) > 1){
    qual.func <- split.qual.par[1]
    if(!qual.func %in% c("median", "mean", "any")){
        message("Error in evaluating quality function in '", min.qual, "': setting default 'median'")
        qual.func <- "median"
    }
    minqual <- suppressWarnings(as.numeric(split.qual.par[2]))
    if(is.na(minqual)){
        message("Error in evaluating quality value in '", min.qual, "': setting default '10'")
        minqual <- 10
    }
}else{
    message("Error in evaluating quality function '", min.qual, "': setting default 'median_10'")
    qual.func <- "median"
    minqual <- 10
}

# Start and end position indicate the genomic range of the predicted intron.
# The name has the format (read-group;type;read-name;mate-status), the bed
# score is the alignment score of the respective alignment. The type is either
# 'R' (in case of a regular, collinear split), 'C' (circular split) or 'B' (backsplice)

sege_circ <- fread(cmd = paste0('grep ";B\\|C;" ', input),
                   header = F, skip = 1)

splicesites.bed <- data.table()
reads_output.dt <- data.table()

if(qual.func == "any"){
    sege_circ <- sege_circ[V5 >= minqual]
}

if(nrow(sege_circ) > 0){

    if(arguments$keep_mates){
        sege_circ[, c("read.group", "type", "read.name",
                      "mate.status"):=(tstrsplit(V4, ";"))][, read.name :=
                                                                paste0(read.name, "/",
                                                                       mate.status)][, `:=`(V4 = NULL,
                                                                                            read.group = NULL,
                                                                                            type = NULL,
                                                                                            mate.status = NULL)]
    }else{
        sege_circ[, c("read.group", "type", "read.name",
                      "mate.status"):=(tstrsplit(V4, ";"))][, `:=`(V4 = NULL,
                                                                   read.group = NULL,
                                                                   type = NULL,
                                                                   mate.status = NULL)]
    }


    ## remove duplicated lines/alignments
    sege_circ <- sege_circ[, .(multi.mapping = .N,
                               map.qual = max(V5)),
                           by = .(chr = V1, left = V2,
                                  right = V3, read.name,
                                  strand = V6)]

    if(qual.func == "any" || qual.func == "median"){
        splicesites.bed <-
            sege_circ[, .(n = .N,
                          map.qual = median(map.qual)),
                      by = .(chr, left, right,
                             strand)][, .(chr, left, right, n, map.qual, strand,
                                          score = n)][order(chr, left,
                                                            right)]
    }

    if(qual.func == "mean"){
        splicesites.bed <-
            sege_circ[, .(n = .N,
                          map.qual = mean(map.qual)),
                      by = .(chr, left, right,
                             strand)][, .(chr, left, right, n, map.qual, strand,
                                          score = n)][order(chr, left,
                                                            right)]
    }

    if(qual.func == "median" || qual.func == "mean"){
        ## filter by quality
        splicesites.bed <- splicesites.bed[map.qual >= minqual]
    }

    ## report selected backsplices' read names
    # reads_output.dt <- sege_circ[, .(chr, left, right, read.name, map.qual, strand)]
    reads_output.dt <- sege_circ[splicesites.bed[, .(chr, left, right, strand)],
                                 on = c("chr", "left", "right",
                                        "strand")][, .(chr, left, right, read.name,
                                                       map.qual, strand)]
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
