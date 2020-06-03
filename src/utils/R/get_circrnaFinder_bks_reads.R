#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

option_list <- list(
    make_option(c("-r", "--chimreads"), action="store", type="character",
                help="Chimeric.out.junction BED filtered and formatted by cf_filterChimout.awk"),
    make_option(c("-c", "--circrnas"), action="store", type="character",
                help="circularRNA_known.txt as output by CIRCexplorer2 annotate, or back_spliced_junction.bed as output by CIRCexplorer2 parse"),
    make_option(c("-o", "--output"), action="store", type="character",
                help="The circRNA read IDs for each circRNA in compressed BED (circular.reads.bed.gz)")
)

parser <-
    OptionParser(usage="%prog -r cf.filtered.Chimeric.out.junction -c filteredJunctions.bed -o circular.reads.bed.gz",
                 option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

orig.file <- arguments$circrnas
orig.est <- fread(orig.file)[, .(orig = sum(V5)), by = .(V1, V2, V3, V6), ]
orig.est$V1 <- as.character(orig.est$V1)

## process the Chimeric.out.junction entries filtered as circrna_finder does:
## see https://github.com/orzechoj/circRNA_finder/blob/master/filterCirc.awk
## Assume chimeric reads pre-processed by the cf_filterChimout.awk form CirComPara

chimout.file <- arguments$chimreads
chimout.junc.bed <- fread(chimout.file, verbose = F,
                          showProgress = F, header = F)
chimout.junc.bed$V1 <- as.character(chimout.junc.bed$V1)

filterd.chimout.junc <-
    merge(chimout.junc.bed,
          orig.est[, .(V1, V2, V3, V6)],
          by = c("V1", "V2", "V3", "V6"),
          all.x = F,
          all.y = T)[, .(V1, V2, V3, V4, V5, V6)]

splitted.filename <- strsplit(arguments$output, ".", fixed = T)[[1]]
if(tail(splitted.filename, 1) == "gz"){
    tmp.outfile <- sub(".gz$", "", arguments$output)
}

fwrite(x = filterd.chimout.junc,
       file = tmp.outfile,
       sep = "\t",
       col.names = F,
       row.names = F)

if(tail(splitted.filename, 1) == "gz"){
    gzip(tmp.outfile, destname = arguments$output)
}
