#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

option_list <- list(
    make_option(c("-r", "--chimreads"), action="store", type="character",
                help="Chimeric.out.junction BED transformed by chimout_junc_to_bed.py"),
    make_option(c("-c", "--circrnas"), action="store", type="character",
                help="circularRNA_known.txt as output by CIRCexplorer2 annotate, or back_spliced_junction.bed as output by CIRCexplorer2 parse"),
    make_option(c("-o", "--output"), action="store", type="character",
                help="The circRNA read IDs for each circRNA in compressed BED (circular.reads.bed.gz)"),
    make_option(c("-g", "--range"), action="store", type="integer", default = 10,
                help="Number of basepairs tolerated in realigning circRNAs from CIRCexplorer2 annotate")
)

parser <-
    OptionParser(usage="%prog -r Chimeric.out.junction.bed -c circularRNA_known.txt -g 10 -o circular.reads.bed.gz",
                 option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

orig.file <- arguments$circrnas
orig.est <- fread(file = orig.file)
annotation <- "bed"
if(ncol(orig.est) > 6){
    annotation <- "annotated"
}

orig.est <- orig.est[, .(V1, V2, V3, V5, V6)]

chimout.junc.file <- arguments$chimreads
chimout.junc.bed <- fread(file = chimout.junc.file, header = F)

if(annotation == "annotated"){
    ext.range <- c(-arguments$range:arguments$range)
    ext.range <- ext.range[ext.range != 0]
    unfix.orig.est <-
        orig.est[V5 != 0, .(start = V2 + ext.range,
                            end = V3 + ext.range),
                 by = .(V1, fixed.V2 = V2, fixed.V3= V3, V6)]

    fixed.chimout.junc <-
        merge(chimout.junc.bed,
              unfix.orig.est,
              by.x = c("V1", "V2", "V3"),
              by.y = c("V1", "start", "end"))[, .(V1, V2 = fixed.V2,
                                                  V3 = fixed.V3, V4, V5, V6 = V6.y)]


    unfixed.chimout.junc <-
        merge(chimout.junc.bed,
              orig.est[V5 == 0],
              by = c("V1", "V2", "V3"))[, .(V1, V2, V3, V4, V5 = V5.y, V6 = V6.y)]

    annotated.chimout.junc <-
        rbindlist(list(fixed.chimout.junc,
                       unfixed.chimout.junc),
                  use.names = T)
}else{
    annotated.chimout.junc <-
        merge(chimout.junc.bed[, .(V1, V2, V3, V4)],
              orig.est[, .(V1, V2, V3, V5, V6)],
              by = c("V1", "V2", "V3"), all.x = F, all.y = T)
}

splitted.filename <- strsplit(arguments$output, ".", fixed = T)[[1]]
if(tail(splitted.filename, 1) == "gz"){
    tmp.outfile <- sub(".gz$", "", arguments$output)
}

fwrite(x = annotated.chimout.junc,
       file = tmp.outfile,
       sep = "\t",
       col.names = F,
       row.names = F)

if(tail(splitted.filename, 1) == "gz"){
    gzip(tmp.outfile, destname = arguments$output)
}
