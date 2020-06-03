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
    make_option(c("-t", "--tolerance"), action="store", type="integer", default = 5,
                help="Number of basepairs tolerated in realigning circRNAs from CIRCexplorer2 annotate"),
    make_option(c("-s", "--stranded"), action="store_true", default = FALSE,
                help="Set if stranded library")
)

parser <-
    OptionParser(usage="%prog -r Chimeric.out.junction -c strandedCircRNACount -g 10 -o circular.reads.bed.gz",
                 option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

orig.file <- arguments$circrnas
orig.est <- fread(file = orig.file)[, .(V1 = as.character(Chr), V2 = Start, V3 = End,
                                        orig = Chimeric.out.junction, V6 = Strand), ]

## process the Chimeric.out.junction entries as DCC does:
## see https://github.com/dieterich-lab/DCC/blob/master/DCC/findcircRNA.py
## findcirc and cigarGenomicDist functions

## read the file and keep only rows in which
## 1. chromosome of the segments is the same
## 2. strand of the segments is the same
chimout.file <- arguments$chimreads
chimout <- fread(chimout.file, verbose = F,
                 showProgress = F, header = F)[V1 == V4 & V3 == V6 & V7 >= 0]
chimout$V1 <- as.character(chimout$V1)

cigar_dist <- function(cigar.str){
    cigar.split <- strsplit(cigar.str, "")[[1]]
    C <- grep('[a-zA-Z]', cigar.split, value = T)
    L <- gregexpr('-?[0-9]+', cigar.str)
    g <- 0
    for(i in 1:length(L[[1]])){
        if(C[i] != 'S' & C[i] != "I"){
            g <- g + as.integer(paste0(cigar.split[L[[1]][i]:(L[[1]][i]+attr(L[[1]],
                                                                             "match.length")[i]-1)],
                                       collapse = ""))
        }
    }
    g
}

chimout[, `:=`(V14.cigar.dist = sapply(V14, cigar_dist),
               V12.cigar.dist = sapply(V12, cigar_dist))]
## check orientation of segment alignement:
## if + strand then first segment comes after second segment, first->second if - strand
## Mind the length of aligned segment, and tolerance (5 bases default)

tolerance <- arguments$tolerance

chimout.junc.bed <-
    rbindlist(list(chimout[V3 == "+" & V11 + tolerance > V5 & (V13 - tolerance) + V14.cigar.dist <= V2],
                   chimout[V3 == "-" & V13 + tolerance > V2 & (V11 - tolerance) + V12.cigar.dist <= V5]),
              use.names = T)

## swap start/end if + strand
chimout.junc.bed[V3 == "+", `:=`(tmp = V2, V2 = V5)][V3 == "+", V5 := tmp]

## prepare BED coordinates:
## Chimeric.out.junction coordinates are 1-based, but refer to the intron first position
## add 1 to start (acceptor intron base) and remove 1 to end (donor intron base)
chimout.junc.bed <-
    chimout.junc.bed[, .(V1, V2 = V2 + 1, V3 = V5 - 1, V4 = V10, V5 = V7, V6 = V3)]

strand <- arguments$stranded
invert.strand <- function(s){
    if(s == "+"){
        s <- "-"
    }else{
        s <- "+"
    }
    s
}

if(strand){
    chimout.junc.bed[, V6 := sapply(V6, invert.strand)]

    filterd.chimout.junc <-
        merge(chimout.junc.bed,
              orig.est[, .(V1, V2, V3, V6)],
              by = c("V1", "V2", "V3", "V6"),
              all.x = F,
              all.y = T)[, .(V1, V2, V3, V4, V5, V6)]
}else{
    filterd.chimout.junc <-
        merge(chimout.junc.bed,
              orig.est[, .(V1, V2, V3)],
              by = c("V1", "V2", "V3"),
              all.x = F,
              all.y = T)[, .(V1, V2, V3, V4, V5, V6)]
}

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
