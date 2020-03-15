#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

option_list <- list(
    make_option(c("-r", "--chimreads"), action="store", type="character",
                help="sample.sngl.bed "),
    make_option(c("-c", "--circrnas"), action="store", type="character",
                help="circularRNA_known.txt as output by CIRCexplorer2 annotate, or back_spliced_junction.bed as output by CIRCexplorer2 parse"),
    make_option(c("-o", "--output"), action="store", type="character",
                help="The circRNA read IDs for each circRNA in compressed BED (circular.reads.bed.gz)")
)

parser <-
    OptionParser(usage="%prog -r sample.sngl.bed -c circularRNA_known.txt -o circular.reads.bed.gz",
                 option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

switch.strand <- function(s){
    if(s == "+"){
        s <- "-"
    }else{
        s <- "+"
    }
    s
}

orig.file <- arguments$circrnas
orig.est <- fread(file = orig.file)
annotation <- "bed"
if(ncol(orig.est) > 6){
    annotation <- "annotated"
}

orig.est <- orig.est[, .(V1, V2, V3, V5, V6)]

seg.bks.reads.file <- arguments$chimreads
seg.bks.reads <-
    fread(seg.bks.reads.file, showProgress = F,
          skip = 1)[grepl(";C|B;", V4)]
seg.bks.reads <-
    seg.bks.reads[, c("read.group", "type", "read.name",
                      "mate.status") :=
                      tstrsplit(V4, ";",
                                type.convert = T)][, read.name :=
                                                       paste0(read.name, "/",
                                                              mate.status)][, `:=`(V4 = NULL,
                                                                                   read.group = NULL,
                                                                                   mate.status = NULL)]

unfixed.seg.bks.reads <-
    merge(seg.bks.reads,
          orig.est[V5 == 0],
          by = c("V1", "V2", "V3",
                 "V6"))[, .(V1, V2, V3,
                            V4 = read.name, V5 = V5.y,
                            V6)]

## circrnas missed
# nrow(orig.est[V5 == 0]) - nrow(unfixed.seg.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)])

## check if the missed circrnas have just a diferent strand
unmatched <-
    merge(orig.est[V5 == 0],
          unfixed.seg.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)],
          by = c("V1", "V2", "V3", "V6"),
          all = T)[is.na(ccp)]

unmatched$V6 <- sapply(unmatched$V6, switch.strand)
strand.switched <-
    merge(unmatched,
          seg.bks.reads,
          by = c("V1", "V2", "V3", "V6"),
          all.y = F)[, .(V1, V2, V3,
                         V4 = read.name, V5 = V5.y,
                         V6)]

## missed circrnas fixed by strand switch
# nrow(strand.switched[, .(ccp = .N), by = .(V1, V2, V3, V6)])

## check if now the matched circrnas correspond to the
## circexplorer-unmodified circrnas
# nrow(orig.est[V5 == 0]) == nrow(rbindlist(list(unfixed.seg.bks.reads,
#                                                strand.switched),
#                                           use.names = T)[, .(ccp = .N), by = .(V1, V2, V3, V6)])

## now, check the circexplorer-MODIFIED circrnas

## number of ce2-modified circrnas
# nrow(orig.est[V5 != 0])

## prepare tolerance intervals for coordinate fix
ext.range <- c(-10:10)
ext.range <- ext.range[ext.range != 0]
unfix.orig.est <-
    orig.est[V5 != 0, .(start = V2 + ext.range,
                        end = V3 + ext.range),
             by = .(V1, fixed.V2 = V2, fixed.V3= V3, V6)]

## search the reads with the extended coordinates
## fixed.seg.bks.reads will be the extended coordinates-found reads
fixed.seg.bks.reads <-
    merge(seg.bks.reads,
          unfix.orig.est,
          by.x = c("V1", "V2", "V3", "V6"),
          by.y = c("V1", "start", "end",
                   "V6"))[, .(V1, V2 = fixed.V2,
                              V3 = fixed.V3, V4 = read.name,
                              V5, V6)]

## check if some circrna is still missed also after coordinate extension search
fixed.match <-
    merge(orig.est[V5 != 0],
          fixed.seg.bks.reads[, .N, by = .(V1, V2, V3, V6)],
          by = c("V1", "V2", "V3", "V6"),
          all = T)#[orig - N != 0]

found <- fixed.match[!is.na(N)]

still.missed <- fixed.match[is.na(N)]

## check if the still missed have a different strand in reads
still.missed$V6 <- sapply(still.missed$V6, switch.strand)

## expand interval of still missed
extended.still.missed <-
    still.missed[, .(start = V2 + ext.range,
                     end = V3 + ext.range),
                 by = .(V1, fixed.V2 = V2, fixed.V3= V3, V6)]

still.missed.strand.switched <-
    merge(seg.bks.reads,
          extended.still.missed,
          by.x = c("V1", "V2", "V3", "V6"),
          by.y = c("V1", "start", "end",
                   "V6"))[, .(V1, V2 = fixed.V2,
                              V3 = fixed.V3, V4 = read.name,
                              V5, V6)]

## missed circrnas fixed by strand switch
# nrow(still.missed.strand.switched[, .(ccp = .N), by = .(V1, V2, V3, V6)])

## bind all reads (all.seg.bks.reads)
# annotated.seg.bks.reads <-
#     rbindlist(list(unfixed.seg.bks.reads,
#                    strand.switched,
#                    fixed.seg.bks.reads,
#                    still.missed.strand.switched),
#               use.names = T)

## check the number of identified circrnas
# nrow(annotated.seg.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)]) == nrow(orig.est)

# ccp.est <- annotated.seg.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)]
ccp.est.reads <-
    rbindlist(list(unf = unfixed.seg.bks.reads,
                   unf.str.swi = strand.switched[, .(V1, V2, V3, V4, V5,
                                                     V6 = sapply(V6, switch.strand))],
                   fix.coo = merge(seg.bks.reads,
                                   unfix.orig.est,
                                   by.x = c("V1", "V2", "V3", "V6"),
                                   by.y = c("V1", "start", "end",
                                            "V6"),
                                   all = F)[, .(V1, V2 = fixed.V2,
                                                V3 = fixed.V3, V4 = read.name,
                                                V5, V6)],
                   fix.coo.str = merge(seg.bks.reads,
                                       extended.still.missed,
                                       by.x = c("V1", "V2", "V3", "V6"),
                                       by.y = c("V1", "start", "end",
                                                "V6"))[, .(V1, V2 = fixed.V2,
                                                           V3 = fixed.V3, V4 = read.name,
                                                           V5, V6 = sapply(V6,
                                                                           switch.strand))]),
              use.names = T,
              idcol = "Set")

splitted.filename <- strsplit(arguments$output, ".", fixed = T)[[1]]
if(tail(splitted.filename, 1) == "gz"){
    tmp.outfile <- sub(".gz$", "", arguments$output)
}

fwrite(x = ccp.est.reads[, .(V1, V2, V3, V4, V5, V6)],
       file = tmp.outfile,
       sep = "\t",
       col.names = F,
       row.names = F)

if(tail(splitted.filename, 1) == "gz"){
    gzip(tmp.outfile, destname = arguments$output)
}
