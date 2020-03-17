#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

option_list <- list(
    make_option(c("-r", "--chimreads"), action="store", type="character",
                help="back_spliced_junction.bed BED file output by CIRCexplorer parse TopHat >= v2.3.5"),
    make_option(c("-c", "--circrnas"), action="store", type="character",
                help="circularRNA_known.txt as output by CIRCexplorer2 annotate, or back_spliced_junction.bed as output by CIRCexplorer2 parse"),
    make_option(c("-o", "--output"), action="store", type="character",
                help="The circRNA read IDs for each circRNA in compressed BED (circular.reads.bed.gz)"),
    make_option(c("-g", "--range"), action="store", type="integer", default = 10,
                help="Number of basepairs tolerated in realigning circRNAs from CIRCexplorer2 annotate")
)

parser <-
    OptionParser(usage="%prog -r back_spliced_junction.bed -c circularRNA_known.txt -o circular.reads.bed.gz",
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

bks.reads.file <- arguments$chimreads
bks.reads <- fread(bks.reads.file, header = F)[, V4 := V7][, V7 := NULL][]
bks.reads <-
    bks.reads[, strsplit(V4, ",", fixed = TRUE),
              by = .(chr = V1, V2, V3, V5,
                     V6)][, c("V4",
                              "th.cig") := tstrsplit(V1, ";",
                                                     fixed = T)][, `:=`(V1 = chr,
                                                                        th.cig = NULL)][, chr := NULL][]

if(annotation == "annotated"){

    unfixed.bks.reads <-
        merge(bks.reads,
              orig.est[V5 == 0],
              by = c("V1", "V2", "V3",
                     "V6"))[, .(V1, V2, V3, V4, V5 = V5.y, V6)]

    chimout.slices <- list(unf = unfixed.bks.reads)

    ## circrnas missed
    # nrow(orig.est[V5 == 0]) - nrow(unfixed.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)])
    ## check if the missed circrnas have just a diferent strand
    unmatched <-
        merge(orig.est[V5 == 0],
              unfixed.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)],
              by = c("V1", "V2", "V3", "V6"),
              all = T)[is.na(ccp)]

    if(nrow(unmatched) > 0){
        unmatched$V6 <- sapply(unmatched$V6, switch.strand)
        strand.switched <-
            merge(unmatched,
                  bks.reads,
                  by = c("V1", "V2", "V3", "V6"),
                  all.y = F)[, .(V1, V2, V3,
                                 V4, V5 = V5.y,
                                 V6)]

        ## missed circrnas fixed by strand switch
        # nrow(strand.switched[, .(ccp = .N), by = .(V1, V2, V3, V6)])
        if(nrow(strand.switched[, .(ccp = .N), by = .(V1, V2, V3, V6)]) > 0){
            chimout.slices$unf.str.swi <-
                strand.switched[, .(V1, V2, V3, V4, V5,
                                    V6 = sapply(V6, switch.strand))]
        }

        ## check if now the matched circrnas correspond to the
        ## circexplorer-unmodified circrnas
        # nrow(orig.est[V5 == 0]) == nrow(rbindlist(list(unfixed.bks.reads,
        #                                                strand.switched),
        #                                           use.names = T)[, .(ccp = .N), by = .(V1, V2, V3, V6)])
    }
    ## now, check the circexplorer-MODIFIED circrnas

    ## number of ce2-modified circrnas
    # nrow(orig.est[V5 != 0])

    if(nrow(orig.est[V5 != 0]) > 0){
        ## prepare tolerance intervals for coordinate fix
        ext.range <- c(-arguments$range:arguments$range)
        ext.range <- ext.range[ext.range != 0]
        unfix.orig.est <-
            orig.est[V5 != 0, .(start = V2 + ext.range,
                                end = V3 + ext.range),
                     by = .(V1, fixed.V2 = V2, fixed.V3= V3, V6)]

        ## search the reads with the extended coordinates
        ## fixed.bks.reads will be the extended coordinates-found reads
        fixed.bks.reads <-
            merge(bks.reads,
                  unfix.orig.est,
                  by.x = c("V1", "V2", "V3", "V6"),
                  by.y = c("V1", "start", "end",
                           "V6"))[, .(V1, V2 = fixed.V2,
                                      V3 = fixed.V3, V4,
                                      V5, V6)]

        fixed.match <-
            merge(orig.est[V5 != 0],
                  fixed.bks.reads[, .N, by = .(V1, V2, V3, V6)],
                  by = c("V1", "V2", "V3", "V6"),
                  all = T)

        if(nrow(fixed.match[!is.na(N)]) > 0){
            ## store fixed coordinates for reads
            chimout.slices$fix.coo <-
                merge(bks.reads,
                      unfix.orig.est,
                      by.x = c("V1", "V2", "V3", "V6"),
                      by.y = c("V1", "start", "end",
                               "V6"),
                      all = F)[, .(V1, V2 = fixed.V2,
                                   V3 = fixed.V3, V4,
                                   V5, V6)]
        }

        ## check if some circrna is still missed also after coordinate extension search
        still.missed <- fixed.match[is.na(N)]

        if(nrow(still.missed) > 0){
            ## check if the still missed have a different strand in reads
            still.missed$V6 <- sapply(still.missed$V6, switch.strand)

            ## expand interval of still missed
            extended.still.missed <-
                still.missed[, .(start = V2 + ext.range,
                                 end = V3 + ext.range),
                             by = .(V1, fixed.V2 = V2, fixed.V3= V3, V6)]

            still.missed.strand.switched <-
                merge(bks.reads,
                      extended.still.missed,
                      by.x = c("V1", "V2", "V3", "V6"),
                      by.y = c("V1", "start", "end",
                               "V6"))[, .(V1, V2 = fixed.V2,
                                          V3 = fixed.V3, V4,
                                          V5, V6)]

            ## missed circrnas fixed by strand switch
            if(nrow(still.missed.strand.switched[, .(ccp = .N), by = .(V1, V2, V3, V6)]) > 0){
                ## store fixed coordinates and strand for reads
                chimout.slices$fix.coo.str <-
                    merge(bks.reads,
                                    extended.still.missed,
                                    by.x = c("V1", "V2", "V3", "V6"),
                                    by.y = c("V1", "start", "end",
                                             "V6"))[, .(V1, V2 = fixed.V2,
                                                        V3 = fixed.V3, V4, V5,
                                                        V6 = sapply(V6, switch.strand))]
            }
        }
    }

    ## bind all reads (all.seg.bks.reads)
    # annotated.bks.reads <-
    #     rbindlist(list(unfixed.bks.reads,
    #                    strand.switched,
    #                    fixed.bks.reads,
    #                    still.missed.strand.switched),
    #               use.names = T)

    ## check the number of identified circrnas
    # nrow(orig.est) - nrow(annotated.bks.reads[, .(ccp = .N), by = .(V1, V2, V3, V6)])

    # chimout.slices <-
    #     list(unf = unfixed.bks.reads,
    #          unf.str.swi = strand.switched[, .(V1, V2, V3, V4, V5,
    #                                            V6 = sapply(V6, switch.strand))],
    #          fix.coo = merge(bks.reads,
    #                          unfix.orig.est,
    #                          by.x = c("V1", "V2", "V3", "V6"),
    #                          by.y = c("V1", "start", "end",
    #                                   "V6"),
    #                          all = F)[, .(V1, V2 = fixed.V2,
    #                                       V3 = fixed.V3, V4,
    #                                       V5, V6)],
    #          fix.coo.str = merge(bks.reads,
    #                              extended.still.missed,
    #                              by.x = c("V1", "V2", "V3", "V6"),
    #                              by.y = c("V1", "start", "end",
    #                                       "V6"))[, .(V1, V2 = fixed.V2,
    #                                                  V3 = fixed.V3, V4, V5,
    #                                                  V6 = sapply(V6, switch.strand))]
    #     )

    annotated.chimout.junc <-
        rbindlist(chimout.slices,
                  use.names = T,
                  idcol = "Set")
}else{

    annotated.chimout.junc <- bks.reads[, .(V1, V2, V3, V4, V6)]
}

splitted.filename <- strsplit(arguments$output, ".", fixed = T)[[1]]
if(tail(splitted.filename, 1) == "gz"){
    tmp.outfile <- sub(".gz$", "", arguments$output)
}

fwrite(x = annotated.chimout.junc[, .(V1, V2, V3, V4, V5, V6)],
       file = tmp.outfile,
       sep = "\t",
       col.names = F,
       row.names = F)

if(tail(splitted.filename, 1) == "gz"){
    gzip(tmp.outfile, destname = arguments$output)
}
