#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-c", "--circrnas_gtf"), action = "store", type = "character",
                help = "A circrnas.gtf file or a text file listing circrnas.gtf file paths to merge"),
    make_option(c("-r", "--min_reads"), action = "store", type = "integer",
                default = 2,
                help = "The minimum detection read threshold for a circRNA  (in at least one sample)"),
    make_option(c("-m", "--min_methods"), action = "store", type = "integer",
                default = 2,
                help = "Keep circRNAs commonly detected by >= m circRNA detection methods (in at least one sample)"),
    make_option(c("-l", "--lin_reads"), action = "store", type = "character",
                default = "",
                help = "The backsplice linear expression read count file or text file listing bks_linear_counts.tab.gz file paths to merge"),
    make_option(c("-o", "--outdir"), action = "store", type = "character",
                default = "./",
                help = "Output directory")
)

parser <- OptionParser(usage = paste0("%prog -c analysis_A/circular_expression/circRNA_collection/circrnas.gtf,",
                                      "analysis_B/circular_expression/circRNA_collection/circrnas.gtf -r 2 -m 2 ",
                                      "-o merged/circRNA_expression"),
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments = F)

circrnas.gtf.files <- arguments$circrnas_gtf
min_reads <- arguments$min_reads
min_methods <- arguments$min_methods
lin_reads <- arguments$lin_reads

## prepare result dir
results.path <- arguments$outdir
dir.create(path = results.path, recursive = T, showWarnings = F)

if(ncol(fread(circrnas.gtf.files, showProgress = F, nrows = 1)) == 1){
    ## case: text file listing circrnas.gtf files
    circrnas.gtf.files <- readLines(circrnas.gtf.files)
}

## read circRNA results: filter low expressed (less than min_reads reads) circRNAs
colClasses <- c("factor", "factor", "character", "integer",
                "integer", "integer", "factor", "character", "character")
circrnas.gtf <- rbindlist(lapply(circrnas.gtf.files, fread, data.table = T,
                                 colClasses = colClasses, showProgress = F),
                          use.names = T)


circrnas.gtf[, `:=` (sample_id = sub('.*sample_id "([^"]*)".*', "\\1", V9),
                     gene_id = paste0(V1, ":", V4, "-", V5))
             ][, V9 := NULL]

## sumup expression by circ_id (within the same detection method and sample)
circrnas.gtf <-
    circrnas.gtf[, .(V6 = sum(V6)),
                 by = .(V1, V2, V4, V5, sample_id, gene_id)]

## save
fwrite(x = circrnas.gtf[, .(sample_id, circ_id = gene_id, chr = V1,
                            start = V4, end = V5, method = V2,
                            reads = V6)][order(sample_id, circ_id, method)],
       file = file.path(results.path, "unfiltered_circrnas.csv"),
       showProgress = F)

# ## circrna ids detected with >= 'min_reads' reads
# xpr.filtered.circ.ids <- circrnas.gtf[V6 >= min_reads, unique(gene_id)]

## compute median read count, counting also methods that detected the circRNA
## with less than 'min_reads' reads
circrna.median.reads <- circrnas.gtf[, .(median.reads = median(V6),
                                         n_methods = .N),
                                     by = .(sample_id, circ_id = gene_id)]

## >= 'min_methods' methods must have detected the circRNA with >= 'min_reads'
meth.filtered.circ.ids <-
    circrnas.gtf[V6 >= min_reads,
                 .(n_methods = .N),
                 by = .(sample_id, gene_id)][n_methods >= min_methods,
                                             unique(gene_id)]

## make wide expression table circ ~ sample
circrna.median.reads.w <-
    dcast(data = circrna.median.reads, #fun.aggregate = sum,
          formula = circ_id ~ sample_id,
          value.var = "median.reads",
          fill = 0)
## attach column reporting if the circrna is reliable
circrna.median.reads.w[, reliable.circrna := 0][circ_id %in% meth.filtered.circ.ids,
                                                reliable.circrna := 1]

## save
sample_ids <-
    colnames(circrna.median.reads.w)[!colnames(circrna.median.reads.w) %in%
                                         c("circ_id", "reliable.circrna")]

fwrite(x = circrna.median.reads.w[order(-reliable.circrna, circ_id),
                                  c("circ_id", "reliable.circrna", sample_ids),
                                  with = F],
       file = file.path(results.path, "ccp_circrna_raw_xpr.csv"),
       showProgress = F)

## make wide n_methods table circ ~ sample
circrna.n.methods.w <-
    dcast(data = circrna.median.reads, #fun.aggregate = sum,
          formula = circ_id ~ sample_id,
          value.var = "n_methods",
          fill = 0)

circrna.n.methods.w <-
    merge(circrnas.gtf[, .N,
                       by = .(circ_id = gene_id,
                              V2)][, .(overall.n.methods = .N),
                                   by = .(circ_id)],
          circrna.n.methods.w,
          by = "circ_id")

## attach column reporting if the circrna is reliable
circrna.n.methods.w[, reliable.circrna := 0][circ_id %in% meth.filtered.circ.ids,
                                             reliable.circrna := 1]

circrna.n.methods.w <-
    circrna.n.methods.w[, c("circ_id", "reliable.circrna",
                            "overall.n.methods", sample_ids),
                        with = F]

## save
fwrite(x = circrna.n.methods.w[order(-reliable.circrna, -overall.n.methods, circ_id)],
       file = file.path(results.path, "ccp_circrna_n_methods.csv"),
       showProgress = F)

## save which methods detected the circrna in each sample and overall
circ.meth.per.sample <-
    circrnas.gtf[, .(circ.methods = paste(sort(unique(V2)), collapse = "|")),
                 by = .(sample_id, circ_id = gene_id)]

overall.circ.meth <-
    circrnas.gtf[, .N,
                 by = .(circ_id = gene_id,
                        V2)][, .(overall.circ.methods = paste(sort(unique(V2)),
                                                              collapse = "|")),
                             by = .(circ_id)]

circ.meth.tab <-
    merge(overall.circ.meth,
          dcast(circ.meth.per.sample,
                formula = circ_id ~ sample_id,
                value.var = "circ.methods",
                fill = ""),
          by = c("circ_id"))

## attach column reporting if the circrna is reliable
circ.meth.tab[, reliable.circrna := 0][circ_id %in% meth.filtered.circ.ids,
                                       reliable.circrna := 1]

circ.meth.tab <-
    circ.meth.tab[, c("circ_id", "reliable.circrna",
                      "overall.circ.methods", sample_ids),
                  with = F][circrna.n.methods.w[order(-reliable.circrna,
                                                      -overall.n.methods, circ_id),
                                                .(circ_id)]]

## save
fwrite(x = circ.meth.tab,
       file = file.path(results.path, "ccp_circrna_methods.csv"),
       quote = T,
       showProgress = F)

## save expression matrix for each circrna method
circ.meth.xpr.mats <-
    lapply(split(circrnas.gtf[, .(circ.method = V2, circ_id = gene_id, sample_id,
                              reads = V6)],
             by = "circ.method"),
       dcast, formula = circ_id ~ sample_id, value.var = "reads", fill = 0)

circ.meth.xpr.mats.files <-
    sapply(names(circ.meth.xpr.mats),
           function(x)fwrite(circ.meth.xpr.mats[[x]],
                             file = file.path(results.path, paste0(x, "_raw_xpr.csv"))))

if(lin_reads != ""){

    ## linear expression
    first.line.dt <- fread(lin_reads, showProgress = F, nrows = 1, header = F)
    if(ncol(first.line.dt) == 1){
        ## case: text file listing circrnas.gtf files
        lin_reads <- readLines(lin_reads)
    }

    first.line.dt <- fread(lin_reads[1], showProgress = F, nrows = 1, header = F)
    if(first.line.dt$V1 == "chr"){
        ## DCC input: samples in various tables need to be merged
        bks_linear_counts.tab.gz.tables <-
            lapply(lin_reads, fread, data.table = T,
                   showProgress = F, header = T)

        melt.dcc.linexp.tab <- function(x){
            ## DCC uses BED coordinates, we need to convert
            x[, circ_id := paste0(chr, ":", start+1, "-",
                                  end)][, `:=`(chr = NULL,
                                               start = NULL,
                                               end = NULL)]
            melt(x,
                 id.vars = "circ_id",
                 variable.name = "sample_id",
                 value.name = "lin.reads")
        }

        bks_linear_counts.tab.gz.m <-
            rbindlist(lapply(bks_linear_counts.tab.gz.tables,
                             melt.dcc.linexp.tab),
                      use.names = T)
    }else{
        ## CCP input
        bks_linear_counts.tab.gz <- rbindlist(lapply(lin_reads, fread, data.table = T,
                                                     showProgress = F),
                                              use.names = T)
        bks_linear_counts.tab.gz[, c("sample_id", "chr") := (tstrsplit(V1, ":"))]
        bks_linear_counts.tab.gz[, `:=`(sample_id = sub("_bks_linear_counts.tab$", "",
                                                        basename(sample_id)),
                                        circ_id = sub('.*"([^"]+)".*', "\\1", V9))]
        bks_linear_counts.tab.gz[, V1 := NULL]
        ## remove strand info from circ_id
        bks_linear_counts.tab.gz[, circ_id := sub(":[-+.]", "", circ_id)]

        bks_linear_counts.tab.gz.m <-
            bks_linear_counts.tab.gz[, .(lin.reads = sum(V10)),
                                     by = .(sample_id, circ_id)]
    }
    ## save
    bks_linear_counts.tab.gz.filename <- file.path(results.path, "ccp_bks_linexp.csv")
    fwrite(x = bks_linear_counts.tab.gz.m,
           file = bks_linear_counts.tab.gz.filename,
           showProgress = F)

    # library(R.utils)
    # gzip(bks_linear_counts.tab.gz.filename,
    #      destname = paste0(bks_linear_counts.tab.gz.filename, '.gz'))

}


