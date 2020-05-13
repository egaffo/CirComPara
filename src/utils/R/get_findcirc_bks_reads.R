#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-s", "--sites_reads"), action = "store", type = "character",
                default = "sites.reads",
                help = "Read file (in FASTA format) as the ouput by -R option of find_circ.py (default: sites.reads)"),
    make_option(c("-c", "--circ_candidates"), action = "store", type = "character",
                default = "circ_candidates.bed",
                help = "BED file with backsplices, as find_circ.py main output (default: circ_candidates.bed)"),
    make_option(c("-r", "--reads_bed"), action = "store", type = "character",
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

sites.reads.file <- arguments$sites_reads
circ.candidates.bed.file <- arguments$circ_candidates
reads.bed <- arguments$reads_bed
read.list <- arguments$read_list

# sites.reads.file <- "/blackhole/enrico/circular/circompara_testing/circompara/test_circompara/analysis/samples/sample_A/processings/circRNAs/find_circ_out/sites.reads"
# circ.candidates.bed.file <- "/blackhole/enrico/circular/circompara_testing/circompara/test_circompara/analysis/samples/sample_A/processings/circRNAs/find_circ_out/circ_candidates.bed"

sites.reads <-
    fread(cmd = paste0('grep ">" ', sites.reads.file),
          sep = " ",
          header = F)[, .(V4 = sub(">", "", V1),
                          read_id = sub("(.*)_([0-9]+)", "\\1", V2))]

circ.candidates.bed <- fread(input = circ.candidates.bed.file,
                             header = F, select = c(1, 2, 3, 4, 5, 6))

bks.reads <- data.table()
all.reads <- data.table()
if(nrow(circ.candidates.bed) > 0 & nrow(sites.reads) > 0){
    bks.reads <- merge(sites.reads, circ.candidates.bed,
                       by = "V4", all.x = F,
                       all.y = T)[, .(V1, V2, V3, read_id, V5, V6)]

    all.reads <-
        bks.reads[, .N, by = read_id][order(-N), .(N, read_id)]
}

## write gzipped file for circular reads
reads_output.gz <- gzfile(reads.bed, "w")
write.table(bks.reads, file = reads_output.gz,
            sep = "\t", col.names = F, row.names = F, quote = F)
close(reads_output.gz)

write.table(all.reads, file = read.list,
            sep = "\t", col.names = F, row.names = F, quote = F)
