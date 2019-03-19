#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-c", "--counts"), action = "store", type = "character",
                help = "CircRNACounts file resulting from DCC"),
    make_option(c("-d", "--coordinates"), action = "store", type = "character",
                help = "CircCoordinates file resulting from DCC"),
    make_option(c("-o", "--output"), action = "store", type = "character",
                default = "strandedCircRNACounts", help = "The output file")
)

parser <- OptionParser(usage = "%prog -c CircRNACounts -d CircCoordinates -o strandedCircRNACounts",
                       option_list = option_list,
                       description = "Attach strand column to DCC's CircRNACount file")
arguments <- parse_args(parser, positional_arguments = F)
circcount.file <- arguments$counts
coordinates.file <- arguments$coordinates
output <- arguments$output

circrnacounts <- fread(file = circcount.file)
coordinates <- fread(file = coordinates.file,
                     col.names = c("Chr", "Start", "End", "info", "score", "Strand"))

if(! "Strand" %in% colnames(circrnacounts)){
    circrnacounts$Strand <- coordinates$Strand
}

write.table(x = circrnacounts[, .(Chr, Start, End, Strand, Chimeric.out.junction)],
            file = output, row.names = F, sep = "\t", quote = F)
