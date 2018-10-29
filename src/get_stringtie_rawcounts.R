#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
    make_option(c("-g", "--transcripts_gtf"), action = "store", type = "character", 
                help = "The _transcripts.gtf file output by StringTie"),
    make_option(c("-f", "--fastqc_data"), action = "store", type = "character", 
                help = paste0("The fastqc_data.txt output by FastQC, or a comma ",
                              "separated list of fastqc_data.txt files if paired-end reads. ",
                              "If a number is given then it will be interpreted as the ",
                              "(already computed) average read length"),
                default = "100"),
    make_option(c("-o", "--outprefix"), action = "store", type = "character", 
                help = "The prefix of output file names", default = "")
                
)

parser <- OptionParser(usage = "%prog -g sample_transcripts.gtf -f sample_R1/fastqc_data.txt,sample_R2/fastqc_data.txt", 
                       description = paste0("Computes the raw read counts for both transcripts and genes ",
                                            "like the prepDE.py script from StringTie. ",
                                            "It will output two files: the transcripts' and the genes' ", 
                                            "raw read counts."),
                       option_list = option_list)
opt <- parse_args(object = parser, print_help_and_exit = T, 
                        positional_arguments = F, 
                        convert_hyphens_to_underscores = F)
 
## This function parses a fastqc_data.txt file and computes the average read length
## TODO: improve by parsing "fastqc_data.txt" and get the mode instead of average read length
get.avg.read.len <- function(fastqc_data.txt.file){
    fastqc_data.txt <- readLines(fastqc_data.txt.file)
    mean(as.numeric(strsplit(strsplit(grep("Sequence length", 
                                           x = fastqc_data.txt, 
                                           value = T), 
                                      "\t")[[1]][2], "-")[[1]]))
}

sample_transcripts.gtf.file <- opt$transcripts_gtf
fastqc_data.txt.files <- opt$fastqc_data

## as in prepDE.py script describeed in StringTie manual, compute raw read counts from
## the 'cov' field. We need the average read length
if(suppressWarnings(is.na(as.numeric(fastqc_data.txt.files)))){
    avg_read_len <- mean(sapply(strsplit(fastqc_data.txt.files, ",")[[1]], get.avg.read.len))
}else{
    avg_read_len <- as.numeric(fastqc_data.txt.files)
}

## give a message on the average read length used
message(paste0("Average read length = ", avg_read_len))

sample_transcripts.gtf <- fread(sample_transcripts.gtf.file, skip = 2)

## get cov, transcript and gene id for each exon and compute exon length
exon.cov <- sample_transcripts.gtf[V3 == "exon", 
                                   .(gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
                                     transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', V9),
                                     exon_number = sub('.*exon_number "([^"]+)".*', '\\1', V9),
                                     len = V5 - V4 + 1,
                                     cov = sub('.*cov "([^"]+)".*', '\\1', V9))]

## transcript raw count is computed by multiplying the cov (= average read count 
## for each nucleotide) by the transcript length (sum of exons' length) and 
## divide by average read length
transcript.cov <- 
    exon.cov[, .(len = sum(len), cov = sum(as.numeric(cov))), 
             by = .(gene_id, transcript_id)][, raw.reads := ceiling(cov*len/avg_read_len), 
                                             by = .(gene_id, transcript_id)][]

write.csv(x = transcript.cov, 
          file = paste0(opt$outprefix, "transcript_expression_rawcounts.csv"), 
          row.names = F)

## gene raw count is computed as the sum of its transcripts' raw counts
## (prepDE.py seems not to merge overlapping exons)
gene.cov <- transcript.cov[, .(raw.reads = sum(raw.reads)), 
                           by = .(gene_id)]

write.csv(x = gene.cov, 
          file = paste0(opt$outprefix, "gene_expression_rawcounts.csv"), 
          row.names = F)
