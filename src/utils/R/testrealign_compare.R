#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))

option_list <- list(
  make_option(c("-l", "--labels"), action="store", type="character",
              help="A comma separated list of labels, one for each testrealign.x output given as input"),
  make_option(c("-i", "--input"), action="store", type="character",
              help="A comma separated list of testrealign.x outputs. List must be ordered according to label list"),
  make_option(c("-o", "--output"), action="store", type="character",
              help="The output file (CSV format)"),
  make_option(c("-r", "--minreads"), action="store", type="integer", default=2,
              help="TODO: Minimum number of reads to consider a valid circular junction"),
  make_option(c("-f", "--fixstart"), action="store", type="logical", default=T,
	      help="Decrease by 1 position the start coordinates: this fixes a possible testrealign bug")
)

parser <- OptionParser(usage="%prog -l S1,S2,S3 -i S1/splicesites.bed,S2/splicesites.bed,S3/splicesites.bed -o testrealign_compared.csv ",
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)
opt <- arguments$options
minreads <- opt$minreads

if(length(arguments$input)>0){

  labels <- unlist(strsplit(arguments$labels, ',', fixed=T))
  inputs <- unlist(strsplit(arguments$input,  ',', fixed=T))
  names(inputs) <- labels
  # colnames <- c("CHR", "POS1", "POS2", "INFO", "SCORE", "STRAND")
  colnames <- c("CHR", "POS1", "POS2", "N", "QUAL", "STRAND", "SCORE")
  splicesites.colclasses <- c("factor", "numeric", "numeric", "numeric",
                              "numeric", "character", "numeric")
  combined.df <- ldply(inputs, function(x){
    a <- read.table(file = x, sep="\t", header = F, skip = 1,
                    col.names = colnames, colClasses = splicesites.colclasses,
                    stringsAsFactors = F); a}
  )


  colnames(combined.df)[1] <- "sampleID"

  if(arguments$fixstart){
	combined.df$POS1 <- combined.df$POS1-1
  }

  #   The information INFO in field 4 shows
  #   three numbers. The first number indicates the number of reads that support the
  #   reported splice site. The second and the third number inform on the total number
  #   of splits that can be observed at POS1 and POS2, respectively.
  #   The letter ’N’ stands for ’normal’, indicating that this is a regular splice junction.
  #   The letter ’C’ stands for ’circular’, indicating that this is a potential backsplice junction.
  #   The last letter is a quality flag ’P’ stands for ’passed’ and indicates that there was only one cor-
  #   responding acceptor-donor pair.
  #   The letter ’M’ indicates that their are ’multiple’ acceptor-donor pairs.
  #   The letter ’F’ is a warning. It indicates that the algorithm
  #   identified sequence similarities between different acceptor or donor sites, i.e. the
  #   junction could be wrong.
  #   The strand information is only meaningful if a strand specific sequencing protocol was used

  # splicesites.info <- data.frame(do.call('rbind', strsplit(as.character(combined.df$INFO), ':', fixed=TRUE)), stringsAsFactors =F)
  # colnames(splicesites.info) <- c("SPLITS", "READS", "READS_POS1", "READS_POS2", "SPLICE_TYPE", "QUALITY")
  # splicesites.info$READS <- as.integer(splicesites.info$READS)
  # splicesites.info$READS_POS1 <- as.integer(splicesites.info$READS_POS1)
  # splicesites.info$READS_POS2 <- as.integer(splicesites.info$READS_POS2)
  # splicesites.ext <- cbind(combined.df, splicesites.info)
  # splicesites.ext$SCORE <- splicesites.ext$READS
  #
  # #passed.circular <- subset(splicesites.ext, SPLICE_TYPE=="C" & QUALITY=="P" & READS>=minreads)
  # passed.circular <- subset(splicesites.ext, SPLICE_TYPE=="C" & QUALITY=="P", select = c(1, 2, 3, 4, 6, 7, 9, 10, 11))
  passed.circular <- combined.df[, c("sampleID", "CHR", "POS1", "POS2", "QUAL", "STRAND", "SCORE")]

  if(length(arguments$output)>0){
    write.csv(passed.circular, file = arguments$output, row.names = F)
  }else{
    write.csv(passed.circular, file = stdout(), row.names = F)
  }
}else{
  cat("Please give a valid input table file\n")
  print_help(parser)
  stop()
}
