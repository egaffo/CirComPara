#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bedr))

option_list <- list(
    make_option(c("-c", "--combined_circrnas"), action = "store", type = "character",
                help="combined_circrnas.gtf.gz"),
    make_option(c("-d", "--cluster_dist"), action = "store", type = "integer",
                default = 5000L,
                help="The maximum distance allowed to cluster intergenic circRNAs"),
    make_option(c("-o", "--outdir"), action = "store", type = "character",
                default = "./",
                help="Output directory were circ_to_genes.tsv and gene_to_circ.tsv will be saved")
)

parser <- OptionParser(usage="%prog -c circular_expression/circRNA_collection/combined_circrnas.gtf.gz -d 5000 -o ../circular_expression/circRNA_collection",
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

gene.annotation.file <- arguments$combined_circrnas #"../circular_expression/circRNA_collection/combined_circrnas.gtf.gz"
max.dist <- arguments$cluster_dist
## prepare result dir
results.path <- arguments$outdir

dir.create(path = results.path, recursive = T, showWarnings = F)

paste.bar <- function(x){
    sapply(x, function(y){
        paste0(sort(y), collapse = "|")}
    )
}

simplify.gene.region <- function(x){
    label <- "unknown"
    if(grepl("exon", x)){
        label <- "exon"
    }else{
        if(grepl("intron", x)){
            label <- "intron"
        }else{
            if(grepl("intergenic", x)){
                label <- "intergenic"
            }
        }
    }
}

simplify.gene.region.list <- function(x){
    label <- "unknown"
    if("exon" %in% x){
        label <- "exon"
    }else{
        if("intron" %in% x){
            label <- "intron"
        }else{
            if("intergenic" %in% x){
                label <- "intergenic"
            }
        }
    }
}

## function to cluster intergenic circrnas
## intergenic.circs.bed is a data.table in BED format
get.circ.cluster <- function(intergenic.circs.bed, max.dist = 5000L){
    fw.strand.interg.circ <-
        intergenic.circs.bed[, paste0(chr, ":", start, "-", end)]
    circ.clusters <- data.table(cluster.region(fw.strand.interg.circ,
                                               distance = max.dist,
                                               check.chr = F, verbose = F))

    circ.clusters[, .(circ_id = index, cluster = regionIndex)]
}

## get gene annotation
gene.annotation <- fread(gene.annotation.file)

gene.annotation[, c("gene_id", "gene_name",
                    "gene_biotype"):=tstrsplit(V18, ";", fixed = T, fill = ".")]
gene.annotation[V12 == ".", V12 := "intergenic"]

gene.annotation[, `:=`(circ_id = sub('.*gene_id "([^"]+)".*', "\\1", V9),
                       gene_id = sub('.*gene_id "([^"]+)".*', "\\1", gene_id),
                       gene_name = sub('.*gene_name "([^"]+)".*', "\\1", gene_name),
                       gene_biotype = sub('.*gene_biotype "([^"]+)".*', "\\1", gene_biotype))]

## collapse/remove circ_id strand
gene.annotation[, circ_id := sub(":[+-]$", "", circ_id)]

circ.id.strand <- gene.annotation[, .(strand = paste0(unique(V7), collapse = "|")),
                                  by = .(circ_id)]

## CIRCRNA TO GENES
## make a table with all circrna ids in one column, one per row, and
## the corresponding genes in a side column, in a bar-separated list.
### first, compose the circ -> gene-list table
circ_to_genelist <-
    unique(gene.annotation[, .(circ_id, gene_id, gene_name,
                               gene_biotype,
                               gene_region = V12,
                               gene_strand = V16)])[, .(gene_ids = list(unique(gene_id)),
                                                        gene_names = list(unique(gene_name)),
                                                        gene_biotypes = list(unique(gene_biotype)),
                                                        gene_region = list(unique(gene_region)),
                                                        gene_strand = list(unique(gene_strand))),
                                                    by = circ_id]
circ_to_genelist[, simple_gene_region := sapply(gene_region, simplify.gene.region.list)]

## Intergenic clusters:
## save intergenic circrnas in a BED file to be input to bedtools cluster
## which will tell how many cluster (putative new genes) express the circrnas
intergenic.circs <- circ_to_genelist[simple_gene_region == "intergenic", .(circ_id)]

## convert into BED format
intergenic.circs.bed <-
    intergenic.circs[, c("chr", "coords") :=
                         tstrsplit(circ_id,
                                   ":")][, c("start",
                                             "end") :=
                                             tstrsplit(coords, "-",
                                                       type.convert = T)][, .(chr,
                                                                              start,
                                                                              end,
                                                                              circ_id)]

## cluster circrnas irrespective of strand
intergenic.circ.clusters <-
    merge(intergenic.circs.bed,
          get.circ.cluster(intergenic.circs.bed, max.dist = max.dist),
          by = "circ_id")
intergenic.circ.clusters[, cluster.coords := paste0("CC", cluster, "_",
                                                    chr, ":", min(start),
                                                    "-", max(end)),
                         by = cluster]

## set gene_id as cluster name for intergenic circrnas
circ_to_genelist <-
    merge(circ_to_genelist, intergenic.circ.clusters[, .(circ_id, cluster.coords)],
          by = "circ_id",
          all.x = T)[!is.na(cluster.coords),
                     gene_names := as.list(cluster.coords)][, .(circ_id, gene_ids,
                                                                gene_names, gene_biotypes,
                                                                gene_region,
                                                                simple_gene_region,
                                                                gene_strand)]

### now, convert the gene lists to characters by separating genes with a bar char
circ_to_genes <-
    circ_to_genelist[, lapply(.SD, paste.bar), by = circ_id]

## remove intergenic notation when a known host genes is listed
remove.point <- function(x)sub("\\|\\.|\\.\\|", "", x)

cols <- c("gene_ids", "gene_names", "gene_biotypes", "gene_strand")
circ_to_genes[simple_gene_region != "intergenic",
                  (cols) := lapply(.SD, remove.point), #by = circ_id,
                  .SDcols = cols]

## attach circ_id strand
circ_to_genes <- merge(circ.id.strand, circ_to_genes, by = "circ_id")
## fix circrna strand by using host gene strand
circ_to_genes[strand == "-|+", strand := gene_strand]
circ.id.strand <- circ_to_genes[, .(circ_id, strand)]

## save the table
fwrite(x = circ_to_genes,
       file = file.path(results.path, "circ_to_genes.tsv"),
       quote = F, sep = "\t", row.names = F)

## GENES TO CIRCRNAS
## make a table with all genes in one column and corresponding circrnas in a
## side column, one per row, indeed duplicating gene ids if necessary
gene_to_circ <-
    unique(gene.annotation[, .(circ_id, gene_id,
                               gene_name, gene_biotype,
                               gene_region = V12,
                               gene_strand = V16)])[, .(gene_region = list(gene_region)),
                                                    by = .(circ_id, gene_id, gene_name,
                                                           gene_biotype, gene_strand)]

gene_to_circ <-
    merge(gene_to_circ,
          intergenic.circ.clusters[, .(circ_id, cluster.coords)],
          by = "circ_id", all.x = T)[!is.na(cluster.coords),
                                     gene_name := cluster.coords][, .(circ_id, gene_id,
                                                               gene_name, gene_biotype,
                                                               gene_region, gene_strand)]

## remove rows that were (erroneously) kept because only one end of the
## backsplice was intergenic
gene_to_circ <- gene_to_circ[gene_name != "."]

## set simple region labels
gene_to_circ[, simple_gene_region := sapply(gene_region, simplify.gene.region.list)]

## convert lists into strings
gene_to_circ <-
    gene_to_circ[, .(gene_region = paste.bar(gene_region)),
                 by = .(gene_id, gene_name, circ_id,
                        simple_gene_region, gene_biotype, gene_strand)]

## attach circrna strand
gene_to_circ <- merge(circ.id.strand,
                      gene_to_circ,
                      by = "circ_id")

## save table
fwrite(x = gene_to_circ,
            file = file.path(results.path, "gene_to_circ.tsv"),
            quote = F, sep = "\t", row.names = F)
