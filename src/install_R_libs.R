#!/usr/bin/env Rscript

packs <- c("optparse", 
           "ggplot2", "ggthemes",
           "scales", "reshape2", 
           "data.table", "plyr", 
           "RSvgDevice")

for(pack in packs){
  install.packages(pack, repos="http://cran.r-project.org")
}

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("ReportingTools")