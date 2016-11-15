#!/usr/bin/env Rscript

#prepacks <- c("BH", "svglite")
#
#install.packages(prepacks, repos="http://cran.r-project.org", dependencies = T)

packs <- c("BH", "svglite", "optparse", "plyr", 
	   "scales", "reshape2", 
           "ggplot2", "ggthemes",
           "RSvgDevice")

install.packages(packs, repos="http://cran.r-project.org", dependencies = T)

#for(pack in packs){
#  install.packages(pack, repos="http://cran.r-project.org", dependencies = T)
#}

source("http://bioconductor.org/biocLite.R")
biocLite(pkgs = c("DESeq2", "ReportingTools"))
