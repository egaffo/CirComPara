#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character",
              help="A comma separated table of raw read counts for different samples as given by DESeqDataSetFromHTSeqCount()"),
  make_option(c("-s", "--samples"), action="store", type="character",
              help="A comma separated table, describing samples and conditions"),
  make_option(c("-o", "--output"), action="store", type="character",
              help="The DESeq2 differential expression tests file (CSV format)"),
  make_option(c("-d", "--reportDirectory"), action="store", type="character",
              help="The directory where to store results' report")
)

parser <- OptionParser(usage="%prog -i rawReads.csv -o differentialExpressionTests.csv -s sampleTable.csv -d DESeq2_HTMLReport",
                       option_list=option_list)
arguments <- parse_args(parser, positional_arguments=F)

if(length(arguments$input) > 0){

  sampleTable <- read.csv(file = arguments$sample)
  rownames(sampleTable) <- as.character(sampleTable[, 1])
  countData   <- read.csv(file = arguments$input, row.names = 1)

  dds <- tryCatch({
      suppressPackageStartupMessages(library(DESeq2))

      ddsHTSeq <- DESeqDataSetFromMatrix(countData = countData,
                                         colData = sampleTable,
                                         design = ~ condition)
      DESeq(ddsHTSeq)


  }, error = function(e){
      writeLines(paste("DESeq failed with the following error:", e, sep = "\n"),
                 arguments$output)
      NULL
  })

  if(!is.null(dds)){
      res <- results(dds)
      write.csv(as.data.frame(res), file = arguments$output)
      pvalueCutoff <- 0.1
      if(any(res$padj <= pvalueCutoff, na.rm = T)){
          ## if no differentially expressed elements the following would fail
          ## result reporting tool
          suppressPackageStartupMessages(library(ReportingTools))

          reportDirectory <- arguments$reportDirectory
          des2Report <- HTMLReport(shortName = "RNAseq_analysis_with_DESeq2",
                                   title = "RNA-seq analysis of differential expression using DESeq2",
                                   reportDirectory = reportDirectory)
          publish(dds, des2Report, pvalueCutoff = pvalueCutoff,
                  annotation.db = "", factor = colData(dds)$condition,
                  reportDir = reportDirectory)
          finish(des2Report)
      }else{
          writeChar(paste0("No differentially expressed elements at ", pvalueCutoff, " FDR"),
                    file.path(arguments$reportDirectory, "RNAseq_analysis_with_DESeq2.html"))
      }

      writeLines(capture.output(sessionInfo()),
                 file.path(dirname(arguments$output), "DESeq2_diffexp.sessionInfo"))
  }
}else{
    cat("Please give a valid input table file\n")
    print_help(parser)
    stop()
}
