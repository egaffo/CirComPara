#!/bin/bash

SEGEMEHL_HOME="/home/enrico/tools/segemehl/"
#export SEGEMEHL_INDEX="/blackhole/enrico/circular/tools/segemehl_indexes/hg38/hg38.idx"
#export GENOME_FASTA="/blackhole/enrico/circular/tools/segemehl_indexes/hg38/chr*.fa"
SRATOOLKIT_HOME="/home/enrico/tools/sra/bin/"
SAMTOOLS_HOME="/home/enrico/tools/samtools/"
FASTX_TOOLKIT="$HOME/tools/fastx_toolkit/bin"
export TRIMMOMATIC="/home/enrico/tools/trimmomatic/trimmomatic-0.32.jar"
export JUNK2_HOME="/home/enrico/tools/junk2"
CUFFLINKS_HOME="/home/enrico/tools/cufflinks"

export PATH="$CUFFLINKS_HOME:$JUNK2_HOME/bin:$FASTX_TOOLKIT:$SEGEMEHL_HOME:$SRATOOLKIT_HOME:$SAMTOOLS_HOME:$PATH"
