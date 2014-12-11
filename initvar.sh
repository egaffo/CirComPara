#!/bin/bash

export SEGEMEHL_HOME="/blackhole/enrico/circular/tools/segemehl_0_1_9/segemehl/"
export SEGEMEHL_INDEX="/blackhole/enrico/circular/tools/segemehl_indexes/hg38/hg38.idx"
export GENOME_FASTA="/blackhole/enrico/circular/tools/segemehl_indexes/hg38/chr*.fa"
export SRATOOLKIT_HOME="/home/enrico/tools/sra/sratoolkit.2.3.5-2-ubuntu64/bin/"
export SAMTOOLS_HOME="/blackhole/enrico/circular/tools/samtools-1.1/"
export FASTX_TOOLKIT="$HOME/tools/fastx_toolkit/bin"
export TRIMMOMATIC="/home/enrico/tools/trimmomatic/trimmomatic-0.32.jar"
export JUNK2_HOME="/home/enrico/tools/junk2"
export CUFFLINKS_HOME="/home/enrico/tools/cufflinks"

export PATH="$CUFFLINKS_HOME:$JUNK2_HOME/bin:$FASTX_TOOLKIT:$SEGEMEHL_HOME:$SRATOOLKIT_HOME:$SAMTOOLS_HOME:$PATH"
