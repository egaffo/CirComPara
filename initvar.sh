#!/bin/bash

export JUNK2_HOME="/home/enrico/tools/junk2"

export TRIMMOMATIC="/home/enrico/tools/trimmomatic/trimmomatic-0.32.jar"

SEGEMEHL_HOME="/home/enrico/tools/segemehl/"
SRATOOLKIT_HOME="/home/enrico/tools/sra/bin/"
SAMTOOLS_HOME="/home/enrico/tools/samtools/"
FASTX_TOOLKIT="$HOME/tools/fastx_toolkit/bin"
CUFFLINKS_HOME="/home/enrico/tools/cufflinks"

export PATH="$CUFFLINKS_HOME:$JUNK2_HOME/bin:$FASTX_TOOLKIT:$SEGEMEHL_HOME:$SRATOOLKIT_HOME:$SAMTOOLS_HOME:$PATH"
