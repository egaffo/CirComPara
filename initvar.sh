#!/bin/bash

export JUNK2_HOME="/home/enrico/tools/junk2"

export TRIMMOMATIC="/home/enrico/tools/trimmomatic/trimmomatic-0.32.jar"

SEGEMEHL_HOME="/home/enrico/tools/segemehl/"
SRATOOLKIT_HOME="/home/enrico/tools/sra/bin/"
SAMTOOLS_HOME="/home/enrico/tools/samtools/"
FASTX_TOOLKIT="$HOME/tools/fastx_toolkit/bin"
CUFFLINKS_HOME="/home/enrico/tools/cufflinks"

export PATH="$CUFFLINKS_HOME:$JUNK2_HOME/bin:$FASTX_TOOLKIT:$SEGEMEHL_HOME:$SRATOOLKIT_HOME:$SAMTOOLS_HOME:$PATH"

## other tools to add to PATH for circpipe
CIRC_TOOLS=/blackhole/circrna/tools

export PATH=$CIRC_TOOLS/samtools:$CIRC_TOOLS/tophat:$CIRC_TOOLS/segemehl:$CIRC_TOOLS/bwa.kit:$CIRC_TOOLS/ciri:$CIRC_TOOLS/bowtie2:$CIRC_TOOLS/find_circ:$CIRC_TOOLS/bedtools2:$PATH
export JUNK2_HOME=$CIRC_TOOLS/junk2

## add CIRCExplorer Pyhton library installation path to PYTHONPATH
export PYTHONPATH="$PYTHONPATH:$CIRC_TOOLS/CIRCexplorer/lib/python2.7/site-packages/"

TOOLS_DIR=/blackhole/circrna/tools
TOPHAT_HOME=$TOOLS_DIR/tophat2
HISAT2_HOME=$TOOLS_DIR/hisat2-2.0.1-beta
STAR=$TOOLS_DIR/star-2.5/STAR/bin/Linux_x86_64
CIRCEXPLORER=$TOOLS_DIR/CIRCexplorer/bin

export PATH=$TOPHAT_HOME:$HISAT2_HOME:$JUNK2_HOME/bin:$STAR:$CIRCEXPLORER:$PATH
