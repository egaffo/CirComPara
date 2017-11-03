---
Title:      CirComPara  
Subtitle:   a multi-method comparative bioinformatics pipeline to detect and study circRNAs from RNA-seq data  
Project:    CirComPara  
Author:     Enrico Gaffo  
Affiliation:    Compgen - University of Padova  
Web:        http://compgen.bio.unipd.it  
Date:       December 21, 2016  
output: 
  html_document: 
    keep_md: no
    number_sections: no
    toc: no
---

# CirComPara

CirComPara is a computational pipeline to detect, quantify, and correlate expression of linear and circular RNAs from RNA-seq data.

<!--TODO: more exhaustive description -->

## Quick install

Execute the following commands to download and install (locally) in your system the scripts and tools required to run CirComPara. 
If something goes wrong with the installation process try to manually install the software as described below.

Download and extract [the latest release of CirComPara][circompara_pack_link], or clone the GIT repository, enter CirComPara directory and run the automatic installer script:  

    git clone http://github.com/egaffo/CirComPara
    cd CirComPara
    ./install_circompara

### Test your installation

NB: in the `sed` string change the `/full/circompara/dir/path` path with your installation directory 

    cd test_circompara
    sed "s@\$CIRCOMPARA@/full/circompara/dir/path@g" vars.py > analysis/vars.py
    sed "s@\$CIRCOMPARA@/full/circompara/dir/path@g" meta.csv > analysis/meta.csv
    cd analysis
    ../../circompara

If you plan to use single-end reads, test with `meta_se.csv` file instead of `meta.csv`.  

If you receive some error messages try to follow instructions in **Installation troubleshooting** section.

### Add  CirComPara to your environment

Once completed the installation, if you do not want to type the whole path to the CirComPara executable each time, you can update your `PATH` environment variable. From the terminal type the following command (replace the `/path/to/circompara/install/dir` string with CirComPara's actual path)   

    export PATH=/path/to/circompara/install/dir:$PATH

Another way is to link CirComPara's main script in your local `bin` directory  

    cd /home/user/bin
    ln -s /path/to/circompara/install/dir/circompara_CirComPara

## CirComPara Docker image

A [Docker image of CirComPara](http://hub.docker.com/r/egaffo/circompara-docker/) is available from DockerHub.

To pull the image:

    docker pull egaffo/circompara-docker
    
You'll find the instructions on how to use the docker image at https://hub.docker.com/r/egaffo/circompara-docker.

# How to use

## Set your analysis project

This section shows how to set your project directory and run the analysis.
To run an analysis usually you want to specify your data (the sequenced reads in FASTQ format) and a reference genome in FASTA format.

### Compose META file

You have to specify read files, sample names and sample experimental condition in a metadata table file. The file format is a comma separated text file with the following header:  

    file,sample,condition
  
Then, each row corresponds to a read file. If you have paired-end sequenced samples write one line per file with the same sample name and condition. 

An example of the metadata table:

file|sample|condition
----|------|---------
/home/user/reads_S1_1.fq|S1|WT
/home/user/reads_S1_2.fq|S1|WT
/home/user/reads_S2_1.fq|S2|MU
/home/user/reads_S2_1.fq|S2|MU


and metadata file content:

    file,sample,condition
    /home/user/reads_S1_1.fq,S1,WT
    /home/user/reads_S1_2.fq,S1,WT
    /home/user/reads_S2_1.fq,S2,MU
    /home/user/reads_S2_1.fq,S2,MU

In the meta file you can also specify the adapter sequences to preprocess the reads, just add an `adapter` column with the adpter file.

file|sample|condition|adapter
----|------|---------|-------
/home/user/reads_S1_1.fq|S1|WT|/home/user/circompara/adapter.fa
/home/user/reads_S1_2.fq|S1|WT|/home/user/circompara/adapter.fa


### Specify the reference genome file

A required parameter is the reference genome. You can either pass the reference genome from the command line

    ./circompara "GENOME_FASTA='/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'"

or by setting the `GENOME_FASTA` parameter in the `vars.py` file; e.g.:

    GENOME_FASTA = '/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'


### Specify options in vars.py

Although parameters can be set from command line (sorrounded by quotes), you can set them in a local `vars.py` file, which must be placed in the analysis directory. Parameters not specified by the user will take defaulkt values.  
Below there is the full list of the parameters:

```
META: The metadata table file where you specify the project samples, etc.
    default: meta.csv

ANNOTATION: Gene annotation file (like Ensembl GTF/GFF)
    default: 

GENOME_FASTA: The FASTA file with the reference genome
    default: 

CIRCRNA_METHODS: Comma separated list of circRNA detection methods to use. Repeated values will be collapsed into unique values. Currently supported: ciri, find_circ, circexplorer2_star, circexplorer2_bwa, circexplorer2_tophat, circexplorer2_segemehl, testrealign (unfiltered segemehl; use of circexplorer2_segemehl is recommended for a better filtering of segemehl predictions). Set an empty string to use all methods available (including deprecated methods). 
    default: ciri,find_circ,circexplorer2_star,circexplorer2_bwa,circexplorer2_segemehl

CPUS: Set number of CPUs
    default: 4

GENEPRED: The genome annotation in GenePred format
    default: 

GENOME_INDEX: The index of the reference genome for HISAT2
    default: 

SEGEMEHL_INDEX: The .idx index for segemehl
    default: 

BWA_INDEX: The index of the reference genome for BWA
    default: 

BOWTIE2_INDEX: The index of the reference genome for BOWTIE2
    default: 

STAR_INDEX: The directory path where to find Star genome index
    default: 

BOWTIE_INDEX: The index of the reference genome for BOWTIE when using CIRCexplorer2_tophat
    default: 

HISAT2_EXTRA_PARAMS: Extra parameters to add to the HISAT2 aligner fixed parameters '--dta --dta-cufflinks --rg-id <SAMPLE> --no-discordant --no-mixed --no-overlap'. For instance, '--rna-strandness FR' if stranded reads are used.
    default: 

BWA_PARAMS: Extra parameters for BWA
    default: 

SEGEMEHL_PARAMS: SEGEMEHL extra parameters
    default: 

TOPHAT_PARAMS: Extra parameters to pass to TopHat
    default: 

STAR_PARAMS: Extra parameters to pass to STAR
    default: 

CUFFLINKS_PARAMS: Cufflinks extra parameters. F.i. '--library-type fr-firststrand' if dUTPs stranded library were used for the sequencing
    default: 

CUFFQUANT_EXTRA_PARAMS: Cuffquant parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA  --multi-read-correct --max-bundle-frags 9999999
    default: 

CUFFDIFF_EXTRA_PARAMS: Cuffdiff parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA  --multi-read-correct
    default: 

CUFFNORM_EXTRA_PARAMS: Extra parameters to use if using Cuffnorm
    default: --output-format cuffdiff

CIRI_EXTRA_PARAMS: CIRI additional parameters
    default: 

PREPROCESSOR: The preprocessing method
    default: trimmomatic

PREPROCESSOR_PARAMS: Read preprocessor extra parameters. F.i. if Trimmomatic, an empty string defaults to MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30 
    default: 

TOGGLE_TRANSCRIPTOME_RECONSTRUCTION: Set True to enable transcriptome reconstruction. Default only quantifies genes and transcripts from the given annotation GTF file
    default: False

DIFF_EXP: Set True to enable differential expression computation for linear genes/transcripts. Only available if more than one sample and more than one condition are given. N.B: differential expression tests for circRNAs is not yet implemented
    default: False

READSTAT_METHODS: Comma separated list of methods to use for read statistics. Currently supported: fastqc,fastx
    default: fastqc

MIN_METHODS: Number of methods that commmonly detect a circRNA to define the circRNA as reliable. If this value exceeds the number of methods specified, it will be set to the number of methods.
    default: 2

MIN_READS: Number of reads to consider a circRNA as expressed
    default: 2

BYPASS_LINEAR: Skip analysis of linear transcripts. This will also skip the analysis of linear-to-circular expression correlation
    default: False
```

## Run the analysis
To trigger the analyses you simply have to call the `./circompara` script in the analysis directory. Remember that if you use the `vars.py` option file, this has to be in the analysis directory. 

    cd /home/user/circrna_analysis
    /home/user/circompara/circompara

### Additional options from the Scons engine:

* *Dryrun*: to see which commands will be executed without actually execute them, use the `-n` option. NB: many commands will be listed, so you should redirect to a file or pipe to a reader like `less`

        /path/to/circompara/dir/circompara -n | less -SR

* *Basic execution*:

        /path/to/circompara/dir/circompara

* *Multitasks*: the `-j` option specifies how many tasks can be run in parallel. Caveat: the '-j * CPUS' value should not be greater than the number of CPU cores available. 
    
        /path/to/circompara/dir/circompara_CirComPara -j4

* *Ignore errors*: keep executing the tasks even when some of them fails. Caveat: this can break downstream analyses

        /path/to/circompara/dir/circompara -i

* *Combine options*: to set multiple options you must sorround them with quotes
    
        /path/to/circompara/dir/circompara_CirComPara "-j4 -i"


## Output files

Results regarding circRNAs are reported in `circrna_analyze` directory with a summary reported in `circrna_analyze/circRNAs_analysis.html` file. Gene expression tables are saved in `cuffdiff` directory.

<!-- Gene/transcript expression estimation and differential expression testsare reported in `cuffdiff` directory. See [Cuffdiff manual][cuffdiff_output] for file format reference. 

CircRNAs' expression levels and gene annotation overlaps are reported in `circRNA_collect_results` directory

**TODO**: explain files
 
Transcript sequences are reported in FASTA format in `transcript_sequences/transcripts.fa` file.

**TODO**

Statistics on alignments are reported in `read_stats_collect/read_stats_collect.txt` file.

-->

# Advanced features

## Stranded libraries

Some tools in CirComPara require special parameters to handle properly stranded reads. CirComPara allows to specify such parameters
Example: include the following parameters if you used the Illumina TruSeq Stranded Total RNA Library Prep Kit with Ribo-Zero Human/Mouse/Rat

    HISAT2_EXTRA_PARAMS = "--rna-strandness FR "
    CUFFLINKS_PARAMS = "--library-type fr-firststrand "

<!-- Experimental
## Fusion genes and fusion circular RNAs
If you want to analyze fusion genes and enable detection of fusion circular RNAs (f-circRNAs) you have to include a `translocation` column to the metadata file. This field specifies the genomic coordinates of the gene pair involved in the fusion. You do not have to specify fusion breakpoints as the transcript structure is inferred by the transcript reconstruction algorithm (currently [Cufflinks][]).   
Gene coordinates must be defined as follow:

    chr:start-end:strand&chr:start-end:strand

with strand being either + or -; note that the two genes' coordinates are separated by a `&` character.
More gene pairs can be specified in one row, you just have to separate the pairs by a `#` character
Below, a metadata file example where the MLL-AF4 fusion gene is specified, as well as its "complementary" translocation:

    file,sample,condition,translocation
    /home/user/reads_S1_1.fq,S1,WT,11:34488-128832:+&4:21033000-21243056:+#4:21033000-21243056:+&11:34488-128832:+
    /home/user/reads_S1_2.fq,S1,WT,
    /home/user/reads_S2_1.fq,S2,MU,
    /home/user/reads_S2_1.fq,S2,MU,

As you can note from the exampple, you do not have to specify the gene pair in each line. Moreover, different samples can be set with different gene pairs.
Fusion gene analysis will then been performed on all samples for all the fusion gene set. 




## Advanced parameters: the vars.py file

Type

    circompara_CirComPara -h
    
to show an help of all parameters

## Advanced features output

**TODO** Output for f-circRNAs ...

### Fusion genes
When enabled, the fusion gene analysis will generate "synthetic" chromosomes ... **TODO**    

 -->

# Appendix

## Installation troubleshooting

### System packages dependencies

In a freshly installed Ubuntu Server 16.04 LTS (x64) you need to install the dependency packages listed below (you need root system rights, or ask your system administrator):

    sudo apt-get install python2.7 python-pip python-numpy zlib1g-dev unzip pkg-config libncurses5-dev default-jre r-base-core libcurl4-openssl-dev libxml2-dev libssl-dev libcairo2-dev pandoc

(pandoc >= 1.12.3 is required, try install latest packages from http://github.com/jgm/pandoc/releases/)  
You also need to upgrade `pip` version (`pip v8.1.1` has an issue with the `--install-option` parameter; tested working with `pip v9.0.1`): 

    pip install --upgrade pip

### Manual installation of dependency tools

To run CirComPara you need several software to be available in your system. If the automatic installation does not work for some reason, try to install the required tools in your system. 
Here there is the list of the tools used in CirComPara with the version that we used during development. We do not list each tool dependencies and if you need support for a specific tool, please refer to the relative software support.

Software|Website|Version
--------|-------|-------:
Ubuntu Linux|http://www.ubuntu.com|Precise (12.04.5 LTS) Server
R|http://cran.r-project.org/|3.2.5 (2016-04-14)
Python|http://www.python.org/|2.7.3
Scons|http://www.scons.org|2.5.0
Trimmomatic|http://www.usadellab.org/cms/?page=trimmomatic|0.36
FASTQC|http://www.bioinformatics.babraham.ac.uk/projects/fastqc/|0.11.5
HISAT2|http://ccb.jhu.edu/software/hisat2/index.shtml|2.0.4
CIRCexplorer|http://github.com/YangLab/CIRCexplorer|1.1.10
STAR|http://github.com/alexdobin/STAR|2.5.2a
CIRI|http://ciri.sourceforge.io/|2.0.2
BWA|http://bio-bwa.sourceforge.net/|0.7.15-r1140
find_circ|http://github.com/marvin-jens/find_circ|1.2
Bowtie2|http://bowtie-bio.sourceforge.net/bowtie2/index.shtml|2.2.9
testrealign|http://www.bioinf.uni-leipzig.de/Software/segemehl/|0.1
Segemehl|http://www.bioinf.uni-leipzig.de/Software/segemehl/|0.2.0-418
Cufflinks|http://cole-trapnell-lab.github.io/cufflinks/|2.2.1
BEDtools|http://bedtools.readthedocs.io|2.26.0
Samtools|http://www.htslib.org/|1.3.1
ggplot2|http://ggplot2.org/|2.2.0
data.table|https://cran.r-project.org/web/packages/data.table/index.html|1.10.0
knitr|http://yihui.name/knitr/|1.14.0



<!-- ## Details on CirComPara architecture and implementation

CirComPara is part of the [circompara][circompara_link] project, which is a collection of scripts consisting of modules that can be assembled to compose computational pipelines for RNA-seq data analysis.   
Each circompara module consists of a script written within the [Scons][scons_link] building tool (say a 'Sconscript') that can implement its function also by wrapping other software. Each module can be used either standalone or included in a pipeline of commands, as it was done for CirComPara.

The core engine is the Scons build tool, which manage the various steps of the analysis.

**TODO**
-->

# How to cite
If you used CirComPara for your analysis, please add the following citation to your references:  
  
Gaffo, E., Bonizzato, A., Kronnie, G. te & Bortoluzzi, S. CirComPara: A Multi‐Method Comparative Bioinformatics Pipeline to Detect and Study circRNAs from RNA‐seq Data. Non-Coding RNA 3, 8 (2017). [http://www.mdpi.com/2311-553X/3/1/8][circompara_article]

[scons_link]: http://scons.org/
[circompara_git_link]: http://github.com/egaffo/CirComPara "circompara Git repository"
[circompara_pack_link]: http://github.com/egaffo/CirComPara/releases/latest "circompara package"
[test_data_link]: http://github.com/egaffo/CirComPara
[Cufflinks]: https://github.com/cole-trapnell-lab/cufflinks
[circompara_link]: http://github.com/egaffo/CirComPara
[cuffdiff_output]:http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/#cuffdiff-output-files
[circompara_article]: http://www.mdpi.com/2311-553X/3/1/8

