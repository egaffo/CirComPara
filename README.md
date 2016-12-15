---
Title:      Circpipe  
Subtitle:   A pipeline to analyze RNA-seq data, with a focus on circular RNAs  
Project:    Junk2/Circpipe  
Author:     Enrico Gaffo  
Affiliation:    Compgen - University of Padova  
Web:        http://compgen.bio.unipd.it  
Date:       September 22, 2016  
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
---

# Circpipe

Circpipe is a computational pipeline to detect and quantify expression of linear and circular RNAs from RNA-seq data.

<!--TODO: more exhaustive description -->

## Quick install

Execute the following commands to download and install (locally) in your system the scripts and tools required to run Circpipe. 
If something goes wrong with the installation process try to manually install the software as described below.

    git clone ssh://compgen.bio.unipd.it/repo/junk2.git
    cd junk2
    ./junk2_install_circpipe

### Test your installation

NB: in the `sed` commands change the `/full/junk2/dir/path` path with your installation directory 

    wget http://compgen.bio.unipd.it/~enrico/junk2_test_data.tar.gz
    tar xf junk2_test_data.tar.gz
    cd test_data
    mkdir analysis
    sed "s@\$JUNK2_HOME@/full/junk2/dir/path@g" vars.py > analysis/vars.py
    sed "s@\$JUNK2_HOME@/full/junk2/dir/path@g" meta_tran.csv > analysis/meta_tran.csv
    cd analysis
    ../../junk2_circpipe

If you receive some error messages try to follow instructions in **Installation troubleshooting** section.

### Add  Circpipe to your environment

Once completed the installation, if you do not want to type the whole path to the Circpipe executable each time, you can update your `PATH` environment variable. From the terminal type the following command (replace the `/path/to/junk2/install/dir` string with Circpipe's actual path)   

    export PATH=/path/to/junk2/install/dir:$PATH

Another way is to link Circpipe's main script in your local `bin` directory  

    cd /home/user/bin
    ln -s /path/to/junk2/install/dir/junk2_circpipe
    
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

**TODO**: describe META file structure and format in more details



### Specify the reference genome file

A required parameter is the reference genome. You can either pass the reference genome from the command line

    ./junk2_circpipe "GENOME_FASTA='/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'"

or by setting the `GENOME_FASTA` parameter in the `vars.py` file; e.g.:

    GENOME_FASTA = '/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'


### Specify options in vars.py

The `vars.py` must be placed in your analysis directory. 

**TODO**: describe vars.py file

## Run the analysis
To trigger the analyses you simply have to call the `junk2_circpipe` script in the analysis directory. Remember that if you use the `vars.py` option file, this has to be in the analysis directory. 

    cd /home/user/circrna_analysis
    /home/user/junk2/junk2_circpipe

### Common options to run the analysis:

* *Dryrun*: to see which commands will be executed without actually execute them, use the `-n` option. NB: many commands will be listed, so you should redirect to a file or pipe to a reader like `less`

        /path/to/junk2/dir/junk2_circpipe -n | less -SR

* *Basic execution*:

        /path/to/junk2/dir/junk2_circpipe

* *Multiprocess*: the `-j` option
    
        /path/to/junk2/dir/junk2_circpipe -j4

* *Ignore errors*:

        /path/to/junk2/dir/junk2_circpipe -i

* *Combine options*:
    
        /path/to/junk2/dir/junk2_circpipe "-j4 -i"


## Output files

Gene/transcript expression estimation and differential expression testsare reported in `cuffdiff` directory. See [Cuffdiff manual][cuffdiff_output] for file format reference. 

CircRNAs' expression levels and gene annotation overlaps are reported in `circRNA_collect_results` directory

**TODO**: explain files
 
Transcript sequences are reported in FASTA format in `transcript_sequences/transcripts.fa` file.

**TODO**

Statistics on alignments are reported in `read_stats_collect/read_stats_collect.txt` file.

# Advanced features

## Stranded libraries

Some tools in Circpipe require special parameters to handle properly stranded reads. Circpipe allows to specify such parameters
Example: include the following parameters if you used the Illumina TruSeq Stranded Total RNA Library Prep Kit with Ribo-Zero Human/Mouse/Rat

    HISAT2_EXTRA_PARAMS = "--rna-strandness FR "
    CUFFLINKS_PARAMS = "--library-type fr-firststrand "

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


<!--TODO -->

## Advanced parameters: the vars.py file

Type

    junk2_circpipe -h
    
to show an help of all parameters

## Advanced features output

**TODO** Output for f-circRNAs ...

### Fusion genes
When enabled, the fusion gene analysis will generate "synthetic" chromosomes ... **TODO**    

<!--TODO -->

# Appendix

## Installation troubleshooting

### System packages dependencies

In a freshly installed Ubuntu Server 16.04 LTS (x64) you need to install the dependency packages listed below (you need root system rights, or ask your system administrator):

    sudo apt-get install python2.7 python-pip python-numpy zlib1g-dev unzip pkg-config libncurses5-dev default-jre r-base-core
    libcurl4-openssl-dev libxml2-dev libssl-dev libcairo2-dev pandoc
(pandoc >= 1.12.3 is required, try install latest packages from https://github.com/jgm/pandoc/releases/)
You also need to upgrade `pip` version (`pip v8.1.1` has an issue with the `--install-option` parameter; tested working with `pip v9.0.1`) 

    pip install --upgrade pip

### Manual installation of the tools

To run Circpipe you need several software to be available in your system. If the automatic installation does not work for some reason, try to install the required tools in your system. 
Here there is the list of the tools used in Circpipe. We do not list each tool dependencies and if you need support for a specific tool, please refer to the relative software support. 

Tool|Project page
----|------------
Ubunutu Linux|www.ubuntu.com
Python v2.7|www.python.org
R|cran.r-project.org
Scons|www.scons.org
...|...

**TODO**
* list tools with tested version and link


## Details on Circpipe architecture and implementation

Circpipe is part of the [JUNK2][junk2_link] project, which is a collection of scripts consisting of modules that can be assembled to compose computational pipelines for RNA-seq data analysis.   
Each Junk2 module consists of a script written within the [Scons][scons_link] building tool (say a 'Sconscript') that can implement its function also by wrapping other software. Each module can be used either standalone or included in a pipeline of commands, as it was done for Circpipe.

The core engine is the Scons build tool, which manage the various steps of the analysis.

**TODO**

[scons_link]: http://scons.org/
[junk2_git_link]: ssh://compgen.bio.unipd.it/repo/junk2.git "JUNK2 Git repository"
[junk2_pack_link]: http://compgen.bio.unipd.it/junk2/latest/junk2.tar.gz "JUNK2 package"
[test_data_link]: http://a-link/test_data.tar.gz
[Cufflinks]: https://github.com/cole-trapnell-lab/cufflinks
[junk2_link]: http://github.com/egaffo/junk2
[cuffdiff_output]:http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/#cuffdiff-output-files

