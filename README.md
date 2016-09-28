---
Title:      Circpipe  
Subtitle:   A pipeline to analyze RNA-seq data, with a focus on circular RNAs  
Project:    Junk2/Circpipe  
Author:     Enrico Gaffo  
Affiliation:    Compgen - University of Padova  
Web:        http://compgen.bio.unipd.it  
Date:       September 22, 2016  
---

# Circpipe

A pipeline to analyze RNA-seq data, with a focus on circular RNAs.    

Currently, Circpipe is part of the JUNK2 script collection...

The pipeline relies on the [Scons][scons_link] build tool ...

**TODO**

<!--TODO: more exhaustive description -->


##Table of contents


1. [Quick start](#quick-start)
2. [How to install](#how-to-install)
3. [How to use](#how-to-use)
4. [How to interpret the output](#the-output)

<a name="quick-start"></a>
##Quick start

Download and installation
    
    git clone ssh://compgen.bio.unipd.it/repo/junk2.git
    cd junk2
    ./junk2_install_circpipe

Test installation

    wget http://compgen.bio.unipd.it/~enrico/junk2_test_data.tar.gz
    tar xf junk2_test_data.tar.gz
    cd test_data
    mkdir analysis
    sed "s@\$JUNK2_HOME@/full/junk2/dir/path@g" vars.py > analysis/vars.py
    sed "s@\$JUNK2_HOME@/full/junk2/dir/path@g" meta_tran.csv > analysis/meta_tran.csv
    cd analysis
    ../../junk2_circpipe

To run an analysis usually you want to specify your data \(the sequenced reads in FASTQ format\) and a reference genome in FASTA format.

1. *Specify read files in a metadata file.* To set up your analysis compose a comma separated metadata file with the following header and fields:  

        file,sample,condition  

    If you have paired-end sequenced samples write one line per file with the same sample name and condition. E.g: ``meta.csv`` \(as default\)
    
        file,sample,condition
        /home/user/reads_S1_1.fq,S1,WT
        /home/user/reads_S1_2.fq,S1,WT
        /home/user/reads_S2_1.fq,S2,MU
        /home/user/reads_S2_1.fq,S2,MU

2. *Specify the reference genome file.* You can either pass the reference genome from the command line:
        
        ./junk2_circpipe "-i -j2 GENOME_FASTA='/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'"

    or by setting the `GENOME_FASTA` parameter in the `vars.py` file; e.g.:

        GENOME_FASTA = '/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'


<a name="how-to-install"></a>
##How to install
Just clone the [JUNK2 git repository][junk2_git_link] and install the required tools.  
To clone via SSH you can run \(you must have an account to access the server\):    

    git clone ssh://compgen.bio.unipd.it/repo/junk2.git

Or you can download the latest version package from [here][junk2_pack_link] and then extract it

    wget [junk2_pack_link]
    tar xf junk2.tar.gz

<!-- TODO: 
* make repo and package available without authentication
-->

###Automatic tool installation

You can run the auto installation script to download and install all the required third party packages in the local directory. 
    
    cd junk2
    ./junk2_install_circpipe  

NB: if some packages fail to install try to install them manually and then link the executable to the `junk2/bin` directory or let them be available in your system.


Once completed the installation, if you do not want to type the whole path to the Circpipe executable each time, you can update your `PATH` environment variable  

    export PATH=/path/to/junk2/install/dir:$PATH

or you can link it to your local bin directory  

    cd /home/user/bin
    ln -s /path/to/junk2/install/dir/junk2_circpipe

###Software dependencies

To run Circpipe you need several software to be available in your system.    
Here there is the list of the tools used in Circpipe. We do not list each tool dependencies and do not provide support for them...sorry, we have no time for this. 

The core engine is the Scons build tool, which manage the various steps of the analysis.

* Scons   
* R
* ...

**TODO**
<!--TODO 
* list tools with tested version and link
-->


###Test Circpipe
After installation, to test Circpipe download the test data package from [here][test_data_link] and follow the instructions in the README file.

    cd junk2
    wget [test_data_link]
    tar xf test_data.tar.gz
    cd test_data
    ../../junk2_circpipe

This will produce ....

**TODO**

<!--TODO 
* describe results check
-->

<a name="how-to-use"></a>
##How to use
Set project directory:

###1. compose META file

**TODO**: describe META file structure

###2. specify options in vars.py

**TODO**: describe vars.py file

###3. Run the analysis
To trigger the analyses you simply have to call the `junk2_circpipe` script in the analysis directory. Remember that if you use the `vars.py` option file, this has to be in the analysis directory. 

Below we list some common options to run the analysis:

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


###Fusion genes and fusion circular RNAs
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

###Advanced parameters: the vars.py file

###Advanced metadata


<a name="the-output"></a>
##How to interpret the output
Results tables ...     
CircRNAs in `circRNA_collect_results` ....     
**TODO**

###Fusion genes
When enabled, the fusion gene analysis will generate "synthetic" chromosomes ... **TODO**    

<!--TODO -->

[scons_link]: http://scons.org/
[junk2_git_link]: ssh://compgen.bio.unipd.it/repo/junk2.git "JUNK2 Git repository"
[junk2_pack_link]: http://compgen.bio.unipd.it/junk2/latest/junk2.tar.gz "JUNK2 package"
[test_data_link]: http://a-link/test_data.tar.gz
[Cufflinks]: https://github.com/cole-trapnell-lab/cufflinks
