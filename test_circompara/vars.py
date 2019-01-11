META            = "meta.csv"
GENOME_FASTA    = '../annotation/CFLAR_HIPK3.fa'
ANNOTATION      = '../annotation/CFLAR_HIPK3.gtf' 
CPUS            = "2"
PREPROCESSOR    = "trimmomatic"
PREPROCESSOR_PARAMS = "MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:35 AVGQUAL:30"
CIRCRNA_METHODS = "dcc,ciri,findcirc,testrealign,"\
                  "circexplorer2_star,circexplorer2_segemehl,"\
                  "circexplorer2_bwa,circexplorer2_tophat,circrna_finder"

FIX_READ_HEADER = 'True'
#MIN_READS = 2 #default

## aligners' custom parameters
HISAT2_EXTRA_PARAMS = "--rna-strandness RF" # stranded libraries
## parameters from CIRI
BWA_PARAMS = ['-T', '19', '-c', '1']
## parameters from CIRCexplorer2
SEGEMEHL_PARAMS = ['-M','1']
TOPHAT_PARAMS = ['--zpacker','pigz', '--max-multihits', '1']
## parameters used in DCC manual example
STAR_PARAMS = ['--outFilterMultimapNmax', '1', 
               '--outSJfilterOverhangMin', '15', '15', '15', '15',
               '--alignSJoverhangMin', '15',
               '--alignSJDBoverhangMin', '15',
               '--seedSearchStartLmax', '30',
               '--outFilterScoreMin', '1',
               '--outFilterMatchNmin', '1',
               '--outFilterMismatchNmax', '2',
               '--chimSegmentMin', '15',
               '--chimScoreMin', '15',
               '--chimScoreSeparation', '10',
               '--chimJunctionOverhangMin', '15'] 

## pre-computed index and annotation files 
#GENOME_INDEX = "../indexes/hisat2/CFLAR_HIPK3"
#SEGEMEHL_INDEX = "../indexes/segemehl/CFLAR_HIPK3.idx"
#BWA_INDEX = "../indexes/bwa/CFLAR_HIPK3"
#BOWTIE2_INDEX = "../indexes/bowtie2/CFLAR_HIPK3"
#BOWTIE_INDEX = "../indexes/bowtie/CFLAR_HIPK3"
#STAR_INDEX = "../indexes/star/CFLAR_HIPK3"
#GENEPRED = "../annotation/CFLAR_HIPK3.genePred.wgn"

LIN_COUNTER = 'ccp' #'dcc'

DCC_EXTRA_PARAMS = ['-fg', '-M', '-Nr', 1, 1, '-F', '-ss']
TESTREALIGN_PARAMS = ['-q', '1']
CE2_PARAMS = ['--no-fix']
FINDCIRC_EXTRA_PARAMS = ['--best-qual', '0']
