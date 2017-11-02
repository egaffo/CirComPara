META            = "meta.csv"
GENOME_FASTA    = '$CIRCOMPARA/test_circompara/annotation/CFLAR_HIPK3.fa'
ANNOTATION      = '$CIRCOMPARA/test_circompara/annotation/CFLAR_HIPK3.gtf' 
CPUS            = "3"
PREPROCESSOR    = "trimmomatic"
PREPROCESSOR_PARAMS = "MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:35 AVGQUAL:30"
CIRCRNA_METHODS = "ciri,findcirc,testrealign,"\
                  "circexplorer2_star,circexplorer2_segemehl,"\
                  "circexplorer2_bwa,circexplorer2_tophat"
TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = 'False'
DIFF_EXP = 'False'
DESEQ = 'False'
CUFFNORM_EXTRA_PARAMS = "--output-format cuffdiff"

## aligners' custom parameters
#BWA_PARAMS = ['-T', '19', '-c', '2']
#SEGEMEHL_PARAMS = ['-M','1']
#TOPHAT_PARAMS = ['--zpacker','pigz', '--max-multihits', '2']
#STAR_PARAMS = ['--outFilterMultimapNmax', '2' ] 

