Import('*')


## NEEDS: 'gffread' from Cufflinks packages
#transcript_sequences_annotation: the GFF annotation (e.g. merge.gtf from cuffmerge)
#transcript_sequences_genome: the genome in FASTA format (single multifasta preferred,not compressed)
try:
    env = env_transcript_fastas.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('GENOME_FASTA', 'Genome sequence in multi-FASTA format', '')
    vars.Add('ANNOTATION', 'Gene/transcript annotation in GFF/GTF/BED format', '')
    ## TODO
    #vars.Add('UNSTRANDED', '', '')
    
    env = Environment(ENV = os.environ, variables = vars)
    
    Help(vars.GenerateHelpText(cmdline_env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

env.SetDefault(UNSTRANDED = False)

fastas_cmd = 'gffread -w ${TARGET} -g ${SOURCES[0]} ${SOURCES[1]}'
fastas = env.Command('transcripts.fa', 
                     [env['GENOME_FASTA'], env['ANNOTATION']], 
                     fastas_cmd)

if env['UNSTRANDED']:
    ## TODO: remove -s parameter from bedtools nuc command
    pass

## calculate GC content and other statistics with Bedtools nuc
nuc_cmd = 'bedtools nuc -fi ${SOURCES[0]} -bed ${SOURCES[1]} -s | gzip -c > ${TARGETS[0]}'
nuc = env.Command('nuc_content.tsv.gz', 
                  [env['GENOME_FASTA'], env['ANNOTATION']],
                  nuc_cmd)

Return('fastas nuc')

