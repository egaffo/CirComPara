'''
This SConscript performs Cufflinks [1] on a single sample.

[1] Trapnell, C. et al. 
    
    Differential gene and transcript expression analysis of RNA-seq 
    experiments with TopHat and Cufflinks. 
    
    Nature Protocols 7, 562-578 (2012).

Software requirements inherited from the called SConscripts:
 * Cufflinks v2.2.1

Variables to export when calling from a SConscript:
 * env
 * mapping_file
 * sample_name
 * cufflinks_annotation
 * cufflinks_cpus
 * cufflinks_params

'''
import os

Import('*')

try:
    env = env.Clone()
    mapfile     = mapping_file
    mapfile_basename = sample_name
    annotation  = cufflinks_annotation
    CPUS        = cufflinks_cpus
    CUFFLINKS_PARAMS = cufflinks_params
    TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = cufflinks_toggle_transcriptome_reconstruction
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('ANNOTATION', 'The GFF/GTF file with gene annotation (e.g. from Ensembl)', 'exons.gtf') 
    vars.Add('ALIGNMENTS', 'The read alignment file in SAM/BAM format', 'sample.bam')
    vars.Add('SAMPLE', 'The sample name', 'sample')
    vars.Add('CUFFLINKS_PARAMS', '''Cufflinks extra parameters. '''\
             '''F.i. '--library-type fr-firststrand' if dUTPs stranded library were used '''\
             '''for the sequencing''', '')
    vars.Add('TOGGLE_TRANSCRIPTOME_RECONSTRUCTION', 'Set True to enable transcriptome '\
	     'reconstruction. Default only quantifies genes and transcripts from the given '\
	     'annotation GTF file', 'False')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    mapfile_basename = env['SAMPLE']
    mapfile = env['ALIGNMENTS']
    annotation = env['ANNOTATION']
    CPUS = env['CPUS']
    CUFFLINKS_PARAMS = env['CUFFLINKS_PARAMS']
    TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = False
    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION'] == 'True':
    	TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = True

if os.path.splitext(annotation)[1] == '.gz':
    annotation = env.Command(os.path.basename(os.path.splitext(annotation)[0]), 
                             annotation, 'zcat $SOURCE > $TARGET')

infer_transcripts_param = '-G'
if TOGGLE_TRANSCRIPTOME_RECONSTRUCTION: infer_transcripts_param = '-g'

cufflinks_cmd = "cufflinks -q --no-update-check " + CUFFLINKS_PARAMS + " $(-p " + CPUS +\
                "$) --seed 0 " + infer_transcripts_param + " ${SOURCES[1]} ${SOURCES[0]} -o " +\
		Dir('.').path

cufflinks_targets = ['transcripts.gtf', 
                     'isoforms.fpkm_tracking', 
                     'genes.fpkm_tracking', 
                     'skipped.gtf']

cufflinks = env.Command(cufflinks_targets, [mapfile, annotation], cufflinks_cmd)

Return('cufflinks')
