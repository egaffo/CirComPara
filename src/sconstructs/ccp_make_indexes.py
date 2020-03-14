'''
This Sconscript generates genome index file for each read aligner specified in
input. Also, it alllows conversion of gene annotation from GTF to GenePred
format.
'''

import os 

## GET PROGRAM ARGUMENTS
vars = Variables('vars.py')
vars.Add('INDEXES', 'A comma separated list of read aligner programs for which'\
                    ' the genome index has to be built. Also set GENEPRED if'\
                    'you want to convert annotation from GTF to genePred format', 
         'HISAT2,BWA,BOWTIE,BOWTIE2,SEGEMEHL,STAR,GENEPRED')
vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
                   'Comma separated.', 'genome.fa')
vars.Add('HISAT2_EXTRA_PARAMS', 'Extra parameters for htseq2-build', '')
vars.Add('BWA_EXTRA_PARAMS', 'Extra parameters for bwa index', '')
vars.Add('BOWTIE2_EXTRA_PARAMS', 'Extra parameters for bowtie2-build', '')
vars.Add('BOWTIE_EXTRA_PARAMS', 'Extra parameters for bowtie-build', '')
vars.Add('SEGEMEHL_EXTRA_PARAMS', 'Extra parameters for segemehl.x -x', '')
vars.Add('STAR_EXTRA_PARAMS', 'Extra parameters for STAR --genomeGenerate', '')
vars.Add('ANNOTATION', 'GTF annotation file. Needed to create genePred file.',\
                       '')

env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                  variables=vars)
Help(vars.GenerateHelpText(env))
unknown = vars.UnknownVariables()
if unknown:
    print "Run sample: unknown variables", unknown.keys()
    Exit(1)

env['SCONSCRIPT_HOME'] = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

## BUILD GENOME INDEXES FOR READ ALIGNERS
env['GENOME'] = File(env['GENOME']).abspath
env_build_indexes = env.Clone()
indexes = SConscript('ccp_build_indexes.py',
                            src_dir = env['SCONSCRIPT_HOME'],
                            variant_dir = Dir(".").path, duplicate = 0,
                            exports = '''env_build_indexes ''')

## COMPOSE STRINGS TO APPEND TO THE VARS.PY FILE
lines = [Value('GENOME_FASTA="' + env['GENOME'] + '"')]

if "HISAT2" in env['INDEXES'].split(","):
    HISAT2_INDEX	= os.path.abspath(str(indexes['HISAT2'][0]))[0:-6] #remove the .1.ht2 suffix
    lines.append(Value('GENOME_INDEX="' + HISAT2_INDEX + '"'))

if "SEGEMEHL" in env['INDEXES'].split(","):
    SEGEMEHL_INDEX	= os.path.abspath(str(indexes['SEGEMEHL'][0]))
    lines.append(Value('SEGEMEHL_INDEX="' + SEGEMEHL_INDEX + '"'))

if "BWA" in env['INDEXES'].split(","):
    BWA_INDEX	= os.path.abspath(str(indexes['BWA'][0]))[0:-4] #remove the .amb suffix
    lines.append(Value('BWA_INDEX="' + BWA_INDEX + '"'))

if "BOWTIE2" in env['INDEXES'].split(","):
    BOWTIE2_INDEX	= os.path.abspath(str(indexes['BOWTIE2'][0]))[0:-6] #remove .1.bt2 suffix
    lines.append(Value('BOWTIE2_INDEX="' + BOWTIE2_INDEX + '"'))

if "BOWTIE" in env['INDEXES'].split(","):
    BOWTIE_INDEX	= os.path.abspath(str(indexes['BOWTIE'][0]))[0:-7] #remove .1.bt2 suffix
    lines.append(Value('BOWTIE_INDEX="' + BOWTIE_INDEX + '"'))

if "STAR" in env['INDEXES'].split(","):
    STAR_INDEX	= os.path.dirname(os.path.abspath(str(indexes['STAR'][0]))) #index dir
    lines.append(Value('STAR_INDEX="' + STAR_INDEX + '"'))

## CREATE genePred FILE FOR CIRCexplorer IF ONLY GTF ANNOTATION WERE PROVIDED
annotation_dir = "annotation"
if "GENEPRED" in env['INDEXES'].split(","):
    if not env['ANNOTATION'] == '':
        genePred_targets = [os.path.join(annotation_dir, t) for t in ['genePred.transcripts.info',
                                                               '${SOURCES[0].filebase}.genePred',
                                                               '${SOURCES[0].filebase}.genePred.wgn']]
        genePred = env.Command(genePred_targets, File(env['ANNOTATION']),
                               ['gtfToGenePred -infoOut=${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}',
                                'cut -f9 ${TARGETS[0]} | grep -v geneName | '\
                                'paste - ${TARGETS[1]} | sed "s_^\t\([^\t]*\)\t\(.*\)_\1\t\1\t\2_" '\
                                '> ${TARGETS[2]}']
                              )
        lines.append(Value('GENEPRED="' + os.path.abspath(str(genePred[2])) + '"'))
        lines.append(Value('ANNOTATION="' + File(env['ANNOTATION']).abspath + '"'))
    else:
        if not GetOption('help'):
            print "ERROR: annotation in GTF format is required to generate genePred"\
                  "file. Please, set the ANNOTATION parameter."

## SAVE FULL PATHS TO A FILE
def writeVars(target, source, env):    
    with open(target[0].path, "w") as var:
        var.write("\n".join([str(p) for p in source]))
        var.write("\n")
    return None

write_var = env.Command("annotation_vars.py", lines, writeVars)

Depends(write_var, indexes.values())
if "GENEPRED" in env['INDEXES'].split(","):
    if not env['ANNOTATION'] == '':
        Depends(write_var, genePred)

Clean('.', annotation_dir)
