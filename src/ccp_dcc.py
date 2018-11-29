'''
This SConscript detects and annotates circular RNAs from RNA-seq data according to the
DCC [1] protocol.

References:
1. Cheng, J., Metge, F. & Dieterich, C. Specific identification and
quantification of circular RNAs from sequencing data. Bioinformatics 32,
1094–1096 (2016).
'''

import os 

Import('*')

try:
    env = dcc_env.Clone()

except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars.Add('FUSION_FILE', 'The file to be parsed: '\
                            'Chimeric.out.junction for STAR', 
            'Chimeric.out.junction')
    vars.Add('EXTRA_PARAMS', 'Extra parameters to be passed to DCC', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

out_dir = 'dcc'
results = []

#bks_reads_cmd = '''zcat ${SOURCES[0]} | samtools view -Su - | '''\
#                '''bedtools intersect -s -bed -abam stdin -b ${SOURCES[1]} | '''
#
### use sed to write a BED file with 'single nucleotide' intervals
### NB: stop position is not included in BED format (intervals are 
### 0-based and 'half-open'), so the stop position must be incremented by 1.
### Only the first six tab-separated fields are considered
#bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
#                    '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
#                    '''\\1\\t\\2\\t$$((\\2+1))\\t\\4\\t\\5\\t\\6\\n'''\
#                    '''\\1\\t\\3\\t$$((\\3+1))\\t\\4\\t\\5\\t\\6"/e' '''\
#                    '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''
#
#bks_reads_cmd = '''bedtools intersect -s -bed -abam ${SOURCES[0]} -b ${SOURCES[1]} | '''

dcc_cmd = 'DCC -D -O $TARGET.dir $SOURCE'
dcc_targets = os.path.join(out_dir, 'CircRNACount')
dcc_sources = [File(env['FUSION_FILE'])]
dcc = env.Command(dcc_targets, 
                  dcc_sources, 
                  dcc_cmd)

## as default, use the fusion_junction.bed file 
## to collect backsplice reads
circ_bed = dcc[0]

results.append(dcc)

#bed = env.Command([os.path.join(out_dir, "${SAMPLE}.sn.circ.bed")],
#                  [circ_bed],
#                  bed_cmd)
#
#results.append(bed)
#
#bks_reads_cmd = bks_reads_cmd + \
#                    '''cut -f 4 | sort | uniq -c | '''\
#                    '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
#                    '''sort -k1,1nr > ${TARGETS[0]} '''
#    
#bks_reads = env.Command([os.path.join(out_dir, 
#                                      env['SAMPLE'] + bks_sfx + '.bks.reads')], 
#                        [env['ALIGNMENTS'], bed], 
#                        bks_reads_cmd)
#
#results.append(bks_reads)

Return('results')

Clean('.', out_dir)

