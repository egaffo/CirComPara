'''
This SConscript detects and annotates circular RNAs from RNA-seq data according to the
DCC [1] protocol.

References:
1. Cheng, J., Metge, F. & Dieterich, C. Specific identification and
quantification of circular RNAs from sequencing data. Bioinformatics 32,
1094-1096 (2016).
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
    vars.Add('ALIGNMENTS', 'The Star alignment file. It is used to '\
                           'extract the backsplice read names', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

out_dir = 'dcc'
#env['TMP_DIR'] = os.path.join(out_dir, "_tmp_DCC")
env['TMP_DIR'] = '_tmp_DCC'

bks_reads_cmd = '''zcat ${SOURCES[0]} | samtools view -Su - | '''\
                '''bedtools intersect -s -bed -abam stdin -b ${SOURCES[1]} | '''

## use sed to write a BED file with 'single nucleotide' intervals
## NB: stop position is not included in BED format (intervals are 
## 0-based and 'half-open'), so the stop position must be incremented by 1.
## Only the first six tab-separated fields are considered
bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                    '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
                    '''\\1\\t\\2\\t$$((\\2+1))\\t\\4\\t\\5\\t\\6\\n'''\
                    '''\\1\\t\\3\\t$$((\\3+1))\\t\\4\\t\\5\\t\\6"/e' '''\
                    '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''

bks_reads_cmd = '''bedtools intersect -s -bed -abam ${SOURCES[0]} -b ${SOURCES[1]} | '''

## TODO: implement paired-end reads DCC
dcc_cmd = 'DCC $EXTRA_PARAMS $( -T $CPUS $) -D -O ${TARGETS[0].dir} '\
          '$( -t ${TARGETS[0].dir}' + os.path.sep +'$TMP_DIR $) $SOURCE'
dcc_targets = [os.path.join(out_dir, f) for f in ['CircRNACount',
                                                  'CircCoordinates']]
dcc_sources = [File(env['FUSION_FILE'])]
dcc = env.Command(dcc_targets, 
                  dcc_sources, 
                  dcc_cmd)

## as default, use the CircCoordinates file 
## to collect backsplice reads
circ_bed = dcc[1] ## CircCoordinates file is in BED format

#results.append(dcc)

bed = env.Command([os.path.join(out_dir, "${SAMPLE}.sn.circ.bed")],
                  [circ_bed],
                  bed_cmd)

#results.append(bed)

bks_reads_cmd = bks_reads_cmd + \
                    '''cut -f 4 | sort | uniq -c | '''\
                    '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
                    '''sort -k1,1nr > ${TARGETS[0]} '''
    
bks_reads = env.Command([os.path.join(out_dir, 
                                      env['SAMPLE'] + '.dcc.bks.reads')], 
                        [env['ALIGNMENTS'], bed], 
                        bks_reads_cmd)

#results.append(bks_reads)

results = {'CIRCRNA_COUNT': dcc[0],
           'CIRC_COORDINATES': dcc[1],
           'CIRC_SN_BED': bed,
           'BKS_READS': bks_reads}

Return('results')

Clean('.', out_dir)

