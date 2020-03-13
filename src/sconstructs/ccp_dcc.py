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

## use sed to write a BED file with 'single nucleotide' intervals
## NB: stop position is not included in BED format (intervals are 
## 0-based and 'half-open'), so the stop position must be incremented by 1.
## Only the first six tab-separated fields are considered
bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                    '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
                    '''\\1\\t\\2\\t$$((\\2+1))\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6\\n'''\
                    '''\\1\\t$$((\\3-1))\\t\\3\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6"/e' '''\
                    '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''

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
circrnas = dcc[0]

## when -F filtering option is used the results loose strand info
## so we have to retrieve it from the coodinates file
if '-F' in env['DCC_EXTRA_PARAMS']:
    circrnas = env.Command(os.path.join(out_dir, 'strandedCircRNACount'),
                           [dcc[0], dcc[1]],
                           'dcc_fix_strand.R -c ${SOURCES[0]} -d ${SOURCES[1]} -o ${TARGET}')
    circ_bed = circrnas


bed = env.Command([os.path.join(out_dir, "${SAMPLE}.sn.circ.bed")],
                  [circ_bed],
                  bed_cmd)

### get backsplice read ids
circ_reads_cmd = '''get_dcc_bks_reads.R -r ${SOURCES[0]} -c ${SOURCES[1]} -o $TARGET'''

## TODO: if "-E" or "--endTol" in env['EXTRA_PARAMs']
#circ_reads_cmd = circ_reads_cmd + ' -t ' + get_endTol_value

if (not '-N' in env['EXTRA_PARAMS']) or '-ss' in env['EXTRA_PARAMS']:
    circ_reads_cmd = circ_reads_cmd + ' -s'

circ_reads = env.Command([os.path.join(out_dir,
                                      '${SAMPLE}.circular.reads.bed.gz')], 
                         [env['ALIGNMENTS'], circ_bed], 
                         circ_reads_cmd)
    
## collect all read ids
bks_reads_cmd = '''zcat ${SOURCES[0]} | cut -f 4 | sort | uniq -c | '''\
                '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
                '''sort -k1,1nr > ${TARGETS[0]} '''
    
bks_reads = env.Command([os.path.join(out_dir, 
                                      env['SAMPLE'] + '.dcc.bks.reads')], 
                        [circ_reads, bed], 
                        bks_reads_cmd)

results = {'CIRCRNAS': circrnas[0],
           'CIRC_COORDINATES': dcc[1],
           'CIRC_SN_BED': bed,
           'BKS_READS': bks_reads[0],
           'CIRC_READS': circ_reads}

Return('results')

Clean('.', out_dir)

