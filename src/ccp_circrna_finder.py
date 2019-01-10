'''
This SConscript detects and annotates circular RNAs from RNA-seq data according to the
Circrna_finder [1] protocol.

References:
1.Westholm, J. O. et al. Genome-wide Analysis of Drosophila Circular RNAs
Reveals Their Structural and Sequence Properties and Age-Dependent Neural
Accumulation. Cell Reports 9, 1966-1980 (2014).
'''

import os 

Import('*')

try:
    env = cfinder_env.Clone()

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

out_dir = 'circrna_finder'

## use sed to write a BED file with 'single nucleotide' intervals
## NB: stop position is not included in BED format (intervals are 
## 0-based and 'half-open'), so the stop position must be incremented by 1.
## Only the first six tab-separated fields are considered
bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                    '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
                    '''\\1\\t\\2\\t$$((\\2+1))\\t\\4\\t\\5\\t\\6\\n'''\
                    '''\\1\\t\\3\\t$$((\\3+1))\\t\\4\\t\\5\\t\\6"/e' '''\
                    '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''

bks_reads_cmd = '''samtools view -uhS ${SOURCES[0]} | '''\
                '''bedtools intersect -s -bed -abam stdin -b ${SOURCES[1]} | '''

cfinder_cmd = 'postProcessStarAlignment.pl $EXTRA_PARAMS '\
              '--starDir ${SOURCES[0].dir}' + os.path.sep +\
              ' --outDir ${TARGETS[0].dir}' + os.path.sep
cfinder_targets = [os.path.join(out_dir, f) for f in ['filteredJunctions.bed',
                                                      's_filteredJunctions.bed',
                                                      's_filteredJunctions_fw.bed']]
cfinder_sources = [File(env['FUSION_FILE'])]
cfinder = env.Command(cfinder_targets, 
                      cfinder_sources, 
                      cfinder_cmd)

circ_bed = cfinder[1]

bed = env.Command([os.path.join(out_dir, "${SAMPLE}.sn.circ.bed")],
                  [circ_bed],
                  bed_cmd)

bks_reads_cmd = bks_reads_cmd + \
                    '''cut -f 4 | sort | uniq -c | '''\
                    '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
                    '''sort -k1,1nr > ${TARGETS[0]} '''
    
bks_reads = env.Command([os.path.join(out_dir, 
                                      env['SAMPLE'] + '.cfinder.bks.reads')], 
                        [env['ALIGNMENTS'], bed], 
                        bks_reads_cmd)

results = {'CIRCRNAS':          cfinder[0],
           'GTAG_CIRCRNAS':     cfinder[1],
           'FW_GTAG_CIRCRNAS':  cfinder[2],
           'CIRC_SN_BED':       bed,
           'BKS_READS':         bks_reads[0]}

Return('results')

Clean('.', out_dir)

