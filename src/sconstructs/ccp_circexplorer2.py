'''
This SConscript detects and annotates circular RNAs from RNA-seq data according to the
CIRCExplorer2 [1] protocol.

Software dependencies are inherited from the CIRCOMPARA-SConscripts used:
* 

When called from a SConscript it imports the following variables:
* circexplorer2_env

References:
1. Zhang, X.-O.; Dong, R.; Zhang, Y.; Zhang, J.-L.; Luo, Z.; Zhang, J.;
Chen, L.-L.; Yang, L. Diverse alternative back-splicing and alternative
splicing landscape of circular RNAs. Genome Res. 2016, 26, 1277-1287,
doi:10.1101/gr.202895.115.
'''

import os 

Import('*')

try:
    env = circexplorer2_env.Clone()

except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('ALIGNER', 'The aligner used to map backsplices. Supported'\
                        'aligners: STAR,BWA,segemehl,TopHat-Fusion,TopHat-Fusion-PE,'\
                        'Mapsplice. Use TopHat-Fusion-PE if paired-end reads.',
             'STAR')
    vars.Add('FUSION_FILE', 'The file to be parsed, depending on the aligner: '\
                            'Chimeric.out.junction for STAR, RNA_seq_bwa.sam '\
                            'for BWA, splicesites.bed for segemehl/testrealign, '\
                            'tophat_fusion/accepted_hits.bam for TopHat-Fusion, '\
                            'mapsplice_out/fusions_raw.txt for Mapsplice', 
            'Chimeric.out.junction')
    vars.Add('GENEPRED', 'The genome annotation in GenePred format', 'genes.genePred')
    vars.Add('GENOME_FASTA', 'The FASTA file with the reference genome', 'genome.fa')
    vars.Add('ALIGNMENTS', 'The alignment file', '')
    vars.Add('CE2_PARAMS', 'Parameters to pass to CIRCexplorer2 annotate', '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

out_dir = 'CIRCexplorer2_' + env['ALIGNER'].lower()
results = []

env['OPTIONS'] = ''

## use sed to write a BED file with 'single nucleotide' (sn) intervals
## NB: stop position is not included in BED format (intervals are 
## 0-based and 'half-open'), so the stop position must be incremented by 1
## to generate a 1-base interval.
## Moreover, the circrna's end position must be decreased by 1 when it is 
## reported as the sn start.
## Only the first six tab-separated fields are considered
## Finally, in name field set circRNA id defined as chr:start-end:strand
bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                    '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
                    '''\\1\\t\\2\\t$$((\\2+1))\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6\\n'''\
                    '''\\1\\t$$((\\3-1))\\t\\3\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6"/e' '''\
                    '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''

bedsamwawb_parse = '''cut -f4,10 | '''\
                   '''sed -r 's/([^:]+):([0-9]+)-([0-9]+):([^\\t])\\t(.*)'''\
                            '''/\\1\\t\\2\\t\\3\\t\\5\\t0\\t\\4/' '''

#TopHat-Fusion, STAR, MapSplice, BWA, segemel
if env['ALIGNER'].lower() == 'star':
    env.Replace(ALIGNER = 'STAR')
    bks_sfx = ".ce2star"
    chimout_bed_cmd = '''chimoutjunc_to_bed.py -i ${SOURCES[0]} | '''\
                     '''gzip -c > $TARGET'''
    chimout_bed = env.Command([os.path.join(out_dir, 'chimoutjunc.bed.gz')], 
                              [File(env['FUSION_FILE'])], 
                              chimout_bed_cmd)

    env.Replace(ALIGNMENTS = chimout_bed)

    circ_reads_cmd = '''get_ce2_star_bks_reads.R -r ${SOURCES[0]} -c ${SOURCES[1]} -g 10 -o $TARGET'''

if env['ALIGNER'].lower() == 'bwa':
    env.Replace(ALIGNER = 'BWA')
    bks_sfx = ".ce2bwa"

    ## alignment file is gzip'd SAM
    bwa_reads_cmd = '''zcat ${SOURCES[0]} | get_ce2_bwa_circ_reads.py -i - > $TARGET'''
    bwa_reads = env.Command(os.path.join(out_dir, 'unfiltered_ce2_bwa_bks_reads.bed.gz'), 
                            env['FUSION_FILE'], 
                            bwa_reads_cmd)

    env.Replace(ALIGNMENTS = bwa_reads)
    
    circ_reads_cmd = '''get_ce2_bwa_bks_reads.R -r ${SOURCES[0]} -c ${SOURCES[1]} -o $TARGET'''

if env['ALIGNER'].lower() == 'segemehl':
    env.Replace(ALIGNER = 'segemehl')
    bks_sfx = ".ce2seg"

    ## will use segemehl sample.sngl.bed file to get circular read ids
    env.Replace(ALIGNMENTS = env['FUSION_FILE'])

    ## for segemehl >= v0.3.0 modify the input BED file
    fixed_bed_cmd = '''grep ';B\\|C;' ${SOURCES} | cut -f1,2,3,6 | sort | '''\
                    '''uniq -c | sed -r 's/ *([0-9]+) ([^\\t]+)\\t([^\\t]+)'''\
                    '''\\t([^\\t]+)\\t([^\\t]+).*/echo "\\2\\t$$((\\3+1))\\t\\4\\tsplits:'''\
                    '''\\1:\\1:\\1:C:P\\t0\\t\\5"/e' > $TARGET'''

    fixed_bed = env.Command(os.path.join(out_dir + '_tmp', 'fixed_splicesites.bed'),
                            env['FUSION_FILE'],
                            fixed_bed_cmd)
    env.Replace(FUSION_FILE = fixed_bed)

    ## segemehl coordinates seem to refer to the intron and not the
    ## spliced exons. So, CE2 start position must be increased to 
    ## comply with segemehl SAM.
    ## Moreover, strand seems to be always +
    bed_cmd = '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                        '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+).*/echo "'''\
                        '''\\1\\t$$((\\2+1))\\t$$((\\2+2))\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6\\n'''\
                        '''\\1\\t$$((\\3-1))\\t\\3\\t\\1:\\2-\\3:\\6\\t\\5\\t\\6"/e' '''\
                        '''${SOURCES[0]} | sort -k1,1 -k2,2n -k3,3n > $TARGET'''

    ###TODO: should -s (strandness) be considered in bedtools intersect
    ### commmand, since CIRCexplorer annotate may change strand field?
    #circ_reads_cmd = '''grep ";B\\|C;" ${SOURCES[0]} | '''\
    #                 '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
    #                           '''([^;]+);([^;]+);([^;]+);([^\\t]+)\\t'''\
    #                           '''([^\\t]+)\\t([^\\t]+)/echo "'''\
    #                           '''\\1\\t$$((\\2+1))\\t\\3\\t\\6\\t\\8\\t\\9"/e' | '''\
    #                 '''sort | uniq | '''\
    #                 '''bedtools intersect -s -a ${SOURCES[1]} -b stdin -wa -wb | '''\
    #                  + bedsamwawb_parse + ''' | gzip -c > $TARGET'''

    circ_reads_cmd = '''get_ce2_segemehl_bks_reads.R -r ${SOURCES[0]} -c ${SOURCES[1]} -o $TARGET'''

if env['ALIGNER'].lower() == 'tophat':
    env.Replace(ALIGNER = 'TopHat-Fusion')
    bks_sfx = ".ce2th"

if env['ALIGNER'].lower() == 'tophat_pe':
    env.Replace(ALIGNER = 'TopHat-Fusion')
    env.Replace(OPTIONS = '--pe')
    bks_sfx = ".ce2th"

if env['ALIGNER'] == 'TopHat-Fusion':

    ### first, filter TopHat-Fusion alignments to get only fusion alignments
    ### then intersect with backsplices.
    ### 260=256+4 i.e: secondary alignment and not mapped.
    ### Keep SAM header in grep (required by bedtools).
    ### Mind that we want BAM input to bedtools being the B file,
    ### check bedtools intersect note (1) "When a BAM file is used for the A
    ###  file, the alignment is retained if overlaps exist, and exlcuded if an
    ###  overlap cannot be found.  If multiple overlaps exist, they are notreported,
    ###  as we are only testing for one or more overlaps."
    #circ_reads_cmd = '''samtools view -hF 260 ${SOURCES[0]} | grep "^@\\|XF:Z\\|XP:Z" | '''\
    #                 '''samtools view -uS - | bedtools intersect -s -bed '''\
    #                 '''-a ${SOURCES[1]} -b stdin -wa -wb | '''\
    #                 + bedsamwawb_parse + ''' | gzip -c > $TARGET'''
    circ_reads_cmd = '''get_ce2_th_bks_reads.R -r ${SOURCES[0]} -c ${SOURCES[1]} -o $TARGET'''
 
if env['ALIGNER'].lower() == 'mapsplice':
    env.Replace(ALIGNER = 'MapSplice')
    bks_sfx = ".ce2ms"
    ##TODO: circ_reads_cmd
 
CIRCexplorer2_cmd = 'CIRCexplorer2 parse $OPTIONS -b $TARGET -t $ALIGNER $SOURCE'
CIRCexplorer2_targets = os.path.join(out_dir, 'back_spliced_junction.bed')
CIRCexplorer2_sources = [File(env['FUSION_FILE'])]
CIRCexplorer2 = env.Command(CIRCexplorer2_targets, 
                            CIRCexplorer2_sources, 
                            CIRCexplorer2_cmd)

## as default, use the fusion_junction.bed file 
## to collect backsplice reads
circ_bed = CIRCexplorer2[0]

if env['ALIGNER'] == 'TopHat-Fusion':
    env.Replace(ALIGNMENTS = circ_bed)

results.append(CIRCexplorer2)

if not env['GENEPRED'] == '':
    CIRCexplorer2_annotate_targets = os.path.join(out_dir, 'annotate',
                                                  'circularRNA_known.txt')
    CIRCexplorer2_annotate_sources = [File(env['GENEPRED']),
                                      File(env['GENOME_FASTA']),
                                      CIRCexplorer2[0]]
    CIRCexplorer2_annotate_cmd = ' && '.join(['cd ${TARGETS[0].dir}', 
                                 'CIRCexplorer2 annotate $CE2_PARAMS -r '\
                                 '${SOURCES[0].abspath} -g '\
                                 '${SOURCES[1].abspath} -b '\
                                 '${SOURCES[2].abspath} -o $TARGET.file',
                                 'cd ' + Dir('#').abspath])
    CIRCexplorer2_annotate = env.Command(CIRCexplorer2_annotate_targets,
                                         CIRCexplorer2_annotate_sources,
                                         CIRCexplorer2_annotate_cmd)

    ## use the filtered annotate/circularRNA_known.txt BED file
    ## to collect backsplice reads
    circ_bed = CIRCexplorer2_annotate

    results.append(CIRCexplorer2_annotate)

if env['ALIGNER'].lower() in ['star', 'segemehl', 'bwa', 'tophat-fusion']:
    bed = circ_bed
else:
    bed = env.Command([os.path.join(out_dir, "${SAMPLE}.sn.circ.bed")],
                      [circ_bed],
                      bed_cmd)

results.append(bed)

## get backsplice read ids
circ_reads = env.Command([os.path.join(out_dir,
                                      '${SAMPLE}.circular.reads.bed.gz')], 
                         [env['ALIGNMENTS'], bed], 
                         circ_reads_cmd)

env.Replace(ALIGNMENTS = circ_reads)
    
## collect all read ids
bks_reads_cmd = '''zcat ${SOURCES[0]} | cut -f 4 | sort | uniq -c | '''\
                '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
                '''sort -k1,1nr > ${TARGETS[0]} '''
    
bks_reads = env.Command([os.path.join(out_dir, 
                                      env['SAMPLE'] + bks_sfx + '.bks.reads')], 
                        [env['ALIGNMENTS'], bed], 
                        bks_reads_cmd)

results.append(bks_reads)

results = {'CIRCRNAS':          circ_bed[0],
           'CIRC_SN_BED':       bed,
           'BKS_READS':         bks_reads[0],
           'CIRC_READS':        circ_reads[0]}

Return('results')

Clean('.', out_dir)

