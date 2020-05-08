'''
This SConscript performs circRNAs detection on a RNA-seq sample using different
circRNA detection methods.

Software dependencies are inherited from the CIRCOMPARA-SConscripts used:
 * ccp_testrealign.py
 * ccp_ciri.py
 * ccp_findcirc.py
 * ccp_tophat.py

Imports:
 * env
 * sample_cpus
 * sample_genome_fasta
 * sample_annotation
 * sample_raw_reads
 * sample_segemehl_index
 * ciri_bwa_index
 * ciri_bwa_extra_parameters
 * ciri_script
 * ciri_extra_parameters
 * bowtie2_index
 * star_index
 * gene_pred

'''

import os, re, ast
from collections import defaultdict

ccp_testrealign   = 'ccp_testrealign.py'
ccp_segemehl = 'ccp_segemehl.py'

Import('*')

try:
    env = env_sample_circrna_methods.Clone()
   #TODO
    PRE_FILTER  = False

except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('ANNOTATION', 'Gene annotation (Ensembl GTF)', '')
    vars.Add('GENOME_FASTA', 'The FASTA file with the reference genome', 'genome.fa')
    vars.Add('READS', 'RNA-seq reads. Comma separated list if paired-end', 'reads.fa')
    ## parameters for Segemehl/testrealign
    vars.Add('SEGEMEHL_INDEX', '''The .idx index for segemehl''', 'genome.idx')
    ## parameters for CIRI 
    vars.Add('BWA_INDEX', '''The index of the reference genome for BWA''','/path/to/bwa/index')
    vars.Add('BWA_PARAMS','Extra parameters for BWA','')
    vars.Add('CIRI', 'The full path to the CIRI_vx.x.pl perl script', '')
    ## parameters for find_circ
    vars.Add('BOWTIE2_INDEX', '''The index of the reference genome for BOWTIE2''', 
             '/path/to/bowtie2/index')
    vars.Add('BOWTIE_INDEX', '''The index of the reference genome for BOWTIE''', 
             '/path/to/bowtie/index')
   ## parameters for CIRCexplorer
    vars.Add('STAR_INDEX', 'The directory path where to find Star genome index', 
             '/path/to/Star/index/dir')
    vars.Add('GENEPRED', 'The genome annotation in GenePred format', 'genes.genePred')
    vars.Add('CIRI_EXTRA_PARAMS', 'CIRI additional parameters', '')
    vars.Add('DCC_EXTRA_PARAMS', 'DCC additional parameters', '')

    vars.Add('CIRCRNA_METHODS', 'Comma separated list of circRNA detection methods to use. '\
	     'Use all methods available as default', '')
    vars.Add('CIRC_MAPPING', '''By default (SE), linearly unmapped reads are'''\
             '''aligned as single-end reads to search for circRNA backsplices. Set PE '''\
             '''to align as paired-end reads by each circRNA method aligner. You can also '''\
             '''specify each aligner's mode, or just which aligner has to use the PE mode, '''\
             ''' with the syntax for Python dictionaries {'SE':['ALN1','ALN2'],'PE':['ALN3','ALN4','ALNn']} '''\
             '''or simply {'PE':['ALN1','ALN2']} if you want just ALN1 and ALN2 tu align as PE. '''\
             '''Supported aligners are BWA,SEGEMEHL,STAR and TOPHAT. BOWTIE2 is also supported but '''\
             '''it is run only in single-end mode as it serves only Findcirc. ''',
             '''{'SE':['STAR','TOPHAT','BOWTIE2'],'PE':['BWA','SEGEMEHL']}''')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

    env['READS'] = env['READS'].split(',')

    env['CIRCRNA_METHODS'] = [m.lower() for m in env['CIRCRNA_METHODS'].strip().split(',')]

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

default_circ_mapping_val = '''{'SE':['STAR','TOPHAT','BOWTIE2'],'PE':['BWA','SEGEMEHL']}'''
env.SetDefault(CIRC_MAPPING = defaultdict(set, ast.literal_eval(default_circ_mapping_val)))

all_aligners = set(['BOWTIE2','BWA', 'SEGEMEHL','STAR', 'TOPHAT'])

if str(env['CIRC_MAPPING'])== '':
    ## use empty string as alias for default
    env.Replace(CIRC_MAPPING = defaultdict(set, ast.literal_eval(default_circ_mapping_val)))

if not str(env['CIRC_MAPPING']).upper() in ['SE', 'PE']:
    ## mixed mode
    try:
        env.Replace(CIRC_MAPPING = defaultdict(set, ast.literal_eval(str(env['CIRC_MAPPING']))))
    except ValueError as e:
        print e
        print '''Malformed CIRC_MAPPING string. Please, set as a Python dictionary, '''\
                '''e.g. {'SE':['STAR', 'TOPHAT'],'PE':[]}'''
        exit(-1)
else:
    ## all-the-same mode
    env.Replace(CIRC_MAPPING = defaultdict(set, 
                                           {str(env['CIRC_MAPPING']).upper(): all_aligners}))


## By default align in single-end mode. Remove from SE mode the aligners specified as PE
env['CIRC_MAPPING']['PE'] = set(env['CIRC_MAPPING']['PE'])
env['CIRC_MAPPING']['SE'] = all_aligners
env['CIRC_MAPPING']['SE'] -= env['CIRC_MAPPING']['PE']

if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
   'findcirc' in env['CIRCRNA_METHODS']:
    ## if Findcirc has to be run, it requires a single-end alignment mode
    env['CIRC_MAPPING']['SE'].add('BOWTIE2')

results = []

## default values when method is not set
segemap = None
star = None
bwa = None
tophat = None
circexplorer2_bwa = None
circexplorer2_tophat = None
circexplorer2_segemehl = None
circexplorer2_star = None
ciri = None
find_circ = None
testrealign = None
dcc = None
cfinder = None

sample_name = env['SAMPLE']

env.SetDefault(FIX_READ_HEADER = True)

build_dir = 'circRNAs'

env['STRANDED'] = False
strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")
if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
    env['STRANDED'] = True and not env['UNSTRANDED_CIRCS']

## PREPARE CIRCRNA METHODS' INPUT READS
if len(env['CIRC_MAPPING']['PE']) > 0:
    env['READS_PE'] = [File(f) for f in env['READS']]

if len(env['CIRC_MAPPING']['SE']) > 0:
    
    env['READS_SE'] = env['READS']

    if len(env['READS_SE']) > 1:
        ## if stranded reads we should handle the cases before concatenating,
        ## to have single-end reads from the same strand of the transcript.
        ## An issue comes with TopHat2, which distinguish orientation of 
        ## single-end reads: for instance, if HISAT2 is set to RF (as for dUTP
        ## library protocol), then the second-in-pair matehas the same strnad
        ## of the transcript, whilst the first-in-pair mate has the opposite 
        ## direction. However, TopHat2 would have been set to fr-firststrand
        ## for paired end alignments, so we need to: (i) make a strand
        ## consistent set of reads (from the same strand of thje transcript, as
        ## other aligners will assume that), and (ii) tell TopHat2 that the 
        ## (single-end) reads are from the same strand of the transcript,
        ## indeed reversing the strand parameter (into fr-secondstrand)
        if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
            if 'RF' in env['HISAT2_EXTRA_PARAMS'].split():
                ## reverse complement the first mate
                ## to have all reads in the same direction of the transcript
                mate_idx = 0

                ## switch TopHat2 parameter
                th2_lib_type = 'fr-secondstrand'
                ## set DCC strandnedd
                env.AppendUnique(DCC_EXTRA_PARAMS = ['-ss'])

            elif 'FR' in env['HISAT2_EXTRA_PARAMS'].split():
                ## reverse complement second mate
                mate_idx = 1
                ## switch TopHat2 parameter
                th2_lib_type = 'fr-firststrand'

            revcomp_cmd = 'fastq_rev_comp.py -f $SOURCE | gzip -c > $TARGET'
            revcomp = env.Command(os.path.join(build_dir, 
                                               env['SAMPLE'] + 
                                               '_' +
                                               str(mate_idx+1) +
                                               '.revcomp.fq.gz'), 
                                  env['READS_SE'][mate_idx], 
                                  revcomp_cmd) 
            
            ## replace env variable
            env['READS_SE'][mate_idx] = revcomp
            env.Replace(TOPHAT_PARAMS = ['--library-type', th2_lib_type])

    ## collapse paired reads into single file
    ## indeed processing reads as they were single-end
    cat_reads_target = os.path.join(build_dir,  
                                    env['SAMPLE'] + '.unmappedSE.fq.gz')
    
    ## we can concatenate gzipped files, no need to decompress
    cat_cmd = 'cat ${SOURCES} '
    
    if env['FIX_READ_HEADER']:
        ## trim FASTQ header to the first white space, so
        ## that only the read ids are kept in FASTQ headers.
        ## This prevents errors in downstream read processing.
        ## Recommended for SRA read datasets, which usually
        ## report read length and other info in FASTQ header.
        ## Plus, append suffix to read ids in order to discriminate
        ## mate reads after concatenation.
        cat_cmd = '''trim_read_header.py -s '\\' -f ${SOURCES[0]} '''

        if len(env['READS_SE']) > 1:
            cat_cmd = cat_cmd + '''-r ${SOURCES[1]} '''

        cat_cmd = cat_cmd + '''| gzip -c '''

    ## finally, convert the paired-end reads into single-end alike
    unmapped_reads = env.Command(cat_reads_target,
                               env['READS_SE'],
                               cat_cmd + '> $TARGET')

    ## here we update the READS_SE parameter to point to the single-end converted reads
    env.Replace(READS_SE = unmapped_reads)


if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
   any([f in env['CIRCRNA_METHODS'] for f in ['testrealign', 'circexplorer2_segemehl']]):
    ## SEGEMEHL CIRCRNA
    env_segemehl = env.Clone()

    env_segemehl['READS'] = env['READS_PE']
    if 'SEGEMEHL' in env['CIRC_MAPPING']['SE']:
        env_segemehl['READS'] = env['READS_SE']

    segemap = env.SConscript(os.path.join(build_dir, ccp_segemehl),
                             variant_dir = build_dir, 
                             src_dir = SRC_DIR,
                             duplicate = 0, 
                             exports='env_segemehl')

    Depends(segemap['ALIGNMENTS'], env['READS'])

    ## segemehl Sconscript result keys:
    #'ALIGNMENTS'
    #'UNALIGNED_READS'
    #'SINGLE_SPLITS'
    #'MULTI_SPLITS'
    #'TRANS_SPLIT'
    #'MAPPED_READS' 
    results.append([File(f) for f in segemap.iteritems()])


    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or 'testrealign' in \
	   env['CIRCRNA_METHODS']:

        env_testrealign = env.Clone()
        env_testrealign['BED'] = segemap['SINGLE_SPLITS']
        testrealign_dir = 'testrealign'
        testrealign = env.SConscript(os.path.join(build_dir, testrealign_dir, 
                                                  ccp_testrealign),
                                     variant_dir = os.path.join(build_dir, 
                                                                testrealign_dir), 
                                     src_dir = SRC_DIR,
                                     duplicate = 0, 
                                     exports = '''env_testrealign''')
        results.append(testrealign.values())
	
    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
	   'circexplorer2_segemehl' in env['CIRCRNA_METHODS']:
        circexplorer2_env = env.Clone()
        circexplorer2_env['ALIGNER'] = 'segemehl'
        circexplorer2_env['FUSION_FILE'] = segemap['SINGLE_SPLITS']
        circexplorer2_env['ALIGNMENTS'] = segemap['ALIGNMENTS']
        circexplorer2_segemehl = env.SConscript(os.path.join(build_dir,	
                                                             'ccp_circexplorer2.py'),
                                                variant_dir = build_dir, src_dir = SRC_DIR,
                                                duplicate = 0,
                                                exports = '''circexplorer2_env''')
        results.append(circexplorer2_segemehl.values())


if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
   'findcirc' in env['CIRCRNA_METHODS']:
    ## FIND_CIRC
    ccp_findcirc = 'ccp_findcirc.py'

    env_findcirc = env.Clone()

    ## Findcirc needs single-end alignment mode
    env_findcirc['READS'] = env['READS_SE']

    find_circ = env.SConscript(os.path.join(build_dir, ccp_findcirc),
                               variant_dir = build_dir, src_dir = SRC_DIR, 
                               duplicate = 0, 
                               exports = '''env_findcirc''')

    Depends(find_circ.values(), env['READS'])
    
    results.append(find_circ.values())

if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
   any([f in env['CIRCRNA_METHODS'] for f in ['dcc', 
                                              'circexplorer2_star',
                                              'circrna_finder']]):

    ## alignments with STAR
    ## ALIGN WITH STAR, WITH FUSION SEARCH ENABLED
    star_env = env.Clone()
    ## set Star parameters to enable fusion search
    star_env.AppendUnique(STAR_PARAMS = ['--chimSegmentMin', '10', 
                                         '--chimOutType', 'Junctions'])
    
    ## TODO: snippet code for future development
    ## The additional value 'SeparateSAMold' for --chimOutType is required by STARChip
    ## (perhaps only for fusion detection). Unfortunately, in STAR v2.7.3a,
    ## --chimMultimapNmax > 0 only supports 'Junction' as chimOutType value
    #if 'starchip' in env['CIRCRNA_METHODS']:
    #   star_env['STAR_PARAMS'].append('SeparateSAMold')

    star_env['READS'] = env['READS_PE']
    if 'STAR' in env['CIRC_MAPPING']['SE']:
        star_env['READS'] = env['READS_SE']
   
    star = env.SConscript(os.path.join(build_dir, 'ccp_star.py'),
                          variant_dir = build_dir, src_dir = SRC_DIR,
                          duplicate = 0, 
    		      exports = '''star_env''')	

    Depends(star, env['READS'])
    
    results.append(star)
    
    Chimeric_out_junction = star[0][1]
    #Chimeric_out_sam = star[0][6]
    
    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
        'circexplorer2_star' in env['CIRCRNA_METHODS']:

        ## parse STAR alignments with CIRCexplorer2
        circexplorer2_env = env.Clone()
        circexplorer2_env['FUSION_FILE'] = Chimeric_out_junction
        circexplorer2_env['ALIGNER'] = 'star'
        #circexplorer2_env['ALIGNMENTS'] = Chimeric_out_sam #star[0][0]
        circexplorer2_star = env.SConscript(os.path.join(build_dir, 
        						'ccp_circexplorer2.py'),
        			variant_dir = build_dir, src_dir = SRC_DIR,
                                    duplicate = 0,
        			exports = '''circexplorer2_env''')
        
        results.append(circexplorer2_star.values())	
        Depends(circexplorer2_star.values(), star)

    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
        'dcc' in env['CIRCRNA_METHODS']:

        ## parse STAR alignments with DCC
        dcc_env = env.Clone()
        dcc_env['FUSION_FILE'] = Chimeric_out_junction
        dcc_env['ALIGNMENTS'] = Chimeric_out_junction #Chimeric_out_sam #star[0][0]
        #dcc_env.AppendUnique(EXTRA_PARAMS = env['DCC_EXTRA_PARAMS'])
        ## AppendUnique will fail to correctly add multi-value options
        ## such as -Nr 1 1, since it will not repeat identical values
        dcc_env.Append(EXTRA_PARAMS = env['DCC_EXTRA_PARAMS'])

        dcc = env.SConscript(os.path.join(build_dir, 
        						'ccp_dcc.py'),
        			variant_dir = build_dir, src_dir = SRC_DIR,
                                  duplicate = 0,
        			exports = '''dcc_env''')
        
        results.append([File(f) for f in dcc.values()])	
        Depends(dcc.values(), star)

    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
        'circrna_finder' in env['CIRCRNA_METHODS']:

        ## parse STAR alignments with circrna-finder
        cfinder_env = env.Clone()
        ## save external env PATH before modifying it
        old_env_path = env['ENV']['PATH']

        ## circrna-finder uses old version of samtools  ( < v1.0 ) 
        ## we need to modify the external env PATH to run it
        ## unfortunately, modifying the external ENV has side effect
        ## on the system PATH, no matter if it is modified in a cloned environment :(
        ## we will later revert to the original PATH
        ## TODO: Is this trick working also in parallel task execution??
        cfinder_env.PrependENVPath('PATH', os.path.join(env['ENV']['CIRCOMPARA_HOME'],
                                                        'bin', 'samtools_v0'))
        #cfinder_env['FUSION_FILE'] = Chimeric_out_sam
        cfinder_env['FUSION_FILE'] = Chimeric_out_junction
        cfinder_env['ALIGNMENTS'] = Chimeric_out_junction #Chimeric_out_sam
        cfinder_env.Append(EXTRA_PARAMS = env['CFINDER_EXTRA_PARAMS'])

        cfinder = env.SConscript(os.path.join(build_dir, 
        						'ccp_circrna_finder.py'),
        			variant_dir = build_dir, src_dir = SRC_DIR,
                                  duplicate = 0,
        			exports = '''cfinder_env''')
        
        results.append([File(f) for f in cfinder.values()])	
        Depends(cfinder.values(), star)

        ## revert to the old env PATH once finished circrna_finder
        env['ENV']['PATH'] = old_env_path

if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
	any([f in env['CIRCRNA_METHODS'] for f in ['circexplorer2_bwa', 'ciri']]):

    bwa_env = env.Clone()

    bwa_env['READS'] = env['READS_PE']
    if 'STAR' in env['CIRC_MAPPING']['SE']:
        bwa_env['READS'] = env['READS_SE']

    bwa = env.SConscript(os.path.join(build_dir, 'ccp_bwa.py'), 
                          variant_dir = build_dir, src_dir = SRC_DIR, 
                          duplicate = 0, 
                          exports = 'bwa_env')
    
    Depends(bwa, env['READS'])

    results.append(bwa)

    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
       'circexplorer2_bwa' in env['CIRCRNA_METHODS']:

        RNA_seq_bwa_sam = bwa[0]
        ## parse BWA alignments with CIRCexplorer2
        circexplorer2_env = env.Clone()
        circexplorer2_env['FUSION_FILE'] = RNA_seq_bwa_sam
        circexplorer2_env['ALIGNER'] = 'BWA'
        circexplorer2_env['ALIGNMENTS'] = bwa[0]
        circexplorer2_bwa = env.SConscript(os.path.join(build_dir, 
        						'ccp_circexplorer2.py'),
        			variant_dir = build_dir, src_dir = SRC_DIR,
                                duplicate = 0,
        			exports = '''circexplorer2_env''')
        
        results.append(circexplorer2_bwa.values())
        Depends(circexplorer2_bwa.values(), bwa)

    if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
       'ciri' in env['CIRCRNA_METHODS']:

        ## parse BWA alignments with CIRI
        env_ciri = env.Clone()
        env_ciri['BWA_ALIGN'] = bwa[0] 
        ciri = env.SConscript(os.path.join(build_dir, 'ccp_ciri.py'), 
                              variant_dir = build_dir, src_dir = SRC_DIR, 
                              duplicate = 0, 
                              exports = 'env_ciri')
        
        results.append(ciri.values())
        Depends(ciri.values(), bwa)

if env['CIRCRNA_METHODS'] == [''] or env['CIRCRNA_METHODS'] == '' or \
	'circexplorer2_tophat' in env['CIRCRNA_METHODS'] or \
    'circexplorer2_tophat_pe' in env['CIRCRNA_METHODS']:

    ## align reads usign TopHat-Fusion
    env_tophat = env.Clone()

    env_tophat['READS'] = env['READS_PE']
    if 'STAR' in env['CIRC_MAPPING']['SE']:
        env_tophat['READS'] = env['READS_SE']

    env_tophat['TOPHAT_INDEX']  = env['BOWTIE_INDEX']
    env_tophat.AppendUnique(TOPHAT_PARAMS = ['--fusion-search', 
                                             '--keep-fasta-order',
		       	                             '--no-coverage-search',
                                             '--bowtie1'])
    if not env['ANNOTATION'] == '' and \
            not ('--GTF' in env_tophat['TOPHAT_PARAMS'] or '-G' in env_tophat['TOPHAT_PARAMS']):
                env_tophat.AppendUnique(TOPHAT_PARAMS = ['--GTF', env['ANNOTATION']])

    tophat = env.SConscript(os.path.join(build_dir, 'ccp_tophat.py'),
                              variant_dir = build_dir, src_dir = SRC_DIR,
                              duplicate = 0, exports = '''env_tophat''')
    
    Depends(tophat['ALIGNMENTS'], env['READS'])

    results.append(tophat['ALIGNMENTS'])
    results.append(tophat['MAPPED_READS'])

    ## parse TopHat alignments with CIRCexplorer2
    circexplorer2_env = env.Clone()
    circexplorer2_env['FUSION_FILE'] = tophat['ALIGNMENTS']
    circexplorer2_env['ALIGNER'] = 'tophat'
    circexplorer2_env['ALIGNMENTS'] = tophat['ALIGNMENTS']
    if len(env_tophat['READS']) > 1:
        circexplorer2_env.Replace(ALIGNER = 'tophat_pe')
    circexplorer2_tophat = env.SConscript(os.path.join(build_dir, 
    						'ccp_circexplorer2.py'),
    			variant_dir = build_dir, src_dir = SRC_DIR,
                            duplicate = 0,
    			exports = '''circexplorer2_env''')
    
    results.append(circexplorer2_tophat.values())
    Depends(circexplorer2_tophat.values(), tophat['ALIGNMENTS'])

Clean('.', build_dir)

circ_met_dict = {'CE2_BWA': circexplorer2_bwa,
                 'CE2_TH' : circexplorer2_tophat,
                 'CE2_SE' : circexplorer2_segemehl,
                 'CE2_ST' : circexplorer2_star,
                 'CIRI'   : ciri,
                 'FC'     : find_circ,
                 'TR'     : testrealign,
                 'DCC'    : dcc,
                 'CFINDER': cfinder
                }

circ_aln_dict = {'SEGEMEHL_MAP' : segemap,
                 'STAR_MAP'     : star,
                 'BWA_MAP'      : bwa, 
                 'TOPHAT_MAP'   : tophat
                }

results_dict = {'CIRC_METHODS': circ_met_dict,
                'CIRC_ALIGNERS':circ_aln_dict
               }
Return('results results_dict')
