'''
This SConscript performs RNA-seq analysis for each sample specified in the metadata file.
In addition, it merges transcript annotation derived from each samples and assess differential
gene and transcript expression among the samples.
The metadata file given as input must be comma separated and must have the following columns:

file        :the path to the reads file. If paired-end use one row per read file,
             setting the same sample name
sample      :the sample name/ID
condition   :biological condition, used for differential gene expression.
adapter     :read adapter to trim (optional. leave empty string if not applicable)
translocation:sharp (#) separted list of the coordinates involving fusion genes/translocations,
             defined as coords1&coords2 elements (optional)

meta.csv example:

file,sample,condition,adapter,translocation
/home/user/data/reads/SRR445566.fastq.gz,SRR445566,TUMOR,
/home/user/data/reads/SRR534325_1.fastq.gz,SRR534325,CONTROL,/trimmomatic/adapters/TruSeq3-PE-2.fa
/home/user/data/reads/SRR534325_2.fastq.gz,SRR534325,CONTROL,/trimmomatic/adapters/TruSeq3-PE-2.fa
/home/user/data/reads/SRR534326_1.fastq.gz,SRR534325,t(4;1)t(1;X),,4:90000-120000:+&11:20-30000:-#1:20-400:+&X:1000-2045:+
/home/user/data/reads/SRR534326_2.fastq.gz,SRR534325,CONTROL,t(4;1)t(1;X),4:90000-120000:+&11:20-30000:-#1:20-400:+&X:1000-2045:+

'''

import os, csv, itertools, collections, re, errno, ast
from collections import defaultdict

def SymLink(target, source, env):
    try:
        os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.unlink(os.path.abspath(str(target[0])))
            os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))



def get_matching_nodes(nodelist, rexpression, abspath = False):
    files = []        
    for node in Flatten(nodelist):
        path = node.path
        if abspath:
            path = node.abspath
        if re.match(rexpression, path):
            files.append(node)
    return files

## Define an action to write text lines in a file
def writeLines(target, source, env):
    ''' Generate a text file with sources file absolute paths as lines. 
    :param source: the source nodes to be listed in the target file.
    :param target: the output file.'''
    
    with open(target[0].path, 'w') as txtfile:
        for line in source:
            txtfile.write(line.abspath)
            txtfile.write('\n')

    return None

WriteLinesInTxt = Builder(action = writeLines, suffix = '.txt')

## GET PROGRAM ARGUMENTS
variables_file = 'vars.py'
vars = Variables(variables_file)

# basic parameters
vars.Add('META', 
         'The metadata table file where you specify the project samples, etc.', 
         'meta.csv')
vars.Add('ANNOTATION', 
         'Gene annotation file (like Ensembl GTF/GFF)', 
         '')
vars.Add('GENOME_FASTA', 
         'The FASTA file with the reference genome', 
         '')
vars.Add('CIRCRNA_METHODS', 
	 'Comma separated list of circRNA detection methods to use. '\
	 'Repeated values will be collapsed. Currently supported: ciri, find_circ, '\
	 'circexplorer2_star, circexplorer2_bwa, circexplorer2_segemehl, testrealign ('\
	 'circexplorer2_segemehl is actually a filtering of testrealign predictions). '\
	 'Set an empty string to use all methods available (including deprecated methods). ', 
         '')

# non-basic parameters
vars.Add('CPUS', 
         'Set number of CPUs', 
         '4')
vars.Add('GENEPRED', 
         'The genome annotation in GenePred format', 
         '')

## aligners indexes
vars.Add('GENOME_INDEX', 
         '''The index of the reference genome for HISAT2''', 
         '')
vars.Add('SEGEMEHL_INDEX', 
         '''The .idx index for segemehl''', 
         '')
vars.Add('BWA_INDEX', 
         '''The index of the reference genome for BWA''',
         '')
vars.Add('BOWTIE2_INDEX', 
         '''The index of the reference genome for BOWTIE2''',
         '')
vars.Add('BOWTIE_INDEX', 
         '''The index of the reference genome for BOWTIE '''\
         '''when using CIRCexplorer2_tophat''','')

vars.Add('STAR_INDEX', 
         'The directory path where to find Star genome index', 
         '')

## aligners extra parameters
vars.Add('HISAT2_EXTRA_PARAMS', 
         '''Extra parameters to add to the HISAT2 aligner fixed '''\
         '''parameters '--dta --dta-cufflinks --rg-id <SAMPLE> --no-discordant '''\
         '''--no-mixed --no-overlap'. For instance, '--rna-strandness FR' if stranded reads'''\
         ''' are used.''', 
         '')
vars.Add('BWA_PARAMS',
         'Extra parameters for BWA',
         '')
vars.Add('SEGEMEHL_PARAMS', 
         'SEGEMEHL extra parameters', 
         '')
vars.Add('TOPHAT_PARAMS', 
         'Extra parameters to pass to TopHat', 
         '')
vars.Add('STAR_PARAMS', 
         'Extra parameters to pass to STAR', 
         '')
vars.Add('BOWTIE2_PARAMS', 
         'Extra parameters to pass to Bowtie2 in addition to '\
         '-p $CPUS --reorder --score-min=C,-15,0 -q', 
         '')

## linear transcriptome extra parameters
vars.Add('CUFFLINKS_PARAMS', 
         '''Cufflinks extra parameters. '''\
         '''F.i. '--library-type fr-firststrand' if dUTPs stranded library were used '''\
         '''for the sequencing''', 
         '')
vars.Add('CUFFQUANT_EXTRA_PARAMS',
         'Cuffquant parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA '\
         ' --multi-read-correct --max-bundle-frags 9999999', 
         '')
vars.Add('CUFFDIFF_EXTRA_PARAMS',
         'Cuffdiff parameter options to specify. E.g. --frag-bias-correct $GENOME_FASTA '\
         ' --multi-read-correct', 
         '')
vars.Add('CUFFNORM_EXTRA_PARAMS', 
        'Extra parameters to use if using Cuffnorm',
	'--output-format cuffdiff')
vars.Add('STRINGTIE_PARAMS', 
         '''Stringtie extra parameters. '''\
         '''F.i. '--rf' assumes a stranded library fr-firststrand, to be '''\
         '''used if dUTPs stranded library were sequenced''', 
         '')

## circRNA methods' extra parameters
vars.Add('CIRI_EXTRA_PARAMS', 
         'CIRI additional parameters', 
         '')
## DCC extra parameters
vars.Add('DCC_EXTRA_PARAMS', 
         'DCC additional parameters', 
         '')
## CIRCexplorer2 (annotate) parameters
vars.Add('CE2_PARAMS', 
         'Parameters to pass to CIRCexplorer2 annotate', 
         '')
## Segemehl/testrealign extra parameters
vars.Add('TESTREALIGN_PARAMS', 
         'Segemehl/testrealign filtering parameters'\
         '-q indicates the minimum median quality of backsplices ends '\
         '(like the Haarz parameter)', 
         '-q median_1')
## Find_circ parameters
vars.Add('FINDCIRC_EXTRA_PARAMS', 
         'Parameters for find_circ.py. '\
         'Additional parameters: --best-qual INT is used to filter '\
         'find_circ results according to best_qual_left and best_qual_right '\
         'fields >= INT. Default: INT = 40. --filter-tags TAG is used to filter '\
         'lines of find_circ.py output (sites.bed). Repeat it if multiple consecutive '\
         'filter tags has to be applied.', 
         '--best-qual 40 --filter-tags UNAMBIGUOUS_BP --filter-tags ANCHOR_UNIQUE')
## CircRNA_finder parameters
vars.Add('CFINDER_EXTRA_PARAMS', 
         'Parameters for CircRNA_finder', 
         '')

## alternative workflow parameters
vars.Add('PREPROCESSOR', 
         'The preprocessing method', 
         'trimmomatic')
vars.Add('PREPROCESSOR_PARAMS',
         'Read preprocessor extra parameters. F.i. if Trimmomatic, an empty string '\
         'defaults to '\
         'MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30 ',
         '')
vars.Add('LINEAR_EXPRESSION_METHODS', 
         'The method to be used for the linear expression estimates'\
         '/transcriptome reconstruction. To run more methods use a comma separated list. '\
    	 'However, only the first method in the list will be used in downstream processing. '\
	 'Currently supported methods: stringtie,cufflinks,htseq.', 
         'stringtie')
vars.Add('TOGGLE_TRANSCRIPTOME_RECONSTRUCTION', 
         'Set True to enable transcriptome '\
         'reconstruction. Default only quantifies genes and transcripts from the given '\
         'annotation GTF file', 
         'False')
#vars.Add('DIFF_EXP', 'Set the method to and enable differential expression '\
#         'computation for linear genes/transcripts. Current methods '\
#         'supported: cufflinks, ballgown, DESeq2.'\
#     	 'Only available if more than one sample and more than one condition are given. '\
#         'N.B: differential expression tests for circRNAs is not yet implemented', 
#         '')
vars.Add('READSTAT_METHODS', 
         'Comma separated list of methods to use for read statistics. '\
         'Currently supported: fastqc', 
         'fastqc')
vars.Add('MIN_METHODS', 
         'Number of methods that commmonly detect a circRNA to '\
         'define the circRNA as reliable. If this value exceeds the number '\
         'of methods specified, it will be set to the number of methods.', 
         2)
vars.Add('MIN_READS', 
         'Number of reads to consider a circRNA as expressed', 
         2)
vars.Add('BYPASS', 
         'Skip analysis of linear/circular transcripts. This will also skip '\
	 'the analysis of linear-to-circular expression correlation. The circular '\
         'analysis includes the pre-filtering of linearly mapping reads. If you '\
         'want to analyze reads already filtered for linear mappings you should '\
         'set "linear,linmap". Choose among linear and or linmap, circular. '\
         'NB: you still have to set the --rna-strandness '\
         'parameter into the HISAT_EXTRA_PARAMS if you have stranded alignments/reads.',
         '')
vars.Add('LINMAPS', 
         'You can specify here the path to pre-computed files of '\
         'linearly aligned reads. This will skip read pre-processing and linear '\
         'alignments (use jointly to BYPASS linmap to get also circular-to-linear analysis). '\
         'Mind that the alignments must be in BAM format and the .bai mapping file index '\
         'file must be in the same directory. NB: you still have to set the --rna-strandness '\
         'parameter into the HISAT_EXTRA_PARAMS if you have stranded alignments/reads. '\
         'You need to set a Python dict-like string parameter with sampleName and the '\
         'corresponding BAM file. E.g: {"SAMPLE1": "sample1/hisat2.bam", "SAMPLE2": "sample2/hisat2.bam"}', 
         None)
vars.Add('CIRC_MAPPING', 
         '''By default (SE), linearly unmapped reads are'''\
         '''aligned as single-end reads to search for circRNA backsplices. Set PE '''\
         '''to align as paired-end reads by each circRNA method aligner. You can also '''\
         '''specify each aligner's mode, or just which aligner has to use the PE mode, '''\
         ''' with the syntax for Python dictionaries {'SE':['ALN1','ALN2'],'PE':['ALN3','ALN4','ALNn']} '''\
         '''or simply {'PE':['ALN1','ALN2']} if you want just ALN1 and ALN2 tu align as PE. '''\
         '''Supported aligners are BWA,SEGEMEHL,STAR and TOPHAT. BOWTIE2 is also supported but '''\
         '''it is run only in single-end mode as it serves only Findcirc. ''',
         '''{'SE':['STAR','TOPHAT','BOWTIE2'],'PE':['BWA','SEGEMEHL']}''')

vars.Add('LIN_COUNTER', 
         'The method to estimate circRNA-host gene '\
         'linear expression. Available are using the DCC '\
         '[dcc], or the CirComPara [ccp] method', 
         'ccp')
vars.Add('FIX_READ_HEADER', 
         'Trim FASTQ headers to the read ids. '\
         'Recommended when processing SRA datasets', 
         'True')
vars.Add('UNSTRANDED_CIRCS', 
         'Force unstranded circRNAs even if stranded library was used', 
         'False')

## performance options
vars.Add('SAM_SORT_MM', 
         'Value for samtools sort -m option', 
         '768M')

## experimental
#vars.Add('CIRC_DIFF_EXP', '(Experimental) Set True to enable differential expression computation '\
#                          'for circRNAs. Only available if >1 conditions are given.', 
#         'False')
vars.Add('QRE_FIND', 
         '(Experimental) Set True to toggle analysis of QKI response elements sequences', 
         'False')
vars.Add('CCP_COUNTS', 
         'Use the CirComPara merged alignment counts', 
         'False')

## deprecated parameters (legacy)
#vars.Add('CIRI', '(DEPRECATED) The full path to the CIRI_vx.x.pl perl script. '\
#                 'By default, the symlink in CirComPara bin/ directory will be used', 
#         '')


env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                  variables=vars)
Help(vars.GenerateHelpText(env))
unknown = vars.UnknownVariables()
if unknown:
    print "Run sample: unknown variables", unknown.keys()
    Exit(1)

env.Append(BUILDERS = {'WriteLinesInTxt' : WriteLinesInTxt})

env['CIRCOMPARA_HOME'] = env['ENV']['CIRCOMPARA_HOME']
env['VARS'] = File(variables_file)
env['SCONSCRIPT_HOME'] = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')
env['CIRI'] = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'bin', 'CIRI.pl')
env['CCP_RMD_DIR'] = os.path.join(env['CIRCOMPARA_HOME'], 'src', 'utils', 'Rmd')

env.SetDefault(TOPHAT_PARAMS = '')

## the CIRC_MAPPING parameter is handled in the ccp_circrna_methods.py script

if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION'].lower() == 'true':
	env.Replace(TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = True)
else:
    env.Replace(TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = False)

if env['FIX_READ_HEADER'].lower() == 'true':
    env.Replace(FIX_READ_HEADER = True)
else:
    env.Replace(FIX_READ_HEADER = False)

env.SetDefault(UNSTRANDED_CIRCS = False)
if env['UNSTRANDED_CIRCS'].lower() == 'true':
    env.Replace(UNSTRANDED_CIRCS = True)
else:
    env.Replace(UNSTRANDED_CIRCS = False)

env['BYPASS'] = env['BYPASS'].split(',')

env.SetDefault(LINMAPS = None)
if env['LINMAPS']:
    try:
        env.Replace(LINMAPS = defaultdict(set, ast.literal_eval(str(env['LINMAPS']))))
    except ValueError as e:
        print e
        print '''Malformed LINMAPS string. Please, set as a Python dictionary, '''\
                '''e.g. {'SAMPLE1': 'path/to/linmaps.bam'}'''
        exit(-1)

if not 'circular' in env['BYPASS']:
    if env['CIRCRNA_METHODS'] == ['']:
        env.Replace(CIRCRNA_METHODS = ['ciri', 
                                       'circexplorer2_bwa', 
                                       'circexplorer2_segemehl',
                                       'circexplorer2_star', 
                                       'circexplorer2_tophat', 
                                       'findcirc', 
                                       'testrealign',
                                       'dcc',
                                       'circrna_finder'])
    env.Replace(CIRCRNA_METHODS = sorted(set([m.lower() for m in env['CIRCRNA_METHODS'].strip().split(',')])))
else:
    env['CIRCRNA_METHODS'] == ['']

env.Replace(LINEAR_EXPRESSION_METHODS = env['LINEAR_EXPRESSION_METHODS'].strip().split(','))

#if 'circexplorer2_tophat' in CIRCRNA_METHODS:
#	env.AppendUnique(TOPHAT_PARAMS = ['--bowtie1'])

env.SetDefault(DIFF_EXP = False)
#if env['DIFF_EXP'].strip() == '':
#    env.Replace(DIFF_EXP = False)

env.SetDefault(CIRC_DIFF_EXP = False)
#if env['CIRC_DIFF_EXP'].lower() == 'true':
#    env.Replace(CIRC_DIFF_EXP = True)
#else:
#    env.Replace(CIRC_DIFF_EXP = False)

env.SetDefault(CCP_COUNTS = False)
if env['CCP_COUNTS'].strip().lower() == 'true':
    env.Replace(CCP_COUNTS = True)
else:
    env.Replace(CCP_COUNTS = False)
    

env.SetDefault(ORIGINAL_ANNOTATION = env['ANNOTATION'])

## convert relative paths to absolute
if not env['ANNOTATION'] == '':
    env.Replace(ANNOTATION = File(env['ANNOTATION']).abspath)

if not env['GENOME_FASTA'] == '':
    env.Replace(GENOME_FASTA = File(env['GENOME_FASTA']).abspath)

if not env['GENOME_INDEX'] == '':
    env.Replace(GENOME_INDEX = File(env['GENOME_INDEX']).abspath)

if not env['SEGEMEHL_INDEX'] == '':
    env.Replace(SEGEMEHL_INDEX = File(env['SEGEMEHL_INDEX']).abspath)

if not env['BWA_INDEX'] == '':
    env.Replace(BWA_INDEX = File(env['BWA_INDEX']).abspath)

if not env['BOWTIE2_INDEX'] == '':
    env.Replace(BOWTIE2_INDEX = File(env['BOWTIE2_INDEX']).abspath)

if not env['STAR_INDEX'] == '':
    env.Replace(STAR_INDEX = Dir(env['STAR_INDEX']).abspath)

if not env['BOWTIE_INDEX'] == '':
    env.Replace(BOWTIE_INDEX = File(env['BOWTIE_INDEX']).abspath)

if not env['GENEPRED'] == '':
    env.Replace(GENEPRED = File(env['GENEPRED']).abspath)


samples_dir = 'samples'

## GRUB METADATA
env.Replace(META = File(env['META']).abspath)
samples     = defaultdict(list)
adapters    = defaultdict(str)
conditions  = defaultdict(set)
translocations  = set()

with open(env['META']) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        samples[row['sample']].append(os.path.abspath(row['file']))
        try:
            adapters[row['sample']] = row['adapter'].strip()
            ##NB:last sample adapter row overwrites the previous
        except KeyError as ke:
            #adapters[row['sample']] = ''
            print(str(ke) + ' not defined for file ' + row['file'] + ', sample ' +\
                    row['sample'])
            pass
        try:
            conditions[row['condition']].add(row['sample'])
        except KeyError as ke:
            print(str(ke) + ' not defined for file ' + row['file'] + ', sample ' +\
                    row['sample'])
            pass
        try:
            for tr in row['translocation'].split('#'):
                if tr: translocations.add(tr)
        except KeyError as ke:
            #print str(ke) + ' not defined for sample ' + row['sample'] + '. Skipping.'
            pass

env['CONDITIONS'] = conditions

#if len(conditions.keys()) < 2:
#    print str(conditions.keys()) + " is the only "\
#          "condition set for the samples. Differential "\
#          "expression analysis will not be performed."
#    env.Replace(DIFF_EXP = False)

## PREPARE GENOME AND ANNOTATION FOR TRANSLOCATED SAMPLES
#TODO in F-CirComPara

## BUILD READ ALIGNER PROGRAM GENOME INDEXES IF NOT PROVIDED BY THE USER
indexes_dir = 'ccp_dbs'
env_check_indexes = env.Clone()
indexes = SConscript(os.path.join(indexes_dir, 'ccp_check_indexes.py'),
                        src_dir = env['SCONSCRIPT_HOME'],
                        variant_dir = indexes_dir, duplicate = 0,
                        exports = '''env_check_indexes''')

env.Replace(GENOME_INDEX   = indexes['ENV']['GENOME_INDEX'])
env.Replace(SEGEMEHL_INDEX = indexes['ENV']['SEGEMEHL_INDEX'])
env.Replace(BWA_INDEX      = indexes['ENV']['BWA_INDEX'])
env.Replace(BOWTIE2_INDEX  = indexes['ENV']['BOWTIE2_INDEX'])
env.Replace(BOWTIE_INDEX   = indexes['ENV']['BOWTIE_INDEX'])
env.Replace(STAR_INDEX     = indexes['ENV']['STAR_INDEX'])
env.Replace(GENEPRED       = indexes['ENV']['GENEPRED'])

## PROCESS SAMPLES
runs = []
runs_dict = {}
for sample in sorted(samples.keys()):
    
    env_circpipe = env.Clone()
    env_circpipe['SAMPLE'] = sample
    env_circpipe['READS'] = [File(f) for f in samples[sample]]
    env_circpipe['ADAPTER_FILE'] = adapters[sample]
    if env['LINMAPS']:
        env_circpipe['LINMAPS'] = File(env['LINMAPS'][sample])
    if len(env_circpipe['READS']) > 1:
        env_circpipe.Replace(CIRCRNA_METHODS = [re.sub(r'\bcircexplorer2_tophat\b', 
                                                        'circexplorer2_tophat_pe', m) for \
                                                        m in env_circpipe['CIRCRNA_METHODS']])
        
    sample_dir = os.path.join(samples_dir, sample)
    run_sample = SConscript(os.path.join(sample_dir, 'ccp_circpipe.py'),
                            src_dir = env['SCONSCRIPT_HOME'],
                            variant_dir = sample_dir, duplicate = 0,
                            exports = '''env_circpipe get_matching_nodes''')
    runs.append(run_sample[0])
    runs_dict[sample] = run_sample[1]

    ## dependencies from aligners' genome index
    ## TODO: set these dependencies in each Sconscript command by set
    ## proper command sources
    Depends([f for f in runs_dict[sample]['LINEAR_MAPPING_RES'].values() if f], 
            indexes['INDEXES']['HISAT2'])

    if not 'circular' in env['BYPASS']:
        if runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['SEGEMEHL_MAP']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['SEGEMEHL_MAP'].values(), 
                    indexes['INDEXES']['SEGEMEHL']) 

        if runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['STAR_MAP']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['STAR_MAP'], 
                    indexes['INDEXES']['STAR'])

        if runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['BWA_MAP']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['BWA_MAP'], 
                    indexes['INDEXES']['BWA'])

        if runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['TOPHAT_MAP']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_ALIGNERS']['TOPHAT_MAP'].values(), 
                    indexes['INDEXES']['BOWTIE'])

        ## dependencies from processed gene annotation
        if 'circexplorer2_bwa' in env['CIRCRNA_METHODS']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_BWA'].values(),
                    env['GENEPRED'])

        if 'circexplorer2_segemehl' in env['CIRCRNA_METHODS']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_SE'].values(),
                    env['GENEPRED'])

        if 'circexplorer2_star' in env['CIRCRNA_METHODS']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_ST'].values(),
                    env['GENEPRED'])

        if'circexplorer2_tophat' in env['CIRCRNA_METHODS']:
            Depends(runs_dict[sample]['CIRCULAR_EXPRESSION']['CIRC_METHODS']['CE2_TH'].values(),
                    env['GENEPRED'])

env['RUNS_DICT'] = runs_dict

## LINEAR EXPRESSION: 
## - TRANSCRIPTOME RECONSTRUCTION, 
## - NORMALIZED QUANTIFICATION
## - DIFFERENTIAL EXPRESSION
linexp_dir = 'linear_expression'
if not 'linear' in env['BYPASS']:
    env_linear_expression = env.Clone()
    env_linear_expression['SAMPLES'] = samples
    env_linear_expression['RUNS'] = runs
    env_linear_expression['RUNS_DICT'] = runs_dict
    linexp = env.SConscript(os.path.join(linexp_dir,
                                         'ccp_linear_expression.py'),
                            src_dir = env['SCONSCRIPT_HOME'],
                            variant_dir = linexp_dir, duplicate = 0,
                            exports = '''env_linear_expression '''
                                      '''SymLink '''
                                      '''get_matching_nodes''')

    env.Replace(ANNOTATION = linexp['ANNOTATION'])
else:
    ## DO NOT ANALYZE GENE EXPRESSION: USE MOCK EMPTY FILES
    print "BYPASS = " + str(env['BYPASS']) +\
          ": skipping linear transcript analysis"
    
    gene_exp = env.Command(os.path.join(linexp_dir, 'empty_geneexp.csv'), 
    			'', 'touch $TARGET')
    linexp = {'GENE_EXP_ANALYSIS': ['', gene_exp]}

## COLLECT AND REPORT READ PROCESSING STATISTICS
collect_read_stats_dir = 'read_statistics'
env_collect_read_stats = env.Clone()
env_collect_read_stats['RUNS'] = runs
collect_read_stats = env.SConscript(os.path.join(collect_read_stats_dir, 
                                                'ccp_collect_read_stats.py'),
                                    src_dir = env['SCONSCRIPT_HOME'],
                                    variant_dir = collect_read_stats_dir,
                                    duplicate = 0, 
                                    exports = '''env_collect_read_stats '''
                                              '''SymLink '''
                                              '''get_matching_nodes''')

### COLLECT AND REPORT CIRCRNA RESULTS
circexp_dir = 'circular_expression'
if not 'circular' in env['BYPASS']:
    env_circular_expression = env.Clone()
    env_circular_expression['SAMPLES'] = samples
    env_circular_expression['RUNS'] = runs
    env_circular_expression['RUNS_DICT'] = runs_dict
    env_circular_expression['LINEXP'] = linexp
    env_circular_expression['PROCESSING_READ_STATS'] = collect_read_stats['PROCESSING_READ_STATS']
    circexp = env.SConscript(os.path.join(circexp_dir,
                                          'ccp_circular_expression.py'),
                             src_dir = env['SCONSCRIPT_HOME'],
                             variant_dir = circexp_dir, duplicate = 0,
                             exports = '''env_circular_expression '''
                                       '''SymLink '''
                                       '''get_matching_nodes''')

## ADDITIONAL ANALYSIS
qre_dir = 'qre'
if env['QRE_FIND'] == 'True':
    ## ANALYZE GENE SEQUENCES FOR QKI RESPONSE ELEMENTS
    qre_GTF = env['ANNOTATION'] #cuffmerge
    qre_GENOME = env['GENOME_FASTA']
    qre = SConscript(os.path.join(qre_dir, 'ccp_QRE_finder.py'),
                     src_dir = env['SCONSCRIPT_HOME'],
                     variant_dir = qre_dir, duplicate = 0,
                     exports = '''env qre_GTF qre_GENOME''')
    
    Depends(qre, [env['ANNOTATION'], env['GENOME_FASTA']])


## CLEAN DIRS WHEN CLEANING TARGETS
Clean('.', indexes_dir)
Clean('.', samples_dir)
Clean('.', linexp_dir)
Clean('.', circexp_dir)
Clean('.', collect_read_stats_dir)
Clean('.', qre_dir)
