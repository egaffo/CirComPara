import os 

Import('*')

try:
    env = env_collect_read_stats.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('', '', '')

    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

runs = env['RUNS']

## COLLECT AND REPORT READ PROCESSING STATISTICS
## case: Trimmomatic preprocessing
clean_reads_stats_files = get_matching_nodes(runs, '.*preprocess.*trimmomatic\.log')
if len(clean_reads_stats_files) > 0:
    reads_stats_collect_cmd = '''grep --with-filename Input ''' +\
                              ' '.join(f.path for f in clean_reads_stats_files)

else:
    ## case: no read preprocessing
    clean_reads_stats_files = get_matching_nodes(runs,
                                                 '.*read_statistics.*fastqc_data\.txt')
    reads_stats_collect_cmd = '''grep --with-filename "Total Sequences" ''' +\
                              ' '.join(f.path for f in clean_reads_stats_files)

reads_stats_collect_cmd = reads_stats_collect_cmd + ''' > $TARGET '''

read_stats_collect_sources = [clean_reads_stats_files]

if env['LINMAPS']:
    mapped_reads_linear_stats_files = get_matching_nodes(runs, '.*processings.*hisat2\.log')
    if mapped_reads_linear_stats_files:
        reads_stats_collect_cmd = reads_stats_collect_cmd +\
                                  ''' && grep --with-filename "." ''' +\
                                  ' '.join(f.path for f in mapped_reads_linear_stats_files) +\
                                  ''' >> $TARGET'''
        read_stats_collect_sources += mapped_reads_linear_stats_files

read_stats_collect_dir = 'read_stats_collect'
reads_stats_collect = env.Command(os.path.join(read_stats_collect_dir, 'read_stats_collect.txt'), 
                                  read_stats_collect_sources,
                                  reads_stats_collect_cmd)


##COLLECT CIRCRNA READS
# testrealign reads
testrealign_pattern = '.*.testrealign.bks.reads' #'.*segemehl_mapped_reads_count\.txt'
# find_circ reads
find_circ_pattern = '.*.fc.bks.reads' #'.*find_circ_mapped_reads_count\.txt'
# CIRI reads
CIRI_pattern = '.*.ciri.bks.reads' #'.*CIRI_mapped_reads_count\.txt'
# CIRCexplorer2 bwa reads
CIRCexplorerBWA_pattern = '.*.ce2bwa.bks.reads' #'.*CIRCExplorer_mapped_reads_count\.txt'
# CIRCexplorer2 segemehl reads
CIRCexplorerSEG_pattern = '.*.ce2seg.bks.reads'
# CIRCexplorer2 star reads
CIRCexplorerSTR_pattern = '.*.ce2star.bks.reads'
# CIRCexplorer2 tophat reads
CIRCexplorerTPH_pattern = '.*.ce2th.bks.reads'
# DCC reads
dcc_pattern = '.*.dcc.bks.reads'
# CircRNA_finder reads
cfinder_pattern = '.*.cfinder.bks.reads'


## testrealign mapped reads
#testrealign_mappings_pattern = '.*segemehl_mapped_reads_count\.txt'
## BWA mapped reads
#BWA_mappings_pattern = '.*BWA_mapped_reads_count\.txt'
## STAR mapped reads
#STAR_mappings_pattern = '.*STAR_mapped_reads_count\.txt'
## TopHat mapped reads
#tophat_mappings_pattern = '.*tophat_mapped_reads_count\.txt'


testrealign_reads  = get_matching_nodes(runs, testrealign_pattern)
find_circ_reads    = get_matching_nodes(runs, find_circ_pattern)
ciri_reads         = get_matching_nodes(runs, CIRI_pattern)
CIRCexplorerBWA_reads = get_matching_nodes(runs, CIRCexplorerBWA_pattern)
CIRCexplorerSEG_reads = get_matching_nodes(runs, CIRCexplorerSEG_pattern)
CIRCexplorerSTR_reads = get_matching_nodes(runs, CIRCexplorerSTR_pattern)
CIRCexplorerTPH_reads = get_matching_nodes(runs, CIRCexplorerTPH_pattern)
dcc_reads          = get_matching_nodes(runs, dcc_pattern)
cfinder_reads      = get_matching_nodes(runs, cfinder_pattern)

#bwa_mappings          = get_matching_nodes(runs, BWA_mappings_pattern)
#star_mappings         = get_matching_nodes(runs, STAR_mappings_pattern)
#tophat_mappings       = get_matching_nodes(runs, tophat_mappings_pattern)

collect_circrna_maps_counts_sources = [testrealign_reads, 
                                       find_circ_reads, 
                                       ciri_reads, 
                                       CIRCexplorerBWA_reads,
                                       CIRCexplorerSEG_reads,
                                       CIRCexplorerSTR_reads,
                                       CIRCexplorerTPH_reads,
                                       dcc_reads,
                                       cfinder_reads]
				       #bwa_mappings, star_mappings, 
				       #tophat_mappings]

if len(set(Flatten(collect_circrna_maps_counts_sources))) > 0:
    #collect_circrna_maps_counts_cmd = '''tail -n +1 ${SOURCES} > $TARGET '''
    collect_circrna_maps_counts_cmd = '''grep -H "." ${SOURCES} | '''\
                                      '''sed "s/.bks.reads:/\\t/" | gzip -c > $TARGET '''
else:
    collect_circrna_maps_counts_cmd = '''touch $TARGET '''

collect_circrna_maps_counts = env.Command(os.path.join(read_stats_collect_dir, 
                                                       'circrna_maps_counts.txt.gz'), 
                                          collect_circrna_maps_counts_sources, 
                                          collect_circrna_maps_counts_cmd)

if env['BYPASS'] == 'circular':
   env['CIRCRNA_RMD'] = "_empty.Rmd"
else:
   env['CIRCRNA_RMD'] = "_circrna_read_stats.Rmd"

## report read processing statistics
read_stats_report_cmd = '''Rscript -e 'results.dir <- dirname("$TARGET.abspath"); '''\
                        '''read_stats_collect.file <- "${SOURCES[0].abspath}"; '''\
                        '''linear.mapper <- "hisat2"; '''\
                        '''circrna.reads.stats.file <- "${SOURCES[1].abspath}"; '''\
                        '''meta_file <- "${SOURCES[2].abspath}"; '''\
                        '''circrna_read_stats.Rmd <- file.path("$CCP_RMD_DIR", "$CIRCRNA_RMD"); '''\
                        '''rmarkdown::render(input = "$CCP_RMD_DIR/read_statistics.Rmd",'''\
                        '''output_file = "$TARGET.abspath", quiet=T,'''\
                        '''intermediates_dir = dirname("$TARGET.abspath") )' '''

read_stats_report = env.Command([os.path.join(read_stats_collect_dir, f) for 
                                  f in ['read_statistics.html', 
                                        'processing_and_mapped_read_counts.csv']],
                                [reads_stats_collect, collect_circrna_maps_counts, 
                                 File(env['META']).abspath, 
                                 collect_circrna_maps_counts_sources],
                                read_stats_report_cmd)


Clean('.', read_stats_collect_dir)

results = {'READ_MAPS_COUNTS': collect_circrna_maps_counts, 
           'READS_STATS_REPORT': read_stats_report[0],
           'PROCESSING_READ_STATS': read_stats_report[1]}

Return('results')
