'''
This SConscript performs quality selection of RNA-seq FASTQ reads.

Required variables to export when calling from a SConscript:
 * env
 * preprocess_cpus
 * preprocess_preprocessor
 * preprocess_adapter_file
 * preprocess_raw_reads
 * preprocess_params

'''
import os, itertools

Import('*')

try:
    env = env_preprocess.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('PREPROCESSOR', 'The preprocessing method', 'trimmomatic')
    vars.Add('ADAPTER_FILE', 'FASTA file full path of the adapter sequence', '')
    vars.Add('READS', 'RNA-seq reads. Comma separated list if paired-end', 'reads.fa')
    vars.Add('PREPROCESSOR_PARAMS', 
             '''Read preprocessor extra parameters. F.i. if Trimmomatic, an empty string '''\
             '''defaults to '''\
             '''MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30 ''', 
             '')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src')

preprocess_dir = 'preprocess'
preprocessings = []
preprocessing_stats = []
if env['PREPROCESSOR'].lower() == 'trimmomatic':

    trimmomatic_reads = [File(f).abspath for f in env['READS']]
    trimmomatic_cpus = env['CPUS']
    trimmomatic_adapter_file = env['ADAPTER_FILE']
    trimmomatic_extra_params = env['PREPROCESSOR_PARAMS']

    preps = SConscript(os.path.join(preprocess_dir, 'ccp_trimmomatic.py'),
                       src_dir = SRC_DIR, 
                       variant_dir = preprocess_dir, duplicate = 0, 
                       exports = '''env trimmomatic_reads trimmomatic_adapter_file '''
                                 '''trimmomatic_extra_params trimmomatic_cpus''')
    preprocessings.append(preps)

    ## COMPUTE STATISTICS ON PREPROCESSED READS
    files2skip = ['.log']
    for f in itertools.chain(*preprocessings):
        env_read_statistics = env.Clone()
        reads_file_path = File(f).path
        ext = os.path.splitext(os.path.basename(reads_file_path))[1]
        if ext in files2skip: continue
        #read_statistics_readset = reads_file_path
        env_read_statistics['READSET'] = reads_file_path
        #if not os.path.isabs(read_statistics_readset):
            #read_statistics_readset = '#'+read_statistics_readset
        if not os.path.isabs(env_read_statistics['READSET']):
            env_read_statistics['READSET'] = '#' + env_read_statistics['READSET']
        preprocessing_stats.append(SConscript(os.path.join(preprocess_dir, 
                                                           'ccp_read_statistics.py'), 
                                              src_dir = SRC_DIR, 
                                              variant_dir = preprocess_dir, duplicate = 0, 
                                              exports = 'env_read_statistics')
                                  )

    Depends(preprocessing_stats, preprocessings)

elif env['PREPROCESSOR'] == '':
    ## in case of no preprocessor specified, just return the raw reads
    print 'No read preprocessing specified'
    preprocessings = [[File(f) for f in env['READS']]]
else:
    ## in case of preprocessor not supported, just return the raw reads
    print env['PREPROCESSOR'] + ' preprocessor not supported. Reads will not be preprocessed'
    preprocessings.append([File(f) for f in env['READS']])

preprocessed_files = {'READS' : preprocessings,
                      'STATS' : preprocessing_stats}

Clean('.', preprocess_dir)
Return('preprocessed_files')
