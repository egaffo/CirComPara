'''
Find the Quaking Response Element (QRE) sequence for given genes.
Genes must be defined by a GTF (or a BED). 
Requires the genome sequence FASTA.

1. Galarneau, A. & Richard, S. Target RNA motif and target mRNAs of 
the Quaking STAR protein. Nat Struct Mol Biol 12, 691-698 (2005).

Software dependencies:
 * gtf_collapse_features.py 
 * QRE_finder.py

Returns:
    qre.tsv
'''

import os

Import('*')

try:
    env
    GTF         = qre_GTF
    GENOME      = qre_GENOME
except NameError:

    vars = Variables('vars.py')
    vars.Add('GTF', '''The annotation GTF''', 'merged.gtf')
    vars.Add('GENOME', '''The genome multi FASTA''', 'hg38.fa')
   
    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

    GTF     = env['GTF']
    GENOME  = env['GENOME']

qre_cmd = '''gtf_collapse_features.py ${SOURCES[0]} | '''\
          '''bedtools getfasta -name -s -fi ${SOURCES[1]} -bed - -fo stdout | '''\
          '''QRE_finder.py - > $TARGET'''
qre_target = 'qre.tsv'
qre = env.Command(qre_target, [GTF, GENOME], qre_cmd)

Return('qre')
