import os

def SymLink(target, source, env):
    os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

env = Environment(ENV=os.environ, SHELL = '/bin/bash')
tools_dir = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'tools')
ccp_bin_dir = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'bin')
python_lib_dir = os.path.join(tools_dir, 'lib', 'python2.7','site-packages')

#env.PrependENVPath('PYTHONUSERBASE', tools_dir)
env.PrependENVPath('PYTHONUSERBASE', tools_dir)
env.PrependENVPath('PYTHONPATH', python_lib_dir)

## PIP
pip_file = 'get-pip.py'
pip_url = 'https://bootstrap.pypa.io/' + pip_file
pip_targets = [os.path.join(tools_dir, pip_file),
               os.path.join(tools_dir, 'bin', 'pip')]
pip_cmd = ' && '.join(['wget -O ${TARGETS[0]} ' + pip_url, 
                       'python ${TARGETS[0]} --user'])
pip = env.Command(pip_targets, 
                  [], 
                  pip_cmd)


# BIOPYTHON
BIOPYTHON_dir = os.path.join(python_lib_dir, 'Bio')
BIOPYTHON_target = [os.path.join(BIOPYTHON_dir, 'SeqIO', 'FastaIO.py')]
BIOPYTHON = env.Command(BIOPYTHON_target, [pip], 
                        #['pip install --install-option="--prefix=' +\
                        #os.path.abspath(BIOPYTHON_dir) +\
                        #'" biopython'])
                        ['pip install --ignore-installed --user biopython'])

# HTSeq
#HTSeq_dir = os.path.join(tools_dir, 'HTSeq')
HTSeq_dir = os.path.join(python_lib_dir, 'HTSeq')
#HTSeq_target = [os.path.join(HTSeq_dir, 'lib','python2.7','site-packages',
#                'HTSeq', '__init__.py'),
HTSeq_target = [os.path.join(HTSeq_dir, '__init__.py'),
                os.path.join(tools_dir, 'bin', 'htseq-count')]
HTSeq = env.Command(HTSeq_target, [pip], 
                           #['pip install --ignore-installed --install-option="--prefix=' +\
                           #os.path.abspath(HTSeq_dir) +\
                           #'" HTSeq'])
                    ['pip install --ignore-installed --user HTSeq'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), HTSeq[1], SymLink)

# CIRCEXPLORER2
#CIRCEXPLORER2_dir = os.path.join(tools_dir, 'CIRCexplorer2')
CIRCEXPLORER2_dir = os.path.join(python_lib_dir, 'circ2')
CIRCEXPLORER2_target = [os.path.join(tools_dir, 'bin', 'CIRCexplorer2')]#, 
#                       os.path.join(CIRCEXPLORER2_dir, 'bin', 'fast_circ.py')]
CIRCEXPLORER2 = env.Command(CIRCEXPLORER2_target, [pip, HTSeq], 
                           #['pip install --install-option="--prefix=' +\
                           #os.path.abspath(CIRCEXPLORER2_dir) +\
                           #'" -Iv circexplorer2==2.3.3'])
                           ['pip install --ignore-installed --user -Iv circexplorer2==2.3.3'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), CIRCEXPLORER2[0], SymLink)
#env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), CIRCEXPLORER2[1], SymLink)

## DCC
#env_dcc = env.Clone()
#dcc_dir = os.path.join(tools_dir, 'DCC-0.4.6')
#dcc_python_sitepkg = os.path.join(dcc_dir, "lib", "python2.7", "site-packages")
#env_dcc.PrependENVPath('PYTHONPATH', dcc_python_sitepkg)

### PANDAS (DCC dependency)
#env_dcc.PrependENVPath('PYTHONUSERBASE', dcc_dir)

#pandas = env.Command([os.path.join(dcc_python_sitepkg, 'pandas', '__init__.py')], 
pandas_dir = os.path.join(python_lib_dir, 'pandas')
pandas = env.Command([os.path.join(pandas_dir, '__init__.py')], 
                         [pip, HTSeq], 
                         "pip install --ignore-installed --user pandas")

#dcc_dir = os.path.join(python_lib_dir, 'DCC-0.4.6')
dcc_dir = os.path.join(tools_dir, 'DCC-0.4.6')
dcc_tar = 'v0.4.6.tar.gz'
dcc_url = 'https://github.com/dieterich-lab/DCC/archive/' + dcc_tar
dcc_target = [os.path.join(tools_dir, dcc_tar),
		      os.path.join(tools_dir, 'bin', 'DCC'),
              #os.path.join(dcc_python_sitepkg, "site.py")]
              os.path.join(python_lib_dir, 'DCC-0.4.6-py2.7.egg')]
dcc = env.Command(dcc_target, 
			         [pandas], 
          	         ['wget -O ${TARGETS[0]} ' + dcc_url,
			          'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                      'cd ' + dcc_dir +\
                      #' && python setup.py install --prefix ' + dcc_dir,
                      ' && python setup.py install --user',
                      'cd ' + Dir('#').abspath]
			         )

dcc_link = env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), dcc[1], SymLink)

# HISAT2
HISAT2_zip = 'hisat2-2.0.4-Linux_x86_64.zip'
HISAT2_link = 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/' + HISAT2_zip
HISAT2_target = [os.path.join(tools_dir, HISAT2_zip), 
                 os.path.join(tools_dir, 'hisat2-2.0.4', 'hisat2'),
                 os.path.join(tools_dir, 'hisat2-2.0.4', 'hisat2-build')]
#HISAT2_source = HISAT2_link
HISAT2 = env.Command(HISAT2_target, [], ['wget -O ${TARGETS[0]} ' + HISAT2_link,
                                         'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), HISAT2[1], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), HISAT2[2], SymLink)

# BOWTIE2
BOWTIE2_zip = 'bowtie2-2.2.9-linux-x86_64.zip'
BOWTIE2_link = 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/' + BOWTIE2_zip
BOWTIE2_target = [os.path.join(tools_dir, BOWTIE2_zip), 
                  os.path.join(tools_dir, 'bowtie2-2.2.9', 'bowtie2'),
                  os.path.join(tools_dir, 'bowtie2-2.2.9', 'bowtie2-inspect'),
                  os.path.join(tools_dir, 'bowtie2-2.2.9', 'bowtie2-build')]
#BOWTIE2_source = 
BOWTIE2 = env.Command(BOWTIE2_target, [], ['wget -O ${TARGETS[0]} ' + BOWTIE2_link,
                                           'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE2[1], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE2[2], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE2[3], SymLink)

## BWA-MEM
BWAMEM_tar = 'bwakit-0.7.15_x64-linux.tar.bz2'
BWAMEM_link = 'https://sourceforge.net/projects/bio-bwa/files/bwakit/' + BWAMEM_tar
BWAMEM_target = [os.path.join(tools_dir, BWAMEM_tar),
                 os.path.join(tools_dir, 'bwa.kit', 'bwa')]
#BWAMEM_source = 
BWAMEM = env.Command(BWAMEM_target, [], ['wget -O ${TARGETS[0]}  ' + BWAMEM_link,
                                         'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BWAMEM[1], SymLink)

# STAR
STAR_tar = '2.6.1d.tar.gz'
STAR_link = 'https://github.com/alexdobin/STAR/archive/' + STAR_tar
STAR_target = [os.path.join(tools_dir, 'STAR_' + STAR_tar), 
               os.path.join(tools_dir, 'STAR-2.6.1d', 'bin', 'Linux_x86_64_static', 'STAR')]
#STAR_source = 
STAR = env.Command(STAR_target, [], ['wget -O ${TARGETS[0]} ' + STAR_link,
                                     'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), STAR[1], SymLink)

# HTSLIB (required by Segemehl)
HTSLIB_dir = os.path.join(tools_dir, 'htslib-1.9')
HTSLIB_tar = 'htslib-1.9.tar.bz2'
HTSLIB_target = [os.path.join(tools_dir, HTSLIB_tar),
                 os.path.join(HTSLIB_dir, 'lib', 'pkgconfig', 'htslib.pc')]
HTSLIB_link = 'https://github.com/samtools/htslib/releases/download/1.9/' +\
              HTSLIB_tar
HTSLIB = env.Command(HTSLIB_target,
                     [],
                     ['wget -O $TARGET  ' + HTSLIB_link,
                      'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                      ' && '.join(['cd ' + HTSLIB_dir, 
                                   './configure --prefix=`pwd`',
                                   'make', 'make install',
                                   'cd ' + Dir('.').abspath])
                      ])

# SEGEMEHL
SEGEMEHL_tar = 'segemehl-0.3.4.tar.gz'
SEGEMEHL_dir = os.path.join(tools_dir, 'segemehl-0.3.4')
SEGEMEHL_link = 'http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/' + SEGEMEHL_tar
SEGEMEHL_target = [os.path.join(tools_dir, SEGEMEHL_tar),
                   os.path.join(SEGEMEHL_dir, 'segemehl.x'),
                   os.path.join(SEGEMEHL_dir, 'haarz.x')]

env.PrependENVPath('PKG_CONFIG_PATH',
                   os.path.dirname(HTSLIB[1].abspath))
SEGEMEHL = env.Command(SEGEMEHL_target, 
                       HTSLIB, 
                       ['wget -O $TARGET  ' + SEGEMEHL_link,
                       'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                       'cd ${TARGETS[1].dir} && make all',
                       'cd ' + Dir('.').abspath])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), SEGEMEHL[1], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), SEGEMEHL[2], SymLink)

# TRIMMOMATIC
TRIMMOMATIC_zip = 'Trimmomatic-0.38.zip'
TRIMMOMATIC_link = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/' + TRIMMOMATIC_zip
TRIMMOMATIC_target = [os.path.join(tools_dir, TRIMMOMATIC_zip), 
                      os.path.join(tools_dir, 'Trimmomatic-0.38', 'trimmomatic-0.38.jar')]
#TRIMMOMATIC_source = 
TRIMMOMATIC = env.Command(TRIMMOMATIC_target, [], ['wget -O $TARGET  ' + TRIMMOMATIC_link,
                                                   'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), TRIMMOMATIC[1], SymLink)


# FASTQC
FASTQC_zip = 'fastqc_v0.11.5.zip'
FASTQC_link = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/' + FASTQC_zip
FASTQC_target = [os.path.join(tools_dir, FASTQC_zip), 
                 os.path.join(tools_dir, 'FastQC', 'fastqc')]
#FASTQC_source = 
FASTQC = env.Command(FASTQC_target, [], ['wget -O $TARGET  ' + FASTQC_link,
                                         'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]} '\
                                         '&& chmod +x ${TARGETS[1]}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), FASTQC[1], SymLink)


# BEDTOOLS
BEDTOOLS_tar = 'bedtools-2.26.0.tar.gz'
BEDTOOLS_dir = os.path.join(tools_dir, 'bedtools2')
BEDTOOLS_link = 'https://github.com/arq5x/bedtools2/releases/download/v2.26.0/' + BEDTOOLS_tar
BEDTOOLS_target = [os.path.join(tools_dir, BEDTOOLS_tar), 
                   os.path.join(BEDTOOLS_dir, 'bin', 'bedtools')]
#BEDTOOLS_source = 
BEDTOOLS = env.Command(BEDTOOLS_target, [], 
                       ['wget -O $TARGET ' + BEDTOOLS_link, 
                        'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}', 
                        'cd ' + BEDTOOLS_dir + ' && make', 
                        'cd ' + Dir('.').abspath])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BEDTOOLS[1], SymLink)

# SAMTOOLS
SAMTOOLS_tar = 'samtools-1.3.1.tar.bz2'
SAMTOOLS_dir = os.path.join(tools_dir, 'samtools-1.3.1')
SAMTOOLS_link = 'https://github.com/samtools/samtools/releases/download/1.3.1/' + SAMTOOLS_tar
SAMTOOLS_target = [os.path.join(tools_dir, 'samtools-1.3.1.tar.bz2'),
                   os.path.join(SAMTOOLS_dir, 'samtools')]
#SAMTOOLS_source = 
SAMTOOLS = env.Command(SAMTOOLS_target, [], ['wget -O $TARGET  ' + SAMTOOLS_link,
                                             'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                                             'cd ' + SAMTOOLS_dir + ' && '\
                                             'make prefix=' + SAMTOOLS_dir + ' install',
                                             'cd ' + Dir('.').abspath])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), SAMTOOLS[1], SymLink)

# CUFFLINKS
CUFFLINKS_tar = 'cufflinks-2.2.1.Linux_x86_64.tar.gz'
CUFFLINKS_link = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/' + CUFFLINKS_tar
CUFFLINKS_target = [os.path.join(tools_dir, CUFFLINKS_tar), 
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cufflinks'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffcompare'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffdiff'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffmerge'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffnorm'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffquant'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'gffread'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'gtf_to_sam')]
#CUFFLINKS_source = 
CUFFLINKS = env.Command(CUFFLINKS_target, [], ['wget -O $TARGET  ' + CUFFLINKS_link,
                                               'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
for t2link in CUFFLINKS[1:]:
    env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), t2link, SymLink)


# CIRI
CIRI_zip = 'CIRI_v2.0.2.zip'
CIRI_link = 'http://downloads.sourceforge.net/project/ciri/CIRI2/' + CIRI_zip
CIRI_target = [os.path.join(tools_dir, CIRI_zip), 
               os.path.join(tools_dir, 'CIRI_v2.0.2', 'CIRI_v2.0.2.pl')]
#CIRI_source = 
CIRI = env.Command(CIRI_target, [], ['wget -O ${TARGETS[0]} ' + CIRI_link,
                                     'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(ccp_bin_dir, "CIRI.pl"), CIRI[1], SymLink)


# FIND-CIRC
FINDCIRC_tar = 'find_circ.zip'
#FINDCIRC_link = 'http://www.circbase.org/download/' + FINDCIRC_tar
FINDCIRC_link = 'https://github.com/marvin-jens/find_circ/archive/master.zip'
FINDCIRC_target = [os.path.join(tools_dir, FINDCIRC_tar),
                   os.path.join(tools_dir, 'find_circ-master', 'find_circ.py'),
		   os.path.join(tools_dir, 'find_circ-master', 'unmapped2anchors.py')]
#FINDCIRC_source = 
FINDCIRC = env.Command(FINDCIRC_target, [], ['wget -O $TARGET  ' + FINDCIRC_link,
                                             'unzip -u -d ' + tools_dir + ' ${TARGETS[0]}'])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), FINDCIRC[1], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), FINDCIRC[2], SymLink)

# FASTX-TOOLKIT
libgtextutils_tar   = 'libgtextutils-0.7.tar.gz'
libgtextutils_link  = 'https://github.com/agordon/libgtextutils/releases/download/0.7/' +\
                     libgtextutils_tar
libgtextutils_dir   = os.path.join(tools_dir, 'libgtextutils-0.7')
libgtextutils_target = [os.path.join(tools_dir, libgtextutils_tar), 
                        os.path.join(libgtextutils_dir, 'lib', 'libgtextutils.so'),
                        os.path.join(libgtextutils_dir, 'lib', 'pkgconfig', 'gtextutils.pc')]
libgtextutils = env.Command(libgtextutils_target, [], 
                            ['wget -O $TARGET ' + libgtextutils_link, 
                            'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}', 
                            'cd ' + libgtextutils_dir + ' && ./configure --prefix=`pwd`'\
                            ' && make && make install',
                            'cd ' + Dir('#').abspath])

libgtextutils_config_path = os.path.dirname(libgtextutils[2].abspath)

FASTXTOOLKIT_tar = 'fastx_toolkit-0.0.14.tar.bz2'
FASTXTOOLKIT_link = 'https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/' +\
                    FASTXTOOLKIT_tar
FASTXTOOLKIT_dir = os.path.join(tools_dir, 'fastx_toolkit-0.0.14')
FASTXTOOLKIT_target = [os.path.join(tools_dir, FASTXTOOLKIT_tar),
                       os.path.join(FASTXTOOLKIT_dir, 'bin', 'fastx_quality_stats')]
FASTXTOOLKIT = env.Command(FASTXTOOLKIT_target, [Value(libgtextutils_config_path), libgtextutils[2]], 
                           ['wget -O $TARGET ' + FASTXTOOLKIT_link,
                            'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                            'cd ' + FASTXTOOLKIT_dir + ' && ./configure --prefix=' +\
                            FASTXTOOLKIT_dir + ' PKG_CONFIG_PATH=${SOURCES[0]}'\
                            ':$$PKG_CONFIG_PATH '\
                            ' && make && make install', 
                            'cd ' + Dir('#').abspath])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), FASTXTOOLKIT[1], SymLink)

# gtfToGenePred
gtfToGenePred_link = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred'
gtfToGenePred_target = [os.path.join(tools_dir, 'gtfToGenePred')]
gtfToGenePred = env.Command(gtfToGenePred_target, [], 
                            ['wget -O $TARGET ' + gtfToGenePred_link, 
                            Chmod('$TARGET', 0775)])
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), gtfToGenePred, SymLink)

# optparse, ggplot2, DATA.TABLE, plyr, scales, reshape2, ggthemes, RSvgDevice
# Bioconductor: ReportingTools, DESeq2
env['ENV']['R_LIBS'] = os.path.join(tools_dir, "R_libs")
R_libs_targets = [os.path.join(tools_dir, 'R_libs', 'DESeq2', 'R', 'DESeq2')]
R_libs = env.Command(R_libs_targets, [], 'install_R_libs.R')

# BOWTIE v1
# v1.2.1.1, v1.2.1, and v1.2 do not work!
#BOWTIE1_zip = 'bowtie-1.2.1.1-linux-x86_64.zip'
#BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
#		'1.2.1.1/' + BOWTIE1_zip
#BOWTIE1_dir = 'bowtie-1.2.1.1'
#BOWTIE1_zip = 'bowtie-1.2.1-linux-x86_64.zip'
#BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
#               '1.2.1/' + BOWTIE1_zip
#BOWTIE1_dir = 'bowtie-1.2.1'
BOWTIE1_zip = 'bowtie-1.1.2-linux-x86_64.zip'
BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
               '1.1.2/' + BOWTIE1_zip
BOWTIE1_dir = 'bowtie-1.1.2'

BOWTIE1_target = [os.path.join(tools_dir, BOWTIE1_zip), 
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie'),
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie-inspect'),
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie-build')]
BOWTIE1 = env.Command(BOWTIE1_target, 
                      [], 
                      ['wget -O ${TARGETS[0]} ' + BOWTIE1_link,
                       'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]}']
                     )
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE1[1], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE1[2], SymLink)
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), BOWTIE1[3], SymLink)

# TOPHAT2
tophat2_dir = 'tophat-2.1.0.Linux_x86_64' #'tophat-2.1.1.Linux_x86_64'
tophat2_tar = 'tophat-2.1.0.Linux_x86_64.tar.gz' #'tophat-2.1.1.Linux_x86_64.tar.gz' 
tophat2_url = 'https://ccb.jhu.edu/software/tophat/downloads/' + tophat2_tar

tophat2_target = [os.path.join(tools_dir, tophat2_tar), 
                  os.path.join(tools_dir, tophat2_dir, 'tophat2')]

tophat2 = env.Command(tophat2_target, 
                      [], 
                      ['wget -O ${TARGETS[0]} ' + tophat2_url, 
                       'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}']
                     )
env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), tophat2[1], SymLink)

# STRINGTIE
stringtie_dir =	'stringtie-2.0.3.Linux_x86_64'
stringtie_tar = 'stringtie-2.0.3.Linux_x86_64.tar.gz'
stringtie_url = 'http://ccb.jhu.edu/software/stringtie/dl/' + stringtie_tar
stringtie_target = [os.path.join(tools_dir, stringtie_tar),
		    os.path.join(tools_dir, stringtie_dir, 'stringtie')]
stringtie = env.Command(stringtie_target, 
			[], 
			['wget -O ${TARGETS[0]} ' + stringtie_url,
			 'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}']
			)

stringtie_link = env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), stringtie[1], SymLink)

#use cuffcompare instead!
## GFFCOMPARE
#gffcompare_dir =	'gffcompare-0.10.2.Linux_x86_64'
#gffcompare_tar = 'gffcompare-0.10.2.Linux_x86_64.tar.gz'
#gffcompare_url = 'http://ccb.jhu.edu/software/stringtie/dl/' + gffcompare_tar
#gffcompare_target = [os.path.join(tools_dir, gffcompare_tar),
#		    os.path.join(tools_dir, gffcompare_dir, 'gffcompare')]
#gffcompare = env.Command(gffcompare_target,
#                         [],
#                         ['wget -O ${TARGETS[0]} ' + gffcompare_url,
#                         'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}']
#                         )
#
#gffcompare_link = env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), gffcompare[1], SymLink)

# CIRCRNA_FINDER
cfinder_dir = os.path.join(tools_dir, 'circRNA_finder-1.1')
cfinder_tar = 'v1.1.tar.gz'
cfinder_url = 'https://github.com/orzechoj/circRNA_finder/archive/' + cfinder_tar
cfinder_target = [os.path.join(tools_dir, cfinder_tar), 
                  [os.path.join(cfinder_dir, f) for f in
                                                ['postProcessStarAlignment.pl',
                                                  'filterCirc.awk', 
                                                  'filterSpliceSiteCircles.pl', 
                                                  'nrForwardSplicedReads.pl', 
                                                  'starCirclesToBed.pl']]
                 ]
cfinder = env.Command(cfinder_target,
                      [],
                      ['wget -O ${TARGETS[0]} ' + cfinder_url,
                       'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
for t in cfinder:
    env.Command(os.path.join(ccp_bin_dir, "${SOURCE.file}"), t, SymLink)

# SAMTOOLS <v1.0 is required by circrna_finder
samtools0_dir = os.path.join(tools_dir, 'samtools-0.1.20')
samtools0_tar = '0.1.20.tar.gz'
samtools0_url = 'https://github.com/samtools/samtools/archive/' + samtools0_tar
samtools0_target = [os.path.join(tools_dir, samtools0_tar), 
                    os.path.join(samtools0_dir, 'samtools')]
samtools0 = env.Command(samtools0_target, 
                        [], 
                        ['wget -O $TARGET  ' + samtools0_url,
                        'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                        ' && '.join(['cd ' + samtools0_dir, 
                                     'make']) + 
                        '; cd ' + Dir('.').abspath])
env.Command(os.path.join(ccp_bin_dir, 'samtools_v0', "${SOURCE.file}"), samtools0[1], SymLink)


