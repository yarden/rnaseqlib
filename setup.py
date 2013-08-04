from distutils.core import setup, Extension
import glob
import os
import sys

# This forces distutils to place the data files
# in the directory where the Py packages are installed
# (usually 'site-packages'). This is unfortunately
# required since there's no good way to retrieve
# data_files= from setup() in a platform independent
# way.
from distutils.command.install import INSTALL_SCHEMES
for scheme in INSTALL_SCHEMES.values():
        scheme['data'] = scheme['purelib']
        
setup(name = 'rnaseqlib',
      version = '0.1',
      description = 'RNA-Seq analysis pipeline',
      author = 'Yarden Katz',
      author_email = 'yarden@mit.edu',
      data_files = [('examples/rnaseq',
                     ['examples/rnaseq/rnaseq_settings.txt']),
                    ('examples/rnaseq/fastq',
                     ['examples/rnaseq/fastq/sample1_p1.fastq.gz',
                      'examples/rnaseq/fastq/sample1_p2.fastq.gz',
                      'examples/rnaseq/fastq/sample2_p1.fastq.gz',
                      'examples/rnaseq/fastq/sample2_p2.fastq.gz'])],
      packages = ['rnaseqlib',
                  'rnaseqlib.init',
                  'rnaseqlib.genes',
                  'rnaseqlib.rpkm',
                  'rnaseqlib.mapping',
                  'rnaseqlib.cluster_utils',
                  'rnaseqlib.plotting',
                  'rnaseqlib.miso',
                  'rnaseqlib.events',
                  'rnaseqlib.clip',
                  'rnaseqlib.gff',
                  'rnaseqlib.bam',
                  'rnaseqlib.motif',
                  'rnaseqlib.ribo',
                  'rnaseqlib.stats',
                  'rnaseqlib.tests',
                  'rnaseqlib.scripts'],
      # Only non-Python scripts go here
      # distutils always uses forward slashes      
      scripts = ['scripts/gtf2gff3.pl',
                 'scripts/ucsc_table2gff3.pl'],
      # Python scripts are defined as entry points
      entry_points = {
               'console_scripts':
               ['rna_pipeline = rnaseqlib.scripts.rna_pipeline:main',
                # MISO-related scripts
                'misowrap = rnaseqlib.miso.misowrap:main',
                'intersect_events = rnaseqlib.miso.intersect_events:main',
                'miso_helper = rnaseqlib.miso.miso_helper:main',
                # GFF-related scripts
                'gff_annotate_events = rnaseqlib.gff.gff_annotate_events:main',
                'gff_make_annotation = rnaseqlib.gff.gff_make_annotation:main',
                'gff_sanitize = rnaseqlib.gff.gff_sanitize:main',
                'gff_create_db = rnaseqlib.gff.gff_create_db:main',
                'gff_extract_lens = rnaseqlib.gff.gff_extract_lens:main',
                'gff_extract_event_seqs = rnaseqlib.gff.gff_extract_event_seqs:main']
      },
      # Required modules
      install_requires = [
#          "matplotlib >= 1.1.0",
          "matplotlib",
          "numpy >= 1.5.0",
          "pysam >= 0.6.0",
          "misopy >= 0.4.6",
          "pandas >= 0.8.1",
          # pybedtools
          "pybedtools",
          "cutadapt",
          # Khmer: kmer counter 
          #"khmer",
          ],
      platforms = 'ALL',
      keywords = ['bioinformatics', 'sequence analysis',
                  'alternative splicing', 'RNA-Seq'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
      )

