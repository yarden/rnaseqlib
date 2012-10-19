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
#      license = 'MIT License',
      author = 'Yarden Katz,Eric T. Wang,Noah Spies',
#      data_files = [],
      packages = ['rnaseqlib', 'rnaseqlib.init',
                  'rnaseqlib.genes', 'rnaseqlib.qc',
                  'rnaseqlib.mapping', 'rnaseqlib.cluster_utils',
                  'rnaseqlib.plotting', 'rnaseqlib.events',
                  'rnaseqlib.clip', 'rnaseqlib.ribo'],
      # distutils always uses forward slashes      
      scripts = ['scripts/rna_pipeline.py',
                 'scripts/gtf2gff3.pl'],
      # Required modules
      install_requires = [
#          "matplotlib >= 1.1.0",
          "matplotlib",
          "numpy >= 1.5.0",
          "scipy >= 0.9.0",
          "pysam >= 0.6.0",
          "misopy >= 0.4.7",
          "pandas >= 0.8.1"
          ],
      platforms = 'ALL',
      keywords = ['bioinformatics', 'sequence analysis',
                  'alternative splicing', 'RNA-Seq'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
      )

