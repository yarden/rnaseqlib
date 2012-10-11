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
#      scripts = [],
#      data_files = [],
      # Required modules
      install_requires = [
#          "matplotlib >= 1.1.0",
          "matplotlib",
          "numpy >= 1.5.0",
          "scipy >= 0.9.0",
          "pysam >= 0.6.0",
          "pandas"
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

