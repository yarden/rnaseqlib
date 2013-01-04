##
## CLIP utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

def trim_clip_adaptors(fastq_filename,
                       output_dir):
    """
    Trim CLIP adaptors using 'cutadapt'.
    """
    cutadapt_path = utils.which("cutadapt"):
    if cutadapt_path is None:
        print "Could not find \'cutadapt\' on the path."
        print "Please install \'cutadapt\' or make the installed " \
              "version available on path."
    

