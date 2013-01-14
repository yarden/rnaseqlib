##
## FASTX interface utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils


def fastx_collapse_fastq(fastq):
    """
    FASTX collapse FASTQ. Return 
    """
    fastx_collapser = utils.which("fastx_collapser")
    if fastx_collapser is None:
        print "Error: Could not find fastx_collapser."
        sys.exit(1)
    
