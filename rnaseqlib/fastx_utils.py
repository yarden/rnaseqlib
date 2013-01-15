##
## FASTX interface utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils


def fastx_collapse_fastq(fastq_filename, output_dir):
    """
    FASTX collapse FASTQ. Return 
    """
    fastx_collapser = utils.which("fastx_collapser")
    if fastx_collapser is None:
        print "Error: Could not find fastx_collapser."
        return None
    if not os.path.isfile(fastq_filename):
        print "Error: Could not find input fastq %s" \
            %(fastq_filename)
        return None
    output_basename = \
        utils.trim_fastq_ext(os.path.basename(sample_reads_filename))
    collapsed_seq_filename = os.path.join(output_dir,
                                          "%s.collapsed.fastq.gz")
    cat_fastq_cmd = "cat"
    # Handle gzipped input since fastx_collapser does not accept
    # gzipped FASTQ files
    if fastq_filename.endswith(".gz"):
        cut_fastq_cmd = "zcat"
    cat_fastq_cmd += " %s" %(fastq_filename)
    fastx_collapser_cmd = "%s | %s | gzip -c - > %s" \
        %(cat_fastq_cmd,
          fastx_collapser,
          collapsed_seq_filename)
    ret_val = os.system(fastx_collapser_cmd)
    if ret_val != 0:
        print "Error: fastx_collapser command failed."
        return None
    return collapsed_seq_filename
                                          
    
