##
## FASTX interface utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

import rnaseqlib.fasta_utils as fasta_utils
import rnaseqlib.fastq_utils as fastq_utils


def get_fastx_entries(fastx_filename):
    """
    Get entries of FASTQ/FASTA file.
    """
    if (fastx_filename.endswith(".fasta") or \
        fastx_filename.endswith(".fasta.gz")):
        # It's a FASTA file
        entries = fasta_utils.read_fasta(fastx_filename)
    elif (fastx_filename.endswith(".fastq") or \
          fastx_filename.endswith(".fastq.gz")):
        # It's a FASTQ file
        entries = fastq_utils.read_fastq(fastx_filename)
    return entries


def fastx_collapse_fastq(fastq_filename, output_dir, logger):
    """
    FASTX collapse FASTQ. Return 
    """
    fastx_collapser = utils.which("fastx_collapser")
    if fastx_collapser is None:
        logger.critical("Could not find fastx_collapser.")
        return None
    if not os.path.isfile(fastq_filename):
        logger.critical("Could not find input fastq %s" \
                        %(fastq_filename))
        return None
    output_basename = \
        utils.trim_fastq_ext(os.path.basename(fastq_filename))
    collapsed_seq_filename = os.path.join(output_dir,
                                          "%s.collapsed.fasta.gz" \
                                          %(output_basename))
    if os.path.isfile(collapsed_seq_filename):
        logger.info("%s exists, skipping collapsing step." \
                    %(collapsed_seq_filename))
        return collapsed_seq_filename
    cat_fastq_cmd = "cat"
    # Handle gzipped input since fastx_collapser does not accept
    # gzipped FASTQ files
    if fastq_filename.endswith(".gz"):
        cat_fastq_cmd = "zcat"
    cat_fastq_cmd += " %s" %(fastq_filename)
    # Use -Q 33 flag to signal Illumina quality scores to
    # FASTX-Toolkit
    fastx_collapser_cmd = "%s | %s -Q 33 | gzip -c - > %s" \
        %(cat_fastq_cmd,
          fastx_collapser,
          collapsed_seq_filename)
    logger.info("Executing: %s" %(fastx_collapser_cmd))
    ret_val = os.system(fastx_collapser_cmd)
    if ret_val != 0:
        logger.critical("Error: fastx_collapser command failed.")
        return None
    return collapsed_seq_filename
                                          
    
