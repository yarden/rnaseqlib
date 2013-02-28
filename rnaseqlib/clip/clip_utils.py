##
## CLIP utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils


def trim_clip_adaptors(fastq_filename,
                       adaptors_filename,
                       output_dir,
                       logger,
                       min_read_len=5):
    """
    Trim CLIP adaptors using 'cutadapt'.
    """
    logger.info("Trimming CLIP adaptors from: %s" %(fastq_filename))
    cutadapt_path = utils.which("cutadapt")
    if cutadapt_path is None:
        logger.critical("Could not find \'cutadapt\' on the path. " \
                        "Please install \'cutadapt\' or make the installed " \
                        "version available on path.")
    output_basename = \
        utils.trim_fastq_ext(os.path.basename(fastq_filename))
    output_filename = os.path.join(output_dir,
                                   "%s_trimmed.fastq.gz" \
                                   %(output_basename))
    if os.path.isfile(output_filename):
        logger.info("SKIPPING: %s already exists!" \
                    %(output_filename))
        return output_filename
    logger.info("  - Outputting trimmed sequences to: %s" \
                %(output_filename))
    # Load adaptors to pass to 'cutadapt'
    if not os.path.isfile(adaptors_filename):
        logger.critical("Could not find adaptors file %s" \
                        %(adaptors_filename))
        sys.exit(1)
    adaptors_in = open(adaptors_filename, "r")
    # Substitute newlines with spaces
    adaptors = adaptors_in.read().strip().replace("\n", " ")
    adaptors_in.close()
    cutadapt_cmd = "%s %s %s -o %s -m %d -q 3 > %s.log" %(cutadapt_path,
                                                          adaptors,
                                                          fastq_filename,
                                                          output_filename,
                                                          min_read_len,
                                                          output_filename)
    logger.info("Executing: %s" %(cutadapt_cmd))
    t1 = time.time()
    os.system(cutadapt_cmd)
    t2 = time.time()
    logger.info("Trimming took %.2f mins." %((t2 - t1)/60.))
    return output_filename


def collapse_clip_reads(sample, output_dir, logger):
    """
    Collapse CLIP reads. Uses fastx_collapser.
    """
    logger.info("Collapsing CLIP reads for %s" %(sample.label))
    t1 = time.time()
    collapsed_seq_filename = \
        fastx_utils.fastx_collapse_fastq(sample.rawdata.reads_filename,
                                         output_dir,
                                         logger)
    if collapsed_seq_filename is None:
        logger.critical("Collapsing of CLIP reads failed.")
        sys.exit(1)
    t2 = time.time()
    logger.info("Collapsing took %.2f minutes." %((t2 - t1)/60.))
    return collapsed_seq_filename


def check_clip_utils(logger,
                     required_utils=["cutadapt",
                                     "fastx_collapser"]):
    """
    Check that necessary utilities are available.
    """
    logger.info("Checking that utilities required for CLIP are available..")
    for program in required_utils:
        program_path = utils.which(program)
        if program_path is None:
            logger.critical("Could not access: %s" %(program))
            logger.critical("Make %s avaialble and try again." %(program))
            sys.exit(1)
    logger.info("Found CLIP utilities.")
    

def output_clip_clusters(logger, bam_filename, output_filename,
                         cluster_dist=0):
    """
    In contrast to merge, cluster does not flatten the cluster of
    intervals into a new meta-interval; instead, it assigns an unique
    cluster ID to each record in each cluster.

    Return the output filename or None if it went wrong.
    """
    if cluster_dist != 0:
        raise Exception, "Cluster distance must be 0 for now."
    # Find clusters in reads:
    #  (1) Convert BAM file to BED on the fly and sort it
    #  (2) Run clusterBed
    #  (3) Merge the cluster
    # Convert BAM -> BED, sort BED
    bamToBed_cmd = "bamToBed -i %s -split | sortBed -i - " \
        %(bam_filename)
    # Cluster the BED
    clusterBed_cmd = "clusterBed -i - "
    # Merge the clusters while recording read IDs in each cluster
    # Do not make it stranded
    mergeBed_cmd = "mergeBed -i - -nms"
    logger.info("Outputting CLIP clusters...")
    logger.info("  - BAM input: %s" %(bam_filename))
    logger.info("  - Output file: %s" %(output_filename))
    bedtools_cmd = "%s | %s | %s > %s" %(bamToBed_cmd,
                                         clusterBed_cmd,
                                         mergeBed_cmd,
                                         output_filename)
    ret_val = os.system(bedtools_cmd)
    if ret_val != 0:
        return None
