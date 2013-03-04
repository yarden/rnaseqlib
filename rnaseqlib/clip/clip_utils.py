##
## CLIP utilities
##
import os
import sys
import time
import subprocess

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils
import rnaseqlib.mapping.bedtools_utils as bedtools_utils


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


def filter_clusters(logger, clusters_bed_fname, output_dir,
                    num_reads=8,
                    depth=0,
                    min_size=12,
                    max_size=100):
    """
    Filter clusters by number of reads in them and/or depth.
    """
    if not clusters_bed_fname.endswith(".bed"):
        logger.critical("Error: clusters filename %s must end in .bed" \
                        %(clusters_bed_fname))
        sys.exit(1)
    bed_basename = \
        os.path.basename(clusters_bed_fname).rsplit(".bed", 1)[0]
    filtered_clusters_fname = \
        os.path.join(output_dir, "%s.num_reads_%d.depth_%d.bed" \
                     %(bed_basename,
                       num_reads,
                       depth))
    num_passing_filter = 0
    with open(filtered_clusters_fname, "w") as clusters_out:
        with open(clusters_bed_fname, "r") as clusters_in:
            for line in clusters_in:
                fields = line.strip().split("\t")
                cluster_len = int(fields[2]) - int(fields[1]) + 1
                cluster_reads = len(fields[4].split(","))
                # Check that it meets the number of reads filter
                if cluster_reads < num_reads:
                    continue
                # Check that it meets the depth filter
                cluster_depth = cluster_reads / float(cluster_len)
                if cluster_depth < depth:
                    continue
                # Check that it meets the cluster size filter
                if cluster_len < min_size:
                    continue
                if cluster_len > max_size:
                    continue
                clusters_out.write(line)
                num_passing_filter += 1
    if num_passing_filter == 0:
        logger.critical("0 clusters pass filter.")
        sys.exit(1)
    return filtered_clusters_fname
    

def output_clip_clusters(logger, bam_filename, output_filename,
                         cluster_dist=0,
                         skip_junctions=True):
    """
    In contrast to merge, cluster does not flatten the cluster of
    intervals into a new meta-interval; instead, it assigns an unique
    cluster ID to each record in each cluster.

    Return the output filename or None if it went wrong.
    """
    if cluster_dist != 0:
        raise Exception, "Cluster distance must be 0 for now."
    logger.info("Outputting CLIP clusters...")
    logger.info("  - BAM input: %s" %(bam_filename))
    logger.info("  - Output file: %s" %(output_filename))
    logger.info("  - Cluster dist: %d" %(cluster_dist))
    logger.info("  - Skip junctions: %s" %(str(skip_junctions)))
    # Find clusters in reads:
    #  (1) Convert BAM file to BED on the fly and sort it
    #  (2) Run clusterBed
    #  (3) Merge the cluster
    # Convert BAM -> BED, sort BED
    bamToBed_cmd = "bamToBed -i %s -cigar | sortBed -i - " \
        %(bam_filename)
    # If asked to skip junctions, remove them from BED
    if skip_junctions:
        bamToBed_cmd += " | awk \'$7 !~ \"N\" { print $0 }\'"
    # Cluster the BED
    clusterBed_cmd = "%s | clusterBed -i - " %(bamToBed_cmd)
    # Parse the resulting clusters, skipping junction reads if asked
    cluster_proc = subprocess.Popen(clusterBed_cmd, shell=True,
                                    stdin=sys.stdin,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    # Open stream to mergeBed command and write clusters to it
    mergeBed_cmd = "mergeBed -i - -nms -scores collapse"
    output_file = open(output_filename, "w")
    merge_proc = subprocess.Popen(mergeBed_cmd, shell=True,
                                  stdin=subprocess.PIPE,
                                  stdout=output_file)
    # Process each cluster line and feed it to mergeBed
    for line in iter(cluster_proc.stdout.readline, ""):
        fields = line.strip().split("\t")
        cigar = fields[6]
        if skip_junctions and ("N" in cigar):
            # Skip junctions if asked
            continue
        # Make the score be the cluster number
        fields[3], fields[4] = \
            fields[-1], fields[3]
        processed_line = "%s\n" %("\t".join(fields))
        merge_proc.stdin.write(processed_line)
    merge_proc.communicate()
    output_file.close()
    # Post process the merged clusters
    # Rename the current file
    orig_merged_fname = "%s.orig" %(output_filename)
    os.rename(output_filename, orig_merged_fname)
    # Get rid of semi-colon separated cluster names
    with open(orig_merged_fname, "r") as orig_merged:
        with open(output_filename, "w") as output_file:
            for line in orig_merged:
                merged_fields = line.strip().split("\t")
                # Remove semicolon separated read names
                merged_fields[3] = merged_fields[3].split(";")[0]
                merged_line = "%s\n" %("\t".join(merged_fields))
                output_file.write(merged_line)
    # Remove temporary file
    os.remove(orig_merged_fname)
    return output_filename


def intersect_clusters_with_gff(logger,
                                clusters_fname,
                                gff_filenames,
                                output_dir):
    """
    Intersect CLIP clusters with a set of GFF filenames (e.g. GFFs containing
    events.
    """
    event_clusters_fnames = []
    logger.info("Intersecting clusters %s with GFFs.." \
                %(clusters_fname))
    for gff_fname in gff_filenames:
        if not os.path.isfile(gff_fname):
            logger.critical("Cannot find events GFF file %s" \
                            %(gff_fname))
            continue
        gff_label = os.path.basename(utils.trim_gff_ext(gff_fname))
        event_clusters_fname = \
            os.path.join(output_dir,
                         "%s.clusters.bed" %(gff_label))
        logger.info("  - Processing GFF %s" %(gff_label))
        # Intersect the clusters for this sample with each
        # GFF file
        intersectBed_cmd = \
            "%s -a %s -b %s -loj -f 1 > %s" \
            %(bedtools_utils.intersectBed_path,
              clusters_fname,
              gff_fname,
              event_clusters_fname)
        if os.path.isfile(event_clusters_fname):
            logger.info("Found %s, skipping.." %(event_clusters_fname))
            event_clusters_fnames.append(event_clusters_fname)
        logger.info("  - Executing: %s" %(intersectBed_cmd))
        ret_val = os.system(intersectBed_cmd)
        if ret_val != 0:
            logger.critical("Cannot intersect clusters with %s" \
                            %(gff_fname))
        else:
            event_clusters_fnames.append(event_clusters_fname)
    return event_clusters_fnames

