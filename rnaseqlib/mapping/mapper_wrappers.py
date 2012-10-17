##
## Mapping utilities
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.settings

def get_tophat_mapping_cmd(tophat_path,
                           sample,
                           output_dir,
                           settings_info,
                           num_processors=4):
    """
    Get tophat args for mapping for a sample.
    """
    tophat_path = settings_info["mapping"]["tophat_path"]
    index_filename = settings_info["mapping"]["tophat_index"]
    tophat_options = settings_info["mapping"]["tophat_options"]
    tophat_gtf = None
    if "tophat_gtf" in settings_info["mapping"]:
        tophat_gtf = settings_info["mapping"]["tophat_gtf"]
    mapper_cmd = "%s -p %d %s" \
                 %(tophat_path,
                   num_processors,
                   tophat_options)
    if tophat_gtf is not None:
        mapper_cmd += " --GTF %s" %(tophat_gtf)
    # If paired-end, get a pair of files for the sample
    if sample.paired:
        input_files = " ".join([sample.samples[0].reads_filename,
                                sample.samples[1].reads_filename])
        mapper_cmd += " --output-dir %s %s %s" %(output_dir,
                                                 index_filename,
                                                 input_files)
    else:
        input_files = sample.samples[0].reads_filename
    tophat_outfilename = os.path.join(output_dir,
                                      "accepted_hits.bam")
    return mapper_cmd, tophat_outfilename


def get_bowtie_mapping_cmd(bowtie_path,
                           input_filename,
                           genome_index_filename,
                           output_filename,
                           bowtie_options=""):
    """
    Get bowtie args for mapping.
    """
    input_compressed = False
    if input_filename.endswith(".gz"):
        input_compressed = True
    # Make output always gzipped
    if ("--sam" not in bowtie_options):
        # Always output sam
        bowtie_options += " --sam"
    output_filename = "%s.bam" %(output_filename)
    args = {"bowtie_path": bowtie_path,
            "input_filename": input_filename,
            "genome_index_filename": genome_index_filename,
            "output_filename": output_filename,
            "bowtie_options": bowtie_options}
    if input_compressed:
        # Assume the input is compressed. Pass it through
        # bowtie via gzip
        # The "-c" argument to -gunzip is for printing to stdout
        mapper_cmd = "gunzip -c %(input_filename)s | %(bowtie_path)s %(bowtie_options)s " \
                     "%(genome_index_filename)s - | samtools view -bS - -o %(output_filename)s" % args
    else:
        mapper_cmd = "%s %s " %(bowtie_path, bowtie_options)
        mapper_cmd += " %s %s - | samtools view -Sbh - > %s" \
            %(genome_index_filename,
              input_filename,
              output_filename)
    return mapper_cmd, output_filename


def map_bam2gff(bedtools_dir,
                bam_filename,
                gff_filename,
                output_dir,
                settings,
                queue_name=None):
    """
    Map BAM file against intervals in GFF, return results as BAM.

    Uses tagBam utility from bedtools.
    """
    bedtools_dir = settings["paths"]["bedtools_path"]
    # Path to tagBam utility
    tagBam = os.path.abspath(os.path.join(bedtools_dir, "tagBam"))
    gff_basename = os.path.basename(gff_filename)
    bam_basename = os.path.basename(bam_filename)
    output_dir = os.path.join(output_dir, "bam2gff_%s" \
                              %(gff_basename))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    output_filename = os.path.join(output_dir, bam_basename)

    # "-intervals" option embeds the original GFF coordinates
    # in the output BAM file. 
    tagBam_cmd = "%s -i %s -files %s -labels gff -intervals -f 1 > %s" \
                 %(tagBam, bam_filename, gff_filename,
                   output_filename)
    print "Running: %s" %(tagBam_cmd)
    print "  - Output: %s" %(output_filename)

    job_id = "bam2gff_%s_%s" %(bam_basename,
                               gff_basename)
    cluster.run_on_cluster(tagBam_cmd, job_id, output_dir,
                           queue_name=queue_name)

