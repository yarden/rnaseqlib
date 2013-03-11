##
## Basic wrappers for bedtools
##
import os
import sys
import time
import shlex
import subprocess

import numpy as np

from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.coords_utils as coords_utils

import pybedtools

##
## Paths to bedtools programs
##
intersectBed_path = utils.which("intersectBed")
mergeBed_path = utils.which("mergeBed")
tagBam_path = utils.which("tagBam")
coverageBed_path = utils.which("coverageBed")
fastaFromBed_path = utils.which("fastaFromBed")


def sample_intervals(start, end, interval_size,
                     num_intervals=1,
                     exclude_interval=None,
                     max_num_tries=100):
    """
    Sample interval of size 'interval_size' on [start, end)
    using a naive sampling approach.
    """
    intervals = []
    last_start = end - interval_size
    possible_interval_size = last_start - start
    if possible_interval_size < interval_size:
        # Possible interval to sample from too small
        # to get something of size 'interval_size'
        return None
    for n in range(num_intervals):
        interval_accepted = False
        # Sample intervals until one is accepted
        curr_interval = None
        num_tries = 0
        while not interval_accepted:
            # If the number of tries exceeds a limit, raise
            # exception
            if num_tries >= max_num_tries:
                # print "Cannot sample interval on %d, %d " \
                #       "while excluding %d, %d" \
                #                  %(start,
                #                    end,
                #                    exclude_interval[0],
                #                    exclude_interval[1])
                return None
            interval_start = np.random.randint(low=start,
                                               high=last_start)
            interval_end = interval_start + interval_size
            curr_interval = (interval_start, interval_end)
            if exclude_interval is not None:
                # Check that the exclusion interval is not present
                # if we were given one.
                interval_overlap = coords_utils.overlap(curr_interval,
                                                        exclude_interval)
                if interval_overlap <= 0:
                    # No overlap with exclusion interval, so accept it
                    interval_accepted = True
            else:
                interval_accepted = True
            num_tries += 1
        # Accumulate accepted intervals
        intervals.append(curr_interval)
    return intervals
    

def fastaFromBed(logger,
                 genome_seq_fname,
                 bed_fname,
                 output_fname,
                 s=False):
    logger.info("Retrieving FASTA from BED...")
    try:
        bed_file = pybedtools.BedTool(bed_fname)
        result = bed_file.sequence(fi=genome_seq_fname,
                                   bed=bed_fname,
                                   fo=output_fname,
                                   s=s)
        if not os.path.isfile(result.seqfn):
            # If the output FASTA filename is not present,
            # assume the call failed
            logger.critical("fastaFromBed call failed.")
            sys.exit(1)
    except pybedtools.helpers.BEDToolsError as bed_error:
        if "WARNING" in str(bed_error):
            logger.info("BEDTools produced warning: continuing anyway.")
        else:
            logger.critical("BEDTools failed.")
    return output_fname


def intersect_bed_with_gene_regions(bed_filename,
                                    ucsc_tables_dir,
                                    output_bed_fname):
    """
    Take a BED file and intersect it with a set of gene
    regions from the GFF gene annotation file.

    Output a BED file with annotations.
    """
    # Read GFF annotation fname: on the fly
    # concatenate introns, exons, cds_only exons, and
    # the three_prime_utr/five_prime_utr annotations of
    # of the tables into one BED file.
    # Intersect bed with pybedtools
    # ...
    # groupby combine the columns to have the gene name(s)
    # and the region type(s) that it overlaps with
    # ...
    # Output result as BED
    # ...
    return output_bed_fname


def intersect_bam_with_bed(bam_filename,
                           bed_filename,
                           output_filename,
                           frac_overlap=1,
                           unique=False,
                           bed_opts=""):
    """
    Intersect a BAM filename with BED.

    If 'unique' is True, then pass -u option to intersectBed,
    which states that each entry in the BAM is written once if it
    overlaps any BED
    """
    print "Intersecting %s with %s" %(bam_filename,
                                      bed_filename)
    if os.path.isfile(output_filename):
        print "%s exists. Skipping..." %(output_filename)
        return output_filename
    # Intersect given BAM with given BED and write result
    # as a BED file
    intersect_cmd = "intersectBed -abam %s -b %s -f %s -bed -wb " \
        %(bam_filename,
          bed_filename,
          str(frac_overlap))
    if unique:
        print "Looking only at unique entries"
        intersect_cmd += " -u"
    # Other bed options
    intersect_cmd += " %s " %(bed_opts)
    intersect_cmd += "> %s" %(output_filename)
    print "Executing: %s" %(intersect_cmd)
    os.system(intersect_cmd)
    time.sleep(5)
    return output_filename


def get_reads_matching_regions(bam_filename,
                               bed_filename,
                               output_filename,
                               logger,                               
                               frac_overlap=1):
    """
    Get comma-separated list of reads from BAM that match each
    region in the given BED file.

    Proceeds by:

    - First converting BAM to BED
    - Intersecting BAM-as-BED with BED regions (outputting
      result as sorted BED)
    - Grouping the regions by region using groupBy. Uses columns 1-4
      of the BED as entries and groups by the 9th column which
      contains the read IDs from the BAM.

    Note: expects that the input BAM file to be sorted.
    """
    print "Getting reads matching regions."
    args = {"bam_filename": bam_filename,
            "bed_filename": bed_filename,
            "output_filename": output_filename,
            "frac_overlap": str(frac_overlap)}
    if not os.path.isfile(bam_filename):
        logger.critical("BAM file %s does not exist." %(bam_filename))
        return None
    if not os.path.isfile(bed_filename):
        logger.critical("BED file %s does not exist." %(bed_filename))
        return None
    bedtools_cmd = \
      "bamToBed -i %(bam_filename)s -split | sortBed -i - | " \
      "intersectBed -a %(bed_filename)s -b - -sorted -f %(frac_overlap)s -wo | " \
      "groupBy -g 1-4 -c 9 -o collapse > %(output_filename)s" %(args)
    logger.info("Executing: %s" %(bedtools_cmd))
    ret_val = os.system(bedtools_cmd)
    if ret_val != 0:
        logger.critical("bedtools call failed.")
        return None
    return output_filename


def multi_tagBam(bam_filename, intervals_files, intervals_labels,
                 output_filename, logger):
    """
    Call tagBam mapping BAM file to multiple bed/gff files.

    Takes:

    - bam_filename: The BAM file path
    - intervals_files: List of string paths for the interval files
      (BED or GFF file format)
    - intervals_labels: Labels for each of the files in 'intervals_files'
    - output_filename: BAM filename to use as output
    - logger: a logger to log messages to
    """
    num_interval_files = len(intervals_files)
    logger.info("Running tagBam against %d interval files.." \
                %(num_interval_files))
    tagBam = utils.which("tagBam")
    if tagBam is None:
        logger.critical("tagBam not found.")
        return None
    if os.path.isfile(output_filename):
        logger.info("Found %s, skipping." %(output_filename))
        return output_filename
    t1 = time.time()
    args = {"tagBam": tagBam,
            "bam_filename": bam_filename,
            "intervals_files": " ".join(intervals_files),
            "intervals_labels": " ".join(intervals_labels),
            "output_filename": output_filename}
    tagBam_cmd = \
      "%(tagBam)s -i %(bam_filename)s -files %(intervals_files)s " \
      "-labels %(intervals_labels)s -intervals -f 1 > %(output_filename)s" \
      %(args)
    logger.info("Executing: %s" %(tagBam_cmd))
    ret_val = os.system(tagBam_cmd)
    t2 = time.time()
    logger.info("tagBam took %.2f minutes." %((t2 - t1)/60.))
    if ret_val != 0:
        logger.critical("tagBam command failed.")
        return None
    return output_filename


def coverageBed(bam_filename,
                intervals_filename,
                output_filename,
                logger):
    """
    Run coverageBed, mapping BAM against the given intervals
    filename (a GFF or BED). 
    """
    logger.info("Running coverageBed..")
    logger.info("  Mapping %s against %s" %(bam_filename,
                                            intervals_filename))
    logger.info("  Output file: %s" %(output_filename))
    t1 = time.time()
    if not os.path.isfile(intervals_filename):
        logger.critical("Cannot find intervals file %s" %(output_filename))
    if not os.path.isfile(bam_filename):
        logger.critical("Cannot find BAM file %s" %(bam_filename))
    if os.path.isfile(output_filename):
        logger.info("Found %s, skipping.." %(output_filename))
    args = {"bam_filename": bam_filename,
            "intervals_filename": intervals_filename,
            "output_filename": output_filename}
    coverageBed_cmd = \
       "coverageBed -abam %(bam_filename)s -b %(intervals_filename)s -split " \
       "> %(output_filename)s" %(args)
    logger.info("Executing: %s" %(coverageBed_cmd))
    ret_val = os.system(coverageBed_cmd)
    if ret_val != 0:
        logger.critical("coverageBed command failed.")
        return None
    t2 = time.time()
    logger.info("coverageBed took %.2f minutes." %((t2 - t1)/60.))
    return output_filename
    

def count_reads_matching_intervals(bam_filename,
                                   intervals_filename,
                                   output_filename,
                                   logger):
    """
    Count number of reads matching the given intervals
    filename, which could be a BED or a GFF.

    Output the resulting file (a BED) to output filename.
    """
    if not os.path.isfile(intervals_filename):
        logger.critical("Cannot find exons file %s" %(output_filename))
        return None
    # intersect BAM with merged exons, outputting one line per
    # read that overlaps any region in the merged exons BED
    output_filename = get_reads_matching_regions(bam_filename,
                                                 intervals_filename,
                                                 output_filename,
                                                 logger)
    if output_filename is None:
        logger.warning("Cannot access matching regions file %s" \
                       %(output_filename))
        return 0
    if not os.path.isfile(output_filename):
        logger.warning("Matching regions file %s does not exist" \
                       %(output_filename))
        return None
    num_reads = utils.count_lines(output_filename)
    return num_reads


def merge_bed(input_filename, output_filename,
              sort_input=True):
    """
    Merge a BED file. Input filename can be a GFF or a
    BED file.
    """
    merge_cmd = ""
    if os.path.isfile(output_filename):
        print "%s exists. Skipping..." %(output_filename)
        return output_filename
    if sort_input:
        merge_cmd = \
          "sortBed -i %s | mergeBed -i stdin -nms -s -scores sum | sortBed -i - " \
            %(input_filename)
        print "Executing: %s" %(merge_cmd)
        output_file = open(output_filename, "w")
        proc = subprocess.Popen(merge_cmd,
                                stdout=subprocess.PIPE,
                                shell=True)
        sub_bed_delimiter(proc.stdout, output_file,
                          from_delim=";",
                          to_delim=",")
        proc.communicate()
    else:
        raise Exception, "Not implemented."
    return output_filename


def sub_bed_delimiter(input_file, output_file,
                      from_delim=";",
                      to_delim=","):
    """
    Substitute BED delimiter in names field.
    """
    for line in input_file:
        bed_fields = line.strip().split("\t")
        bed_fields[3] = bed_fields[3].replace(from_delim,
                                              to_delim)
        bed_output_line = "%s\n" %("\t".join(bed_fields))
        output_file.write(bed_output_line)
    output_file.close()
        
    
def make_bed_line(chrom, start, end,
                  name, score, strand,
                  delimiter="\t"):
    """
    Return a BED formatted line.
    """
    bed_fields = map(str, [chrom, start, end,
                           name, score, strand])
    bed_line = delimiter.join(bed_fields)
    return bed_line


def sort_bed(input_filename, output_filename,
             gff_to_bed=False,
             attribute_to_use=None):
    """
    Sort BED file with sortBed.

    - If 'gff_to_bed' is given, assumes the output of sortBed
      is in GFF format, and converts it on the fly to BED
    """
    print "Sorting BED %s" %(input_filename)
    if not os.path.isfile(input_filename):
        print "Error: %s does not exist." %(input_filename)
        return None
    if os.path.isfile(output_filename):
        print "%s exists, skipping sortBed" %(output_filename)
        return output_filename
    input_file = open(input_filename, "r")
    output_file = open(output_filename, "w")
    print "  - Output file: %s" %(output_filename)
    if gff_to_bed:
        proc = subprocess.Popen(["sortBed", "-i"],
                                stdin=input_file,
                                stdout=subprocess.PIPE)
        # Iterate through stdout
        convert_gff_to_bed(proc.stdout, output_file,
                           attribute_to_use=attribute_to_use)
    else:
        proc = subprocess.Popen(["sortBed", "-i"],
                                stdin=input_file,
                                stdout=output_file)
    proc.communicate()
    output_file.close()
    print "  - Sorting completed."
    return output_filename
    

def sort_bedfile_inplace(bed_filename):
    """
    Sort the given BED filename and substitute
    the sorted version in its place.
    """
    if not os.path.isfile(bed_filename):
        print "Error: BED file %s does not exist, cannot sort." \
            %(bed_Filename)
        return None
    tmp_filename = "%s.tmp_sorted" %(bed_filename)
    result = sort_bed(bed_filename, tmp_filename)
    if result is None:
        return None
    # Rename tmp file to be original filename
    os.rename(tmp_filename, bed_filename)
    return bed_filename


def output_intervals_as_bed(out_file, chrom, interval_coords, strand,
                            name="exon",
                            score="1"):
    """
    Output intervals as BED.

    Takes:

    - out_file: file handle to the BED
    - chrom: chromosome
    - interval_coords: list of (start, end) coords for parts
    - strand: the strand

    NOTE: Assumes interval_coords are in 0-based format already!
    """
    for part in interval_coords:
        start, end = part
        bed_line = make_bed_line(chrom, start, end,
                                 name, score, strand)
        out_file.write("%s\n" %(bed_line))


def convert_gff_to_bed(input_stream, output_stream,
                       attribute_to_use=None):
    """
    Convert GFF lines from input_stream to
    BED output_stream.
    """
    for line in input_stream:
        if line.startswith("#"):
            # Skip GFF comments
            continue
        gff_fields = line.strip().split("\t")
        chrom = gff_fields[0]
        rec_source = gff_fields[1]
        rec_type = gff_fields[2]
        start = int(gff_fields[3])
        end = int(gff_fields[4])
        strand = gff_fields[6]
        if start > end:
            # Flip coordinates if start > end
            start, end = end, start
        # Convert start to be 0-based
        start = start - 1
        # Convert coordinates back to strings
        start, end = str(start), str(end)
        rec_attributes = gff_fields[-1]
        name = rec_type
        score = "1"
        if rec_attributes != ".":
            # If there's an ID= attribute present,
            # use its value as the BED entry's name
            attributes = utils.parse_attributes(rec_attributes)
            if attribute_to_use is not None:
                if attribute_to_use in attributes:
                    # Use gene_id as ID if found
                    name = attributes[attribute_to_use]
                else:
                    print "WARNING: %s attribute not found" \
                        %(attribute_to_use)
                    name = attributes["ID"]
            elif "ID" in attributes:
                name = attributes["ID"]
        bed_line = make_bed_line(chrom, start, end,
                                 name, score, strand)
        output_stream.write("%s\n" %(bed_line))


def get_bed_lens(bed_filename):
    """
    Get the lengths of all the intervals of the
    BED filename. Returns dictionary mapping coordinate
    string to length.
    """
    region_lens = {}
    total_len = 0
    with open(bed_filename) as bed_file:
        for line in bed_file:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start_coord, end_coord = \
                int(fields[1]), int(fields[2])
            coordinate_id = "%s:%s-%s" \
                %(chrom, start_coord, end_coord)
            if start_coord > end_coord:
                start_coord, end_coord = \
                    end_coord, start_coord
            region_len = end_coord - start_coord + 1
            region_lens[coordinate_id] = region_len
            # Keep track of total length of all regions
            total_len += region_len
    return region_lens, total_len
    

def parse_tagBam_region(regions_field, trans_prefix="ENS"):
    """
    Parse tagBam region. Returns read coordinates
    detected and a set of regions detected.

    Assumes transcripts being with trans_prefix.
    """
    read_categories = regions_field.split(";")
    mapped_parts = utils.flatten([read_cat.split(",") \
                                  for read_cat in read_categories])
    region_coordinates = []
    for part in utils.chunk_list(mapped_parts, ":"):
        region_type, region_coord = part[0].split(":", 1)[0:2]
        curr_transcripts = set([t for t in part if t.startswith(trans_prefix)])
        region_coordinates.append((region_coord, region_type, curr_transcripts))
    # Get coordinates that read fall in
    # Get region types detected
    regions_detected = \
        map(lambda x: x.split(":")[0],
            filter(lambda x: ":" in x, read_categories))
    return region_coordinates, regions_detected

