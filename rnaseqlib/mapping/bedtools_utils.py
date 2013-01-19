##
## Basic wrappers for bedtools
##
import os
import sys
import time
import subprocess

import rnaseqlib
import rnaseqlib.utils as utils


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
                               frac_overlap=1,
                               bed_opts=""):
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
        print "BAM file %s does not exist." %(bam_filename)
        return None
    if not os.path.isfile(bed_filename):
        print "BED file %s does not exist." %(bed_filename)
        return None
    bedtools_cmd = \
      "bamToBed -i %(bam_filename)s | " \
      "intersectBed -a %(bed_filename)s -b - -sorted -f %(frac_overlap)s -wo | " \
      "groupBy -g 1-4 -c 9 -o collapse" %(args)
    print "Executing: %s" %(bedtools_cmd)
    ret_val = os.system(bedtools_cmd)
    if ret_val != 0:
        print "bedtools call failed."
        return None
    return output_filename
    

def count_reads_matching_intervals(bam_filename,
                                   intervals_filename,
                                   output_filename):
    """
    Count number of reads matching the given intervals
    filename, which could be a BED or a GFF.

    Output the resulting file (a BED) to output filename.
    """
    if not os.path.isfile(intervals_filename):
        print "WARNING: Cannot find exons file %s" %(output_filename)
        return None
    # intersect BAM with merged exons, outputting one line per
    # read that overlaps any region in the merged exons BED
    output_filename = get_reads_matching_regions(bam_filename,
                                                 intervals_filename,
                                                 output_filename)
    if output_filename is None:
        return 0
    if not os.path.isfile(output_filename):
        return None
    num_reads = utils.count_lines(output_filename)
    return num_reads


def sort_bed(input_filename, output_filename):
    """
    Sort a BED file. Input filename can be a GFF or
    a BED file.
    """
    print "Sorting BED: %s" %(input_filename)
    if os.path.isfile(output_filename):
        print "%s exists. Skipping..." %(output_filename)
    sort_cmd = "sortBed -i %s > "


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
          "sortBed -i %s | mergeBed -i stdin -nms -s | sortBed -i - > %s" \
            %(input_filename,
              output_filename)
        print "Executing: %s" %(merge_cmd)
        os.system(merge_cmd)
    else:
        raise Exception, "Not implemented."
    return output_filename

    
def make_bed_line(chrom, start, end,
                  name, score, strand,
                  delimiter="\t"):
    """
    Return a BED formatted line.
    """
    bed_fields = map(str, [chrom, start, end, name, score, strand])
    bed_line = delimiter.join(bed_fields)
    return bed_line


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
    sort_bed_cmd = "sortBed -i %s > %s" %(bed_filename,
                                          tmp_filename)
    ret_val = os.system(sort_bed_cmd)
    if ret_val != 0:
        print "Sorting failed."
        return None
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
