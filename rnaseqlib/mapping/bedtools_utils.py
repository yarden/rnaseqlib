##
## Basic wrappers for bedtools
##
import os
import sys
import time

def intersect_bam_with_bed(bam_filename,
                           bed_filename,
                           output_filename,
                           frac_overlap=1):
    """
    Intersect a BAM filename with BED.
    """
    print "Intersecting %s with %s" %(bam_filename,
                                      bed_filename)
    if os.path.isfile(output_filename):
        print "%s exists. Skipping..." %(output_filename)
        return output_filename
    # Intersect given BAM with given BED and write result
    # as a BED file
    intersect_cmd = "intersectBed -abam %s -b %s -f %s -bed > %s" \
        %(bam_filename,
          bed_filename,
          str(frac_overlap),
          output_filename)
    os.system(intersect_cmd) 
    return output_filename


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
        merge_cmd = "sortBed -i %s | mergeBed -i stdin -nms > %s" \
            %(input_filename,
              output_filename)
        print "Executing: %s" %(merge_cmd)
        os.system(merge_cmd)
    else:
        raise Exception, "Not implemented."
    return output_filename


def output_exons_as_bed(out_file, chrom, exon_coords, strand,
                        name="exon",
                        score="1",
                        delimiter="\t"):
    """
    Output exons as BED.

    Takes:

    - out_file: file handle to the BED
    - chrom: chromosome
    - exon_coords: list of (start, end) coords for exons
    - strand
    """
    for exon in exon_coords:
        start, end = exon
        bed_fields = map(str, [chrom, start, end, name, score])
        bed_line = delimiter.join(bed_fields)
        out_file.write("%s\n" %(bed_line))
