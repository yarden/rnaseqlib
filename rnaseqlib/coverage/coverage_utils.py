##
## Utilities related to coverage analysis with BAMs
##
import os
import sys
import time

import numpy as np
import scipy
import scipy.stats
import pandas

import pysam
import pybedtools
import rnaseqlib
import rnaseqlib.stats.stats_utils as stats_utils
import rnaseqlib.utils as utils

from collections import defaultdict


def output_bam_coverage_per_base(bam_fname, gff_fname, output_dir,
                                 num_coverage_fields=11):
    """
    Output BAM file's coverage of the given GFF file.
    Example: BAM file of RNA-Seq sample, with GFF file corresponding
    to exons. Returns the filename.

    Args:
    - bam_fname: BAM filename
    - gff_fname: GFF filename of coordinates
    - output_dir: output directory
    """
    if not os.path.isfile(bam_fname):
        raise Exception, "Not a BAM file: %s" %(bam_fname)
    elif not os.path.isfile(gff_fname):
        raise Exception, "Not a GFF file: %s" %(gff_fname)
    utils.make_dir(output_dir)
    bam_file = pybedtools.BedTool(bam_fname)
    coverage_bed = bam_file.coverage(gff_fname, hist=True)
    # Create BED file with only valid number of fields
    #filtered_coverage_bed = pybedtools.BedTool(valid_bed_lines(coverage_bed))
#    filtered_coverage_bed = \
#      coverage_bed.filter(lambda x: len(x.fields) == num_coverage_fields)
    filtered_coverage_bed = coverage_bed.remove_invalid()
    bam_basename = os.path.basename(bam_fname)
    gff_basename = os.path.basename(gff_fname)
    coverage_output_dir = os.path.join(output_dir, "per_base_coverage")
    utils.make_dir(coverage_output_dir)
    # Output coverage results as BED
    output_fname = \
      os.path.join(coverage_output_dir,
                   "%s.%s.coverage.bed" %(bam_basename,
                                          gff_basename))
    print "Outputting BAM coverage per GFF base.."
    print "  - Output file: %s" %(output_fname)
    if os.path.isfile(output_fname):
        # Don't overwrite results if they exist
        return output_fname
    filtered_coverage_bed.saveas(output_fname)
    return output_fname


            

# def output_bam_coverage_per_exon(bam_fname, exons_gff_fname, output_dir,
#                               group_cols=[9],
#                               op_cols=[11, 11, 11],
#                               ops=["mean", "stdev", "max"],
#                               full=True):
#     """
#     Get BAM file's coverage of the given GFF file, and then
#     group the columns (with groupBy) to get relevant statistics.
#     """
#     print "Getting BAM coverage per exon..."
#     print "  - BAM: %s" %(bam_fname)
#     print "  - Exons GFF: %s" %(exons_gff_fname)
#     print "  - Output dir: %s" %(output_dir)
#     # Get the coverage BED file for the exons (map BAM reads
#     # against the exons in the GFF)
#     coverage_bed_fname = output_bam_coverage_per_base(bam_fname, exons_gff_fname)
#     # Command to run:
#     # "groupBy -i coverage_d.bed -g 9 -c 11,11,11 -o mean,stdev,max -full "
#     grouped_bed = \
#       coverage_bed.groupby(g=group_cols, c=op_cols, o=ops, full=full)
#     output_fname = os.path.basename(bam_fname)


def output_exons_coverage_from_tagBam(bam_fname, output_dir,
                                      interval_label="gff"):
    """
    Count exons coverage from BAM file produced by tagBam.
    """
    utils.make_dir(output_dir)
    bam_basename = os.path.basename(bam_fname)
    output_fname = \
      os.path.join(output_dir, "%s_coverage.txt" %(bam_basename))
    print "Outputting exon coverage statistics from tagBam..."
    print "  - BAM: %s" %(bam_fname)
    print "  - Output file: %s" %(output_fname)
    # Mapping from exons to start positions counts, i.e.
    #  exon1 -> genomic start_pos1 -> # at genomic start position 1
    #        -> genomic start_pos4 -> # at genomic start position 4
    #  exon2 -> ...
    exons_to_start_pos_counts = defaultdict(lambda: defaultdict(int))
    bam_file = pysam.Samfile(bam_fname, "rb")
    for bam_read in bam_file:
        gff_aligned_regions = bam_read.opt("YB")
        parsed_regions = \
          gff_aligned_regions.split("%s:" %(interval_label))[1:]
        # Start position of read
        read_start_pos = bam_read.pos
        # For each exon that the read aligns to, record a +1
        # in the genomic start position of read
        for region in parsed_regions:
            # Parse "chr13:21272363-21272783,exon,.,+" into
            # "chr13:21272363-21272783:+"
            region_fields = region.split(",")
            exon = "%s:%s" %(region_fields[0], region_fields[-1])
            exons_to_start_pos_counts[exon][read_start_pos] += 1
    # Calculate statistics and turn into dataframe
    exon_stats_df = []
    for exon in exons_to_start_pos_counts:
        exon_info = exons_to_start_pos_counts[exon]
        exon_counts = exon_info.values()
        entry = {"exon": exon,
                 "mean": np.mean(exon_counts),
                 "std": np.std(exon_counts),
                 "max": np.max(exon_counts),
                 "min": np.min(exon_counts),
                 "kurtosis": scipy.stats.kurtosis(exon_counts),
                 "cv": stats_utils.coeff_var(exon_counts)}
        exon_stats_df.append(entry)
    exon_stats_df = pandas.DataFrame(exon_stats_df)
    output_cols = ["exon", "mean", "max", "std", "cv", "kurtosis"]
    exon_stats_df.to_csv(output_fname,
                         sep="\t",
                         cols=output_cols,
                         index=False,
                         float_format="%.5f")
    return exon_stats_df


def main():
#    bam_fname = "/lab/solexa_jaenisch/solexa_jaenisch2/yarden/Musashi-seq/ribo-pipeline-output/ribo_oe1/mapping/ribo_oe1_KH2MSI1_DOX/processed_bams/accepted_hits.ribosub.sorted.bam"
    bam_fname = "/home/yarden/jaen/Musashi-seq/ribo-pipeline-output/ribo_oe1/analysis/rpkm/ribo_oe1_KH2MSI1_DOX/bam2gff_const_exons/bam2gff_ensGene.cds_only.const_exons.gff/accepted_hits.ribosub.sorted.bam"
    gff_fname = "/home/yarden/jaen/pipeline_init/mm9/ucsc/exons/ensGene.cds_only.exons.gff"
    output_dir = "/home/yarden/jaen/rnaseqlib/rnaseqlib/coverage/"
#    output_bam_coverage_per_base(bam_fname, gff_fname, output_dir)
#    output_bam_coverage_per_exon(bam_fname, gff_fname, output_dir)
    small_bam_fname = "/home/yarden/jaen/Musashi-seq/ribo-pipeline-output/ribo_oe1/analysis/rpkm/rpkm_trim27_MSI1DOX.bam"
    df = output_exons_coverage_from_tagBam(bam_fname, "./test")
    df.to_csv("./test.txt",sep="\t")
    cols = ["cv", "mean", "max"]
    for c in cols:
        print "mean of %s" %(c)
        print "  %.4f" %(df[c].mean())
        
    ##
    ## TODO:
    ## New strategy: look at tagBam output, where reads are mapped to exons
    ## then make a dictionary mapping exons to their reads?  No.
    ##

if __name__ == "__main__":
    main()
