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
import rnaseqlib.bam.bam_utils as bam_utils
import rnaseqlib.utils as utils

from collections import defaultdict, OrderedDict


# def output_bam_coverage_per_base(bam_fname, gff_fname, output_dir,
#                                  num_coverage_fields=11):
#     """
#     Output BAM file's coverage of the given GFF file.
#     Example: BAM file of RNA-Seq sample, with GFF file corresponding
#     to exons. Returns the filename.

#     Args:
#     - bam_fname: BAM filename
#     - gff_fname: GFF filename of coordinates
#     - output_dir: output directory
#     """
#     if not os.path.isfile(bam_fname):
#         raise Exception, "Not a BAM file: %s" %(bam_fname)
#     elif not os.path.isfile(gff_fname):
#         raise Exception, "Not a GFF file: %s" %(gff_fname)
#     utils.make_dir(output_dir)
#     bam_file = pybedtools.BedTool(bam_fname)
#     coverage_bed = bam_file.coverage(gff_fname, hist=True)
#     # Create BED file with only valid number of fields
#     #filtered_coverage_bed = pybedtools.BedTool(valid_bed_lines(coverage_bed))
# #    filtered_coverage_bed = \
# #      coverage_bed.filter(lambda x: len(x.fields) == num_coverage_fields)
#     filtered_coverage_bed = coverage_bed.remove_invalid()
#     bam_basename = os.path.basename(bam_fname)
#     gff_basename = os.path.basename(gff_fname)
#     coverage_output_dir = os.path.join(output_dir, "per_base_coverage")
#     utils.make_dir(coverage_output_dir)
#     # Output coverage results as BED
#     output_fname = \
#       os.path.join(coverage_output_dir,
#                    "%s.%s.coverage.bed" %(bam_basename,
#                                           gff_basename))
#     print "Outputting BAM coverage per GFF base.."
#     print "  - Output file: %s" %(output_fname)
#     if os.path.isfile(output_fname):
#         # Don't overwrite results if they exist
#         return output_fname
#     filtered_coverage_bed.saveas(output_fname)
#     return output_fname


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



def add_summary_coverage_col_to_df(df, col,
                                   op="max",
                                   na_val="NA"):
    """
    Summarize a coverage column 'col' in the given
    dataframe 'df', e.g. by taking its max across exons.

    Kwargs:
    - op: the operation (e.g. "max") to summarize the coverage
    column across exons.
    - na_val: NA value to assume in input exon values in df
    """
    values = []
    val_to_keep = na_val
    for row_num, row_entry in df.iterrows():
        if not pandas.notnull(row_entry[col]):
            val_to_keep = np.nan
        else:
            curr_values = row_entry[col].split(",")
            curr_values = \
              [np.nan if elt == na_val else float(elt) \
               for elt in curr_values]
            if op == "max":
                val_to_keep = max(curr_values)
        values.append(val_to_keep)
    # Add column to df
    summary_col = "max_%s" %(col)
    df[summary_col] = values
    return df, summary_col


def init_region_start_pos(region, ordered=True):
    """
    Given a region as string, 'chr13:21281961-21281993:+',
    return a dictionary of size length of the region.

    Assumes GFF coordinate conventions.

    Args:
      - region: region as string, e.g. 'chr13:21281961-21281993:+'

    Kwargs:
      - ordered: if True, use OrderedDict

    Returns:
      dictionary mapping start positions to 0
    """
    chrom, start, end, strand = utils.parse_dash_coords(region)
    # Check that start < end
    assert (start <= end), \
      "Start (%d) must be less than end (%d): %s" %(start, end, region)
    # Check strand
    if not (strand == "+" or strand == "-"):
        raise Exception, "Unknown strand symbol %s" %(str(strand))
    # Region length according to GFF
    region_len = end - start + 1
    # Dictionary mapping region starts to counts
    region_pos = {}
    if ordered:
        region_pos = OrderedDict()
    for pos in range(start, end + 1):
        region_pos[pos] = 0
    # Check that the length of the start positions dict equals the
    # length of the region (assuming GFF coordinates)
    pos_dict_len = len(region_pos)
    assert (region_len == pos_dict_len), \
      "Error: Region length (%d) does not match start pos " \
      "dictionary length (%d)" %(region_len, pos_dict_len)
    return region_pos


def get_exons_coverage_from_tagBam(bam_fname,
                                   interval_label="gff",
                                   gff_coords=True):
    """
    Count exons coverage from BAM file produced by tagBam.
    Returns mapping from an exon to its coverage statistics.

    Args:
    - bam_fname: BAM filename produced by tagBam

    Kwargs:
    - interval_label: interval label in optional field of tagBAM
    that identifies the interval from GFF/BED that BAM was mapped
    against
    - gff_coords: if True, then convert the intervals from the BAM
    from BED format (0-based) to GFF based by adding 1 to the start
    coordinate. Note that tagBam produces intervals in BED format
    always.
    """
    # Mapping from exons to start positions counts, i.e.
    #  exon1 -> genomic start_pos1 -> # at genomic start position 1
    #        -> genomic start_pos4 -> # at genomic start position 4
    #  exon2 -> ...

    ###
    ### TODO: need to account for zero coverage positions, which
    ### default dict does not do!
    #exons_to_start_pos_counts = defaultdict(lambda: defaultdict(int))
    exons_to_start_pos_counts = {}
    
    bam_file = pysam.Samfile(bam_fname, "rb")
    num_unmatched = 0
    for bam_read in bam_file:
        gff_aligned_regions = bam_read.opt("YB")
        # Get GFF regions that read aligns to
        parsed_regions = \
          bam_utils.parse_tagBam_opt_field(gff_aligned_regions,
                                           gff_coords=gff_coords)
        # Get read start position
        read_start_pos = bam_read.pos + 1
        for region in parsed_regions:
            # If we've seen this exon region before, use its dictionary
            if region in exons_to_start_pos_counts:
                curr_region_starts = exons_to_start_pos_counts[region]
            else:
                curr_region_starts = init_region_start_pos(region)
                exons_to_start_pos_counts[region] = curr_region_starts
            # Add +1 to each read position that the read overlaps
            # Have a simple read position counter that starts with
            # read start and counts the bases covered by a Match
            # according to the cigar string
            read_counter = read_start_pos
            read_matched_to_interval = False
            for cigar_type, cigar_len in bam_read.cigar:
                # Skip the non-matches (denoted as 0 in pysam)
                # but increment the read_counter
                if cigar_type != 0:
                    read_counter += cigar_len
                    continue
                # It's a match, so record the bases covered
                for match_pos in range(read_counter,
                                       read_counter + cigar_len + 1):
                    # This match portion of read does not land in region
                    if match_pos not in curr_region_starts:
                        # Advance counter
                        read_counter += 1
                        continue
                    curr_region_starts[match_pos] += 1
                    # Record that the read was matched to this interval
                    read_matched_to_interval = True
            if not read_matched_to_interval:
                num_unmatched += 1
#                print "Error: read ", bam_read, " never matched to: ", \
#                      region, "!"
#                print "Current region starts: "
#                print curr_region_starts.keys()
#                print max(curr_region_starts.values())
#                raise Exception, "Test"
 #   print "Total number of reads unmatched: %d" %(num_unmatched)
#    raise Exception, "Test"
    # Calculate statistics and return as dictionary
    # indexed by exons
    print "PRINTING VALUES FOR EXONS: "
    ex = ["chr13:21281961-21281993:+",
          "chr13:21272364-21272783:+"]
    ###
    ### TODO: here, try out several measures for the two exons, including
    ### kurtosis, CV, and some measure using entropy.
    ###
    ### Can use KL divergence to uniform distribution, total entropy,
    ### or JSD
    ###
    for e in ex:
        print "e: %s" %(e)
        print "-" * 10
        # Number of reads covering each base in exon
        exon_reads_per_base = exons_to_start_pos_counts[e].values()
        cv_val = stats_utils.coeff_var(exon_reads_per_base)
        kurtosis_val = scipy.stats.kurtosis(exon_reads_per_base,
                                            fisher=False)
        # Calculate square root of JSD to uniform distribution
        exon_len = len(exon_reads_per_base)
        # Compare to uniform distribution
        uniform_dist = np.array([1 / float(exon_len)] * exon_len)
        # Actual coverage per base: number of reads at base
        # divided by total number of reads
        total_reads = np.sum(exon_reads_per_base)
        observed_dist = np.array(exon_reads_per_base) / float(total_reads)
        sqrt_jsd_val = stats_utils.sqrt_jsd(observed_dist, sqrt_jsd_val)
        print "CV: ", cv_val
        print "kurtosis: ", kurtosis_val
        print "sqrt JSD(observed, uniform): ", jsd_val
        for p in exons_to_start_pos_counts[e]:
            print p, " => ", exons_to_start_pos_counts[e][p]
        print "\n"
    raise Exception, "End Test"
    
    exon_stats_dict = {}
    for exon in exons_to_start_pos_counts:
        exon_info = exons_to_start_pos_counts[exon]
        exon_counts = exon_info.values()
        entry = {"exon": exon,
                 "mean": np.mean(exon_counts),
                 "std": np.std(exon_counts),
                 "max": np.max(exon_counts),
                 "min": np.min(exon_counts),
                 "kurtosis": scipy.stats.kurtosis(exon_counts,
                                                  fisher=False),
                 "cv": stats_utils.coeff_var(exon_counts)}
        exon_stats_dict[exon] = entry
    return exon_stats_dict
#    output_cols = ["exon", "mean", "max", "std", "cv", "kurtosis"]
#    exon_stats_df.to_csv(output_fname,
#                         sep="\t",
#                         cols=output_cols,
#                         index=False,
#                         float_format="%.5f")


def main():
#    bam_fname = "/lab/solexa_jaenisch/solexa_jaenisch2/yarden/Musashi-seq/ribo-pipeline-output/ribo_oe1/mapping/ribo_oe1_KH2MSI1_DOX/processed_bams/accepted_hits.ribosub.sorted.bam"
    bam_fname = "/home/yarden/jaen/Musashi-seq/ribo-pipeline-output/ribo_oe1/analysis/rpkm/ribo_oe1_KH2MSI1_DOX/bam2gff_const_exons/bam2gff_ensGene.cds_only.const_exons.gff/accepted_hits.ribosub.sorted.bam"
    gff_fname = "/home/yarden/jaen/pipeline_init/mm9/ucsc/exons/ensGene.cds_only.exons.gff"
    output_dir = "/home/yarden/jaen/rnaseqlib/rnaseqlib/coverage/"
#    output_bam_coverage_per_base(bam_fname, gff_fname, output_dir)
#    output_bam_coverage_per_exon(bam_fname, gff_fname, output_dir)
    small_bam_fname = "/home/yarden/jaen/Musashi-seq/ribo-pipeline-output/ribo_oe1/analysis/rpkm/rpkm_trim27_MSI1DOX.bam"
    df = get_exons_coverage_from_tagBam(bam_fname, "./test")
    df.to_csv("./test.txt",sep="\t")
    cols = ["cv", "mean", "max"]
    for c in cols:
        print "mean of %s" %(c)
        print "  %.4f" %(df[c].mean())


if __name__ == "__main__":
    main()
