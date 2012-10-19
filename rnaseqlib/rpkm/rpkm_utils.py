##
## Utilities for computing RPKM
##
import os
import sys
import time

from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.exon_utils as exon_utils

import pysam

def output_rpkm(sample,
                output_dir,
                settings_info,
                rna_base):
    """
    Output RPKM tables for the sample.

    Takes as input:

    - sample: a sample object
    - output_dir: output directory
    - settings_info: settings information
    - rna_base: an RNABase object
    """
    # Output RPKM information for all constitutive exon tables in the
    # in the RNA Base
    print "Outputting RPKM for: %s" %(sample.label)
    rpkm_tables = {}
    for table_name, const_exons in rna_base.tables_to_const_exons.iteritems():
        rpkm_output_filename = "%s.rpkm" %(os.path.join(output_dir,
                                                        table_name))
        rpkm_tables[table_name] = rpkm_output_filename
        if os.path.isfile(rpkm_output_filename):
            print "  - Skipping RPKM output, %s exists" %(rpkm_output_filename)
            continue
        # Directory where BAM containing mapping to constitutive
        # exons be stored
        bam2gff_outdir = os.path.join(output_dir,
                                      "bam2gff_const_exons")
        utils.make_dir(bam2gff_outdir)
        # Map reads to GFF of constitutive exons
        exons_bam_fname = exon_utils.map_bam2gff(sample.bam_filename,
                                                 const_exons.gff_filename,
                                                 bam2gff_outdir)
        # Compute RPKMs for sample
        num_mapped = sample.qc.qc_results["num_mapped"]
        if num_mapped == 0:
            print "Error: Cannot compute RPKMs since sample %s has 0 mapped reads." \
                %(sample.label)
            sys.exit(1)
        print "Sample %s has %s mapped reads" %(sample.label, num_mapped)
        read_len = settings_info["readlen"]
        output_rpkm_from_gff_aligned_bam(exons_bam_fname,
                                         num_mapped,
                                         read_len,
                                         rpkm_output_filename)
    return rpkm_output_filename
    
    
def output_rpkm_from_gff_aligned_bam(bam_filename,
                                     num_total_reads,
                                     read_len,
                                     output_filename):
    """
    Given a BAM file aligned by bedtools (with 'gff' field),
    compute RPKM for each region, incorporating relevant
    optional fields from gff.
    """
    bam_file = pysam.Samfile(bam_filename, "rb")
    source = os.path.basename(output_filename)

    print "Computing RPKM from BAM aligned to GFF..."
    print "  - BAM: %s" %(bam_filename)
    print "  - Output filename: %s" %(output_filename)

    loaded_gff = False
    ref_gff_recs = None
    last_chrom = None

    # Map of gff region to read counts
    region_to_count = defaultdict(int)
    region_to_len = defaultdict(int)

    for bam_read in bam_file:
        curr_chrom = bam_file.getrname(bam_read.tid)
        try:
            # Read aligns to region of interest
            gff_aligned_regions = bam_read.opt("YB")
            parsed_regions = gff_aligned_regions.split("gff:")[1:]
            # Compile region counts and lengths
            for region in parsed_regions:
                region_to_count[region] += 1

                # Get region length
                coord_field = region.split(",")[0].split(":")[1]
                region_start, region_end = coord_field.split("-")
                region_start, region_end = int(region_start), \
                                           int(region_end)
                region_len = region_end - region_start
                region_to_len[region] = region_len
        except KeyError:
            gff_aligned_region = None
#    print "region_to_count: ", region_to_count
#    print "region to len: ", region_to_len

#    output_rpkm_as_gff(source, region_to_count, region_to_len,
#                       num_total_reads, output_filename)
