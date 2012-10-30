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

import pandas

import pysam


def load_sample_rpkms(sample,
                      rna_base):
    """
    Load RPKM tables for a sample. Return
    """
    rpkm_tables = {}
    for table_name in rna_base.rpkm_table_names:
        rpkm_table = None
        rpkm_filename = os.path.join(sample.rpkm_dir,
                                     "%s.rpkm" %(table_name))
        if os.path.isfile(rpkm_filename):
            # Insert the sample name into the header
            fieldnames = ["gene_id",
                          "rpkm_%s" %(sample.label),
                          "counts_%s" %(sample.label),
                          "exons"]
            # Load each table as a DataFrame
            rpkm_table = pandas.read_csv(rpkm_filename,
                                         sep="\t",
                                         names=fieldnames,
                                         # Skip the current header
                                         skiprows=1)
            print "rpkm_table: ", rpkm_table
            # Add gene_symbol and gene_desc columns
            # to RPKM DataFrame
            gene_table = rna_base.gene_tables[table_name.split(".")[0]]
            gene_symbols = [gene_table.genes_to_names[gid] \
                            for gid in rpkm_table["gene_id"]]
            gene_descs = [gene_table.genes_to_desc[gid] \
                          for gid in rpkm_table["gene_id"]]
            rpkm_table["gene_symbol"] = gene_symbols
            rpkm_table["gene_desc"] = gene_descs
        else:
            print "WARNING: Cannot find RPKM filename %s" %(rpkm_filename)
        rpkm_tables[table_name] = rpkm_table
    return rpkm_tables
    

def output_rpkm(sample,
                output_dir,
                settings_info,
                rna_base,
                logger):
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
            logger.info("  - Skipping RPKM output, found %s" %(rpkm_output_filename))
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
        num_mapped = int(sample.qc.qc_results["num_mapped"])
        if num_mapped == 0:
            print "Error: Cannot compute RPKMs since sample %s has 0 mapped reads." \
                %(sample.label)
            sys.exit(1)
        print "Sample %s has %s mapped reads" %(sample.label, num_mapped)
        read_len = settings_info["readlen"]
        logger.info("Outputting RPKM from GFF aligned BAM")
        output_rpkm_from_gff_aligned_bam(exons_bam_fname,
                                         num_mapped,
                                         read_len,
                                         const_exons,
                                         rpkm_output_filename)
    return rpkm_output_filename
    
    
def output_rpkm_from_gff_aligned_bam(bam_filename,
                                     num_mapped,
                                     read_len,
                                     const_exons,
                                     output_filename,
                                     rpkm_header=["gene_id",
                                                  "rpkm",
                                                  "counts",
                                                  "exons"],
                                     na_val="NA"):
    """
    Given a BAM file aligned by bedtools (with 'gff' field),
    compute RPKM for each region, incorporating relevant
    optional fields from gff.

    Takes as input:

     - bam_filename: the BAM file
     - num_mapped: number of mapped reads to normalize to
     - read_len: read length
     - const_exons: Constitutive exons object
     - output_filename: output filename
    """
    bam_file = pysam.Samfile(bam_filename, "rb")
    print "Computing RPKM from BAM aligned to GFF..."
    print "  - BAM: %s" %(bam_filename)
    print "  - Output filename: %s" %(output_filename)
    # Map of gff region to read counts
    region_to_count = defaultdict(int)
    for bam_read in bam_file:
        try:
            # Read aligns to region of interest
            gff_aligned_regions = bam_read.opt("YB")
            parsed_regions = gff_aligned_regions.split("gff:")[1:]
            # Compile region counts and lengths
            for region in parsed_regions:
                region_chrom, coord_field = region.split(",")[0].split(":")[0:2]
                # Region internally converted to 0-based start, so we must add 1
                # to get it back
                region_start, region_end = map(int, coord_field.split("-"))
                region_start += 1
                region_str = "%s:%s-%s" %(region_chrom,
                                          str(region_start),
                                          str(region_end))
                # Count reads in region
                region_to_count[region_str] += 1
        except KeyError:
            gff_aligned_region = None
    # For each gene, find its exons. Sum their counts
    # and length to compute RPKM
    rpkm_table = []
    for gene_info in const_exons.genes_to_exons:
        gene_id = gene_info["gene_id"]
        exons = gene_info["exons"]
        if exons == na_val:
            continue
        parsed_exons = exons.split(",")
        # Strip the strand of the exons
        strandless_exons = map(lambda x: x[0:-2], parsed_exons)
        curr_counts = [region_to_count[s_exon] for s_exon in strandless_exons]
        sum_counts = sum(curr_counts)
        curr_lens = [const_exons.exon_lens[exon] for exon in parsed_exons]
        sum_lens = sum(curr_lens)
        assert(len(curr_counts) == len(curr_lens)), \
            "Error: sum_counts != sum_lens in RPKM computation."
        gene_rpkm = compute_rpkm(sum_counts, sum_lens, num_mapped)
        # RPKM entry for gene
        rpkm_entry = {"rpkm": gene_rpkm,
                      "gene_id": gene_id,
                      "counts": sum_counts,
                      "exons": exons}
        rpkm_table.append(rpkm_entry)
    rpkm_df = pandas.DataFrame(rpkm_table)
    rpkm_df.to_csv(output_filename,
                   cols=rpkm_header,
                   na_rep=na_val,
                   sep="\t",
                   # 4-decimal point RPKM format
                   # Not compatible with current pandas versions
                   #float_format="%.4f",
                   index=False)
    return output_filename


def compute_rpkm(region_count,
                 region_len,
                 num_total_reads):
    """
    Compute RPKM for a region.
    """
    # Get length of region in KB
    region_kb = region_len / float(1e3)

    # Numerator of RPKM: reads per kilobase
    rpkm_num = (region_count / region_kb)

    # Denominator of RPKM: per M mapped reads
    num_reads_per_million = num_total_reads / float(1e6)

    rpkm = (rpkm_num / num_reads_per_million)
    return rpkm
            
