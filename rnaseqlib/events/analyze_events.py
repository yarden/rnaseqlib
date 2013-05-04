##
## Compare new events with old
##
import os
import sys
import time

from collections import defaultdict

import misopy
import misopy.Gene as gene_utils

new_se_gff = \
    os.path.expanduser("~/jaen/new-gff-events/splicingDB/mm9/commonshortest/SE.gff3")
old_se_gff = \
    os.path.expanduser("~/jaen/gff-events/mm9/SE.mm9.gff3")


def index_exons(gff_fname):
    gff_genes = gene_utils.load_genes_from_gff(gff_fname)
    exons = defaultdict(bool)
    for gene_id in gff_genes:
        gene_obj = gff_genes[gene_id]["gene_object"]
        se = gene_obj.isoforms[0].parts[1]
        # Index the exon by chromosome
        exons[(gene_obj.chrom, se.start, se.end, gene_obj.strand)] = True
    return exons


def get_gff_fnames(genome, event_type):
    old_gff_fname = \
        os.path.expanduser("~/jaen/gff-events/mm9/%s.mm9.gff3" \
                           %(event_type))
    new_gff_fname = \
        os.path.expanduser("~/jaen/new-gff-events/splicingDB/mm9/commonshortest/%s.gff3" \
                           %(event_type))
    return old_gff_fname, new_gff_fname


def compare_exons(genome, event_type):
    print "Comparing exons..."
    print "  GENOME: %s" %(genome)
    print "  EVENT TYPE: %s" %(event_type)
    # Old and new GFF files
    first_gff_fname, second_gff_fname = get_gff_fnames(genome, event_type)
    first_exons = index_exons(first_gff_fname)
    second_exons = index_exons(second_gff_fname)
    # Compute the overlap between exons
    first_exon_set = set(first_exons.keys())
    second_exon_set = set(second_exons.keys())
    # Number common
    num_intersect = len(first_exon_set.intersection(second_exon_set))
    # Number in first but not second
    num_in_first_not_second = len(first_exon_set.difference(second_exon_set))
    # Number in second but not first
    num_in_second_not_first = len(second_exon_set.difference(first_exon_set))
    print "%d exons in both" %(num_intersect)
    print "%d exons in old but not new" %(num_in_first_not_second)
    print "%d exons in new but not old" %(num_in_second_not_first)
    # Now compute partial overlaps: out of exons that are in first
    # but not second,
    unique_to_first = first_exon_set.difference(second_exon_set)
    unique_to_first_starts = {}
    unique_to_first_ends = {}
    for f_exon in unique_to_first:
        f_exon_chrom = f_exon[0]
        f_exon_start = f_exon[1]
        f_exon_end = f_exon[2]
        unique_to_first_starts[(f_exon_chrom, f_exon_start)] = True
        unique_to_first_ends[(f_exon_chrom, f_exon_end)] = True
    # Check which second exons share a start or end with the
    # exons unique to first
    num_second_share_start_with_first = 0
    num_second_share_end_with_first = 0
    for s_exon in second_exon_set:
        # Check for shared start
        if (s_exon[0], s_exon[1]) in unique_to_first_starts:
            num_second_share_start_with_first += 1
        if (s_exon[0], s_exon[2]) in unique_to_first_ends:
            num_second_share_start_with_first += 1
    print "num new share start with old: %d" \
          %(num_second_share_start_with_first)
    print "num new share end with old: %d" \
          %(num_second_share_end_with_first)


def main():
    for genome in ["mm9"]:
        for event_type in ["SE", "A3SS", "A5SS", "RI", "MXE"]:
            compare_exons(genome, event_type)
    

if __name__ == "__main__":
    main()
