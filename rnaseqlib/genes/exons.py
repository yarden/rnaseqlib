##
## Utilities for working with exons
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.gff_utils as gff_utils

def exon_overlaps_mRNA(gff_in, exon, mRNA, base_diff):
    """
    Return True if exon is in the mRNA, allowing 'base_diff'-many
    bases to be off at the start or end of the exon.
    
    - exon: gff record of exon
    - mRNA: gff record of mRNA
    """
    mRNA_id = mRNA.get_id()
    if mRNA_id not in gff_in.exons_by_mRNA:
        return False
    for curr_exon in gff_in.exons_by_mRNA[mRNA.get_id()]:
        start_diff = abs(curr_exon.start - exon.start)
        end_diff = abs(curr_exon.end - exon.end)
        # Ensure same stranded-ness
        if curr_exon.strand != exon.strand: continue
        # If the exon start differences and exon end
        # differences are 'base_diff' or less, then
        # the exon is considered overlapping the mRNA
        if (start_diff <= base_diff) and \
           (end_diff <= base_diff):
            return True
    return False
    

def const_exons_from_mRNAs(gff_in, mRNAs,
                           base_diff=0):
    """
    optional:

    - all_constitutive: flag to treat all exons as constitutive
    """
    const_exons = []
    # Get first mRNA's exons
    gene_id = mRNAs[0].get_parent()
    mRNA_id = mRNAs[0].get_id()
    if mRNA_id not in gff_in.exons_by_mRNA:
        # If mRNA has no exons, skip
        return const_exons

    exons = gff_in.exons_by_mRNA[mRNA_id]

    for exon in exons:
        const_exon = True
        # Skip exons that don't meet size requirement
        exon_len = exon.end - exon.start + 1

        # Exon flagged constitutive unless its overlap diff
        # with *all* exons in the mRNAs is greater than base_diff
        for mRNA in mRNAs[1:]:
            curr_mRNA_id = mRNA.get_id()
            if not exon_overlaps_mRNA(gff_in, exon, mRNA,
                                      base_diff):
                const_exon = False
                break
        # If exon is constitutive, add the parent gene information
        # as a field
        exon.attributes['GeneParent'] = [gene_id]
        # Record exon if it is constitutive
        if const_exon: const_exons.append(exon)
    return const_exons


def get_const_exons(gff_filename, output_filename,
                    base_diff=5):
    """
    Get constitutive exons for GFF filename.

    - base_diff: Number of bases +/- that can be omitted when
      for an exon to be considered constitutive.
    """
    print "Getting constitutive exons from: %s" %(gff_filename)
    dir_name = os.path.dirname(output_filename)
    if not os.path.isdir(dir_name):
        utils.make_dir(dir_name)
    print "Loading GFF file..."
    gff_db = gff_utils.GFFDatabase(from_filename=gff_filename,
                                   reverse_recs=True)
    print "Done loading."
    gff_out = gff_utils.GFFWriter(open(output_filename, "w"))
    
    for gene, mRNAs in gff_db.mRNAs_by_gene.iteritems():
        # Get constitutive exons from the current set
        # of mRNAs
        const_exons = const_exons_from_mRNAs(gff_in, mRNAs)
        for exon_rec in const_exons:
            # Write exons to file
            gff_out.write_rec(exon_rec)
        
