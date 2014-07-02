##
## Utilities for working with PhyloP files from
## UCSC
##
import os
import sys
import time
import glob

import rnaseqlib
import rnaseqlib.utils as utils
import pybedtools


def load_phylop_bed(chrom, phylop_dir):
    """
    Load PhyloP BED file for a given chromosome
    from the given directory, and return it as a
    BedTool (pybedtools).
    """
    if not os.path.isdir(phylop_dir):
        raise Exception, "Not a PhyloP dir %s" %(phylop_dir)
    fname_wildcard = os.path.join(phylop_dir, "%s*.bed" %(chrom))
    # Get the BED file for the chromosome
    candidate_fnames = glob.glob(fname_wildcard)
    # Ignore the '_random' chromosome files
    candidate_fnames = filter(lambda f: "random" not in f,
                              candidate_fnames)
#    candidate_fnames = ["./test.bed"]
    if len(candidate_fnames) == 0:
        raise Exception, "Could not find PhyloP *.bed files for %s" \
                         %(chrom)
    if len(candidate_fnames) > 1:
        raise Exception, "More than one PhyloP *.bed file for %s" \
                         %(chrom)
    phylop_fname = candidate_fnames[0]
    print "Loading %s" %(phylop_fname)
    # Load PhyloP BED
    phylop_bed = pybedtools.BedTool(phylop_fname)
    return phylop_bed


def get_phylop_for_gene(gene_rec,
                        phylop_dir,
                        species="mm9",
                        window=None,
                        window_format="srcwinnum",
                        window_op="mean"):
    """
    Get phyloP score for gene. Assumes there are *.bed
    versions of the PhyloP wigs from UCSC.

    Args:
    - gene_rec: gene record (read by pybedtools), which
    has a chr, start, end and strand fields for the gene.
    - phylop_dir: directory of PhyloP *.bed files.

    Kwargs:
    - window: if given, bin the PhyloP scores into windows
      and get the mean
    - window_format: format argument to 'windowmaker' of bedtools
    - window_op: operation on counts to do when binning. Set to 'mean'
      by default
    """
    gene_fields = gene_rec.fields
    chrom = gene_fields[0]
    phylop_bed = load_phylop_bed(chrom, phylop_dir)
    # Intersect it with the current gene record. Get
    # only PhyloP scores that fall directly within
    # gene interval
    print "Getting PhyloP regions for gene..."
    t1 = time.time()
    gene_bedtool = pybedtools.BedTool(str(gene_rec), from_string=True)
    if window is not None:
        # Create windows for the gene region
        binned_gene_bed = \
          pybedtools.BedTool().window_maker(b=gene_bedtool,
                                            w=window,
                                            i=window_format)
        # Count the features of interest
        region_phylop = binned_gene_bed.map(b=phylop_bed,
                                            o=window_op,
                                            s=True)
    else:
        region_phylop = phylop_bed.intersect(gene_bedtool, f=1)
    print "Region phylop: ", region_phylop
    t2 = time.time()
    print "  - Took %.2f seconds" %(t2 - t1)
    return region_phylop


    

