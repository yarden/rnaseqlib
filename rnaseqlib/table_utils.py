##
## Utilities for manipulating init tables
##
import os
import sys
import time

from collections import defaultdict

import pandas

import pybedtools

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.tables as tables

    
def output_utr_table(tables_dir,
                     utr_gff_fname,
                     output_fname,
                     choice_rule="longest"):
    """
    Output a UTR table (one UTR per gene) given a
    UTR GFF file. Possible rules for choosing the
    UTR ('choice_rule'):

      - longest, uses longest UTR
      - shortest, uses shortest UTR

    Outputs a GFF file.
    """
    print "Outputting UTR table from %s" %(utr_gff_fname)
    print "  - Output file: %s" %(output_fname)
    if not os.path.isfile(utr_gff_fname):
        raise Exception, "Cannot find %s" %(utr_gff_fname)
    # Load table
    table_fname = os.path.join(tables_dir, "ensGene.kgXref.combined.txt")
    table_df = pandas.read_table(table_fname, sep="\t")
    trans_to_gene = {}
    # Map transcripts to genes
    for row, entry in table_df.iterrows():
        trans_to_gene[entry["name"]] = entry["name2"]
    # Mapping from gene ID to a dictionary mapping each
    # UTR to its length
    genes_to_utr_lens = defaultdict(lambda: defaultdict(int))
    print "Computing lengths of UTRs.."
    gff_utrs = pybedtools.BedTool(utr_gff_fname)
    # Compute lengths of UTRs for each gene
    for entry in gff_utrs:
        # Get transcript that UTR belongs to
        trans_id = entry.attrs["Parent"]
        # Get UTR id
        utr_id = entry.attrs["ID"]
        # Get the gene it corresponds to
        gene_id = trans_to_gene[trans_id]
        # Compute length of UTRs
        # Length of UTR
        utr_len = len(entry)
        genes_to_utr_lens[gene_id][utr_id] = utr_len
    # Select UTR for each gene
    for gene in genes_to_utr_lens:
        all_utrs = genes_to_utr_lens[gene].items()
        utr_lens = [curr_utr[1] for curr_utr in all_utrs]
        print "UTR LENS OF GENE: ", utr_lens, " gene: ", gene
        if choice_rule == "longest":
            print "UTR LENS: ", utr_lens
            utr_indx = utils.max_item(utr_lens)[0]
            print " CHOSE: ", all_utrs[utr_indx]
        else:
            raise Exception, "Unsupported choice rule %s" %(choice_rule)
    print "genes_to_utr_lens: "
    print genes_to_utr_lens
        

def main():
    # Load tables from the init dir.
    tables_dir = os.path.expanduser("~/jaen/test/mm9/ucsc/")
    utr_fname = os.path.expanduser("~/jaen/pipeline_init/mm9/ucsc/utrs/test.gff")
    output_fname = "./testtable.gff"
    output_utr_table(tables_dir, utr_fname, output_fname)


if __name__ == "__main__":
    main()
    
