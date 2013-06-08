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
    # Compute lengths of UTRs
    # Mapping from gene ID to a dictionary mapping each
    # UTR to its length
    genes_to_utr_lens = defaultdict(lambda: defaultdict(int))
    print "Computing lengths of UTRs.."
    gff_utrs = pybedtools.BedTool(utr_gff_fname)
    for entry in gff_utrs:
        # Get transcript that UTR belongs to
        trans_id = entry.attrs["Parent"]
        # Get the gene it corresponds to
        gene_id = trans_to_gene[trans_id]
        print "GENE_ID: %s" %(gene_id), " for ", trans_id
        # Compute length of transcript
        
        

def main():
    # Load tables from the init dir.
    tables_dir = os.path.expanduser("~/jaen/test/mm9/ucsc/")
    utr_fname = os.path.expanduser("~/jaen/pipeline_init/mm9/ucsc/utrs/test.gff")
    output_fname = "./testtable.gff"
    output_utr_table(tables_dir, utr_fname, output_fname)


if __name__ == "__main__":
    main()
    
