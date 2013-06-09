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
    gene_to_chosen_utr = {}
    for gene in genes_to_utr_lens:
        all_utrs = genes_to_utr_lens[gene].items()
        utr_lens = [curr_utr[1] for curr_utr in all_utrs]
        if choice_rule == "longest":
            utr_indx = utils.max_item(utr_lens)[0]
            chosen_utr = all_utrs[utr_indx]
            gene_to_chosen_utr[gene] = chosen_utr
        else:
            raise Exception, "Unsupported choice rule %s" %(choice_rule)
    # Now select the relevant entries for outputting. Also
    # add relevant information about genes/length
    gff_utrs = pybedtools.BedTool(utr_gff_fname)
    gff_out = open(output_fname, "w")
    for entry in gff_utrs:
        # Current UTR id
        curr_utr_id = entry.attrs["ID"]
        # Current UTR's transcript
        curr_utr_trans = entry.attrs["Parent"]
        # Get the current UTR's gene
        curr_utr_gene = trans_to_gene[curr_utr_trans]
        # If this UTR is the chosen UTR, output it
        if gene_to_chosen_utr[curr_utr_gene][0] == curr_utr_id:
            # Look up the gene ID it belongs to
            curr_gene_id = trans_to_gene[curr_utr_trans]
            entry.attrs["gene_id"] = curr_gene_id
            entry.attrs["utr_len"] = \
                str(gene_to_chosen_utr[curr_utr_gene][1])
            gff_out.write("%s" %(str(entry)))
    gff_out.close()
    
        

def main():
    # Load tables from the init dir.
    tables_dir = os.path.expanduser("~/jaen/test/mm9/ucsc/")
    utr_fname = os.path.expanduser("~/jaen/pipeline_init/mm9/ucsc/utrs/ensGene.3p_utrs.gff")
    #output_fname = "./testtable.gff"
    #output_utr_table(tables_dir, utr_fname, output_fname)
    output_fname = "./smallutrs.gff"
    utrs = pybedtools.BedTool(output_fname)
    


if __name__ == "__main__":
    main()
    
