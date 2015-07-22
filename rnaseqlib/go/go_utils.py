##
## GO related utilities
##
import os
import sys
import time

import rnaseqlib

def get_genes_to_go_terms(gene_assoc_fname,
                          ensembl_to_mgi_fname,
                          mgi_col=1,
                          go_term_col=4):
    """
    Return a mapping (dictionary) from gene IDs to their GO terms.

    E.g. map ENSMUSG... -> GO:XXXX

    Args:
    - gene_assoc_fname: gene association filename from GO (*.mgi file)
    - ensembl_to_mgi_fname: mapping from Ensembl ID to MGI IDs 
    - mgi_col: the column number (0-based) in gene_assoc_fname
      that contains the MGI ID.
    """
    gene_to_go_terms = {}
    # Load mapping from MGI to Ensembl
    mgi_to_ensembl = {}
    with open(ensembl_to_mgi_fname) as mapping_in:
        # Skip header
        mapping_in.readline()
        for line in mapping_in:
            fields = line.split("\t")
            mgi_to_ensembl[fields[1]] = fields[0]
    not_found = 0
    with open(gene_assoc_fname) as go_in:
        for line in go_in:
            if line.startswith("!"):
                # Skip GO comments
                continue
            fields = line.split("\t")
            if len(fields) < mgi_col:
                continue
            mgi_id = fields[mgi_col]
            go_term = fields[go_term_col]
            # Retrieve Ensembl ID for current MGI ID
            if mgi_id not in mgi_to_ensembl:
                not_found += 1
                continue
            curr_gene_id = mgi_to_ensembl[mgi_id]
            gene_to_go_terms[curr_gene_id] = go_term
    print "Total of %d genes not found in MGI -> Ensembl" %(not_found)
    return gene_to_go_terms


def main():
    pass

if __name__ == "__main__":
    main()
