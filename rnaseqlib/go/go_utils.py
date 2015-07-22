##
## GO related utilities
##
import os
import sys
import time

import rnaseqlib

def get_genes_to_go_terms(gene_assoc_fname,
                          gene_id_col=10,
                          go_term_col=4,
                          gene_id_prefix="ENS"):
    """
    Return a mapping (dictionary) from gene IDs to their GO terms.

    E.g. map ENSMUSG... -> GO:XXXX

    Args:
    - gene_assoc_fname: gene association filename from GO (*.mgi file)
    - gene_id_col: the column number (0-based) in gene_assoc_fname
      that contains the gene ID.
    """
    gene_to_go_terms = {}
    with open(gene_assoc_fname) as go_in:
        for line in go_in:
            if line.startswith("!"):
                # Skip GO comments
                continue
            fields = line.split("\t")
            if len(fields) < gene_id_col:
                continue
            gene_id = fields[gene_id_col]
            if not gene_id.startswith(gene_id_prefix):
                # Skip gene entries that do not start with required ID
                continue
            gene_to_go_terms[gene_id] = fields[go_term_col]
    return gene_to_go_terms


def main():
    pass

if __name__ == "__main__":
    main()
