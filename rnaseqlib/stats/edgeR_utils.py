import os
import sys
import csv
import collections
import time

import rnaseqlib
import rnaseqlib.stats.count_diffexp as count_diffexp

def run_edgeR(counts_fname, sample1_name, sample2_name, output_fname):
    """
    Run EdgeR on counts dataframe.

    Compares 'sample1_name' to 'sample2_name' and outputs the results
    to output_fname.
    """
    print "Running edgeR (comparing %s vs. %s)" %(sample1_name,
                                                  sample2_name)
    print "  - Output file: %s" %(output_fname)
    t1 = time.time()
    data, groups, sizes, conditions, genes = \
      count_diffexp.edger_matrices(counts)
    probs = count_diffexp.run_edger(data, groups, sizes, genes)
    with open(output_fname, "w") as outfile:
        count_diffexp.write_outfile(outfile, genes, conditions, counts, probs)
    t2 = time.time()
    print "EdgeR took %.2f seconds" %(t2 - t1)
