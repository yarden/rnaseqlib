import os
import sys
import csv
import collections
import time

import pandas

import rnaseqlib
import rnaseqlib.stats.count_diffexp as count_diffexp

import numpy
import rpy2.robjects as robjects
robjects.r('''
  library(edgeR)
''')
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
base = importr("base")

def run_edgeR(counts_fname, sample1_name, sample2_name, output_fname,
              gene_id_col="GeneID",
              delimiter="\t"):
    """
    Run edgeR on counts dataframe.
    """
    robjects.r('''
        library(edgeR)
    ''')
    base = importr("base")
    # Assume two groups (pairwise comparison)
    groups = robjects.r('''1:2''')
    read_delim = robjects.r["read.delim"]
    counts_file_params = {"sep": delimiter,
                          "row.names": gene_id_col}
    counts = read_delim(counts_fname, **counts_file_params)
    print "counts: ", counts
    # can include lib.sizes here if needed
    params = {'group' : groups}
    dgelist = robjects.r.DGEList(counts, **params)
    print "DGELIST: ", dgelist
    # # perform Poisson adjustment and assignment as recommended in the manual
    # robjects.globalenv['dP'] = dgelist
    # print "testing..."
    # robjects.r('''
    #   msP <- de4DGE(dP, doPoisson = TRUE)
    #   dP$pseudo.alt <- msP$pseudo
    #   dP$common.dispersion <- 1e-06
    #   dP$conc <- msP$conc
    #   dP$common.lib.size <- msP$M
    # ''')
    # print "are we here?"
    # #dgelist = robjects.globalEnv['dP']
    # dgelist = robjects.globalenv['dP']
    # de = robjects.r.exactTest(dgelist)
    # tags = robjects.r.topTags(de, n=len(genes))
    # tag_table = tags[0]
    # indexes = [int(t) - 1 for t in tag_table.rownames()]
    # # can retrieve either raw or adjusted p-values
    # #pvals = list(tags.r['p.value'][0])
    # pvals = list(tag_table.r['adj.p.val'][0])
    

def main():
    if len(sys.argv) < 1:
        raise Exception, "Need counts filename to be passed."
    counts_fname = sys.argv[1]
    if not os.path.isfile(counts_fname):
        raise Exception, "Cannot find file %s" %(counts_fname)
    output_fname = "./edger_results.txt"
    sample1_name = "Condition 1"
    sample2_name = "Condition 2"
    run_edgeR(counts_fname, sample1_name, sample2_name, output_fname)
    

if __name__ == "__main__":
    main()
    

# def run_edgeR(counts_fname, sample1_name, sample2_name, output_fname,
#               gene_id_col="GeneID"):
#     """
#     Run EdgeR on counts dataframe.

#     Compares 'sample1_name' to 'sample2_name' and outputs the results
#     to output_fname.
#     """
#     print "Running edgeR (comparing %s vs. %s)" %(sample1_name,
#                                                   sample2_name)
#     print "  - Output file: %s" %(output_fname)
#     counts = \
#       count_diffexp.read_count_file(counts_fname,
#                                     conditions_to_load=[sample1_name,
#                                                         sample2_name])
#     # Get only subset of data that is relevant 
#     t1 = time.time()
#     data, groups, sizes, conditions, genes = \
#       count_diffexp.edger_matrices(counts)
#     probs = count_diffexp.run_edger(data, groups, sizes, genes)
#     with open(output_fname, "w") as outfile:
#         count_diffexp.write_outfile(outfile, genes, conditions, counts, probs)
#     t2 = time.time()
#     print "EdgeR took %.2f seconds" %(t2 - t1)
#     in_file.close()
