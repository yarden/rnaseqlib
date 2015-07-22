import os
import sys
import csv
import collections
import time

import pandas
import pandas.rpy.common as com

import rnaseqlib
import rnaseqlib.stats.count_diffexp as count_diffexp
import rnaseqlib.stats.rpy2_utils as rpy2_utils

import numpy
import rpy2.robjects as robjects
robjects.r('''
  library(edgeR)
''')
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

base = importr("base")

def run_edgeR(counts_fname, sample1_name, sample2_name,
              dispersion=0.1,
              gene_id_col="GeneID",
              delimiter="\t"):
    """
    Run edgeR on counts dataframe. Emulates basic edgeR call:
    
    > y <- DGEList(counts=x,group=group)
    > y <- calcNormFactors(y)
    > y <- estimateCommonDisp(y)
#    > y <- estimateTagwiseDisp(y)
    > et <- exactTest(y)
    > topTags(et)
    """
    robjects.r('''
        library(edgeR)
    ''')
    base = importr("base")
    # Assume two groups (pairwise comparison)
    groups = robjects.r('''1:2''')
    read_delim = robjects.r["read.delim"]
#    counts_file_params = {"sep": delimiter,
#                          "row.names": gene_id_col}
#    counts = read_delim(counts_fname, **counts_file_params)
    counts = pandas.read_table(counts_fname, sep=delimiter)
    row_names = list(counts[gene_id_col].values)
    # Select only relevant samples
    counts = counts[[sample1_name, sample2_name]]
    #counts = rpy2_utils.df_to_r(counts)
    counts = com.convert_to_r_dataframe(counts)
    counts.rownames = row_names
    col_names = [curr_col for curr_col in counts.colnames]
    for curr_sample in [sample1_name, sample2_name]:
        if curr_sample not in col_names:
            raise Exception, "No %s in counts file." %(curr_sample)
    # can include lib.sizes here if needed
    params = {'group' : groups}
    y = robjects.r.DGEList(counts, **params)
    y = robjects.r.calcNormFactors(y)
    y = robjects.r.estimateCommonDisp(y)
    et = robjects.r.exactTest(y, dispersion=dispersion)
    tags = robjects.r.topTags(et)
    tags_df = tags[0]
    result = {"exactTest": com.convert_robj(et),
              "topTags": com.convert_robj(tags_df),
              "y": y}
    return result
    

def main():
    if len(sys.argv) < 1:
        raise Exception, "Need counts filename to be passed."
    counts_fname = sys.argv[1]
    if not os.path.isfile(counts_fname):
        raise Exception, "Cannot find file %s" %(counts_fname)
    output_fname = "./edger_results.txt"
    sample1_name = "Condition1"
    sample2_name = "Condition2"
    results = \
      run_edgeR(counts_fname, sample1_name, sample2_name, output_fname)
    print results["exactTest"]

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
