##
## Utilities for clustering motifs
##
import os
import sys
import time

import pandas

import numpy as np

import cogent
from cogent.align.algorithm import sw_align

import rnaseqlib
import rnaseqlib.utils as utils


class MotifCluster:
    def __init__(self, kmers):
        """
        Store a list of kmers.
        """
        self.kmers = kmers
        # Compute kmer len
        self.kmer_len = None


    def cluster_by_seq(self, method="sw"):
        """
        Cluster the motifs by sequence.
        """
        result = None
        if method == "sw":
            # Smith-Waterman clustering
            result = self.cluster_by_sw()
        return result

            
    def cluster_by_sw(self):
        """
        Cluster sequences pairwise by Smith-Waterman alignment.

        Returns distance matrix.
        """
        # Make pdist matrix with ij entry corresponding
        # to alignment between sequence i and sequence j
        score_matrix = []
        for kmer_i in self.kmers:
            score_row = []
            for kmer_j in self.kmers:
                alignment = sw_align(kmer_i, kmer_j,
                                     return_score=True)
                sw_score = alignment[1]
                score_row.append(sw_score)
            score_matrix.append(score_row)
        score_matrix = np.array(score_matrix)
        return score_matrix


def main():
    kmers = ["TGTAT", "CGTAT", "TTAGT", "TCTAT", "TCTAC"]
    motif_clust = MotifCluster(kmers)
    #result = motif_clust.cluster_by_seq()
    #print "Result: ", result
    data = np.transpose(np.array(kmers))
    print "DATA: "
    print data
    dist_func = lambda x, y: sw_align(x, y, return_score=True)
    linkage_method = "average"
    hclust = clustering.hierarchical_clut(np.array(kmers),
                                          dist_func,
                                          linkage_method)
    print "hclust:", hclust


if __name__ == "__main__":
    main()
