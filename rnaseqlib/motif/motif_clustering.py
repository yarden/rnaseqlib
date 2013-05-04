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
import rnaseqlib.stats.stats_utils as stats_utils
import rnaseqlib.stats.clustering as clustering
import rnaseqlib.utils as utils


class MotifCluster:
    def __init__(self, kmers):
        """
        Store a list of kmers.
        """
        self.kmers = kmers
        # Compute kmer len
        self.kmer_len = None
        # Levenshtein distance function
        self.edit_dist_func = \
            lambda x, y: stats_utils.leven_dist(x[0], y[0])
        #self.edit_dist_func = \
        #    lambda x, y: stats_utils.leven_dist(x, y)


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


    def cluster_by_edit(self, kmers, linkage_method):
        """
        Cluster sequences by edit distances.

        Parameters:
        -----------
        kmers : flat list of kmers
        linkage_method : determines linkage function for hierarchical
        clustering ('average', 'single', ...).

        """
        # Nest kmers to create matrix for clustering purposes.
        kmers = [kmers]
        data = np.array(kmers)
        hclust = clustering.hierarchical_clust(data,
                                               self.edit_dist_func,
                                               linkage_method)
        return hclust


def output_global_alignment(kmers_fname, output_dir):
    """
    Output a global alignment (*.aln) for a set of kmers.

    Using clustawl for now.

    Parameters:
    -----------
    kmers_fname : filename of FASTA file containing kmers
    output_dir : output directory
    """
    utils.make_dir(output_dir)
    output_fname = \
        os.path.join(output_dir,
                     "%s.aln" %(os.path.basename(kmers_fname)))
    if os.path.isfile(output_fname):
        print "Alignment filename %s exists. Skipping..." \
              %(output_fname)
    clustalw_cmd = \
        "clustalw -INFILE=%s -OUTFILE=%s -PIM" %(kmers_fname,
                                                 output_fname)
    print "Executing: %s" %(clustalw_cmd)
    t1 = time.time()
    os.system(clustalw_cmd)
    t2 = time.time()
    print "Global alignment took %.2f minutes." %((t2 - t1)/60.)
    return output_fname
    

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


if __name__ == "__main__":
    main()
