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
        self.kmers = kmers
        # Compute kmer len
        self.kmer_len = None


    def cluster_by_seq(self, method="sw"):
        """
        Cluster the motifs by sequence.
        """
        if method == "sw":
            # Smith-Waterman clustering
            self.cluster_by_sw()

            
    def cluster_by_sw(self):
        """
        Cluster sequences pairwise by Smith-Waterman alignment.
        """
        # Make matrix with ij entry corresponding
        # to alignment between sequence i and sequence j
        pass
                
