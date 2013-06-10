##
## Cogent-based sequence counter
##
import os
import sys
import time

import pandas

import numpy as np

import cogent
from cogent.core.usage import DinucUsage


import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.motif.dinuc_freq as dinuc_freq
import rnaseqlib.fasta_utils as fasta_utils

def overlap_count(s, sub):
    """
    Count overlapping occurrences of sub in s.
    """
    count = start = 0
    while True:
        start = s.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


class SeqCounter:
    """
    Sequence counter for a given FASTA file.
    """
    def __init__(self, fasta_fname):
        self.fasta_fname = fasta_fname
        self.seqs = fasta_utils.read_fasta(self.fasta_fname)


    def count(self, subseq):
        """
        Count occurrences of subseq across all sequences.
        """
        overlapping_counts = []
        for curr_seq in self.seqs:
            overlapping_counts.append(overlap_count(curr_seq[1], subseq))
        return overlapping_counts


    def obs_over_exp_counts(self, subseq):
        """
        Get observed over expected ratio of coutns (non-log!) of
        subseq in all sequences.
        """
        counts_and_ratio = []
        t1 = time.time()
        num_seqs = 0
        for curr_seq in self.seqs:
            seq_name = curr_seq[0]
            # Observed counts for subseq
            obs_count = overlap_count(curr_seq[1], subseq)
            # Expected counts for subseq based on
            # dinucleotide frequencies
            dinuc_freq_obj = dinuc_freq.DinucFreqs(curr_seq[1])
            exp_count = dinuc_freq_obj.get_expected_num(subseq)
            if exp_count == 0:
                # If expected count is 0, don't divide by it, just
                # consider the ratio 0
                ratio = 0
            else:
                ratio = float(obs_count) / float(exp_count)
            # Collect raw counts and ratio in order:
            # sequence name, obs count, exp count, obs / exp ratio
            counts_and_ratio.append([seq_name, obs_count, exp_count, ratio])
            num_seqs += 1
            #####
            if num_seqs == 100:
                print "QUITTING EARLY!"
                break
        t2 = time.time()
        print "Counting occurrences in %d sequences took %.2f seconds" \
              %(num_seqs, (t2 - t1))
        col_names = ["header", "obs_count", "exp_count", "ratio"]
        counts_and_ratio = \
            pandas.DataFrame(np.array(counts_and_ratio),
                             columns=col_names)
        return counts_and_ratio
            

    def expected_dinuc_count(self, subseq):
        """
        Get expected number of occurrences of subseq
        for each seq
        """
        pass


    def __repr__(self):
        return self.__str__()


    def __str__(self):
        return "SeqCounter(fasta=%s)" %(self.fasta_fname)
        



