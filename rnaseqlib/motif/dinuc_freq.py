##
## Wrapper for dinucleotide usage of pycogent
##
## Yarden Katz
##
import os
import sys
import time

from collections import defaultdict

import numpy as np

import cogent
from cogent.core.usage import DinucUsage

import rnaseqlib
import rnaseqlib.utils as utils

class DinucFreqs:
    """
    Dinucleotide frequencies. Wrapper for pycogent class.
    """
    def __init__(self, seq,
                 overlapping=True,
                 normalize=True):
        self.seq = seq.upper()
        self.overlapping = overlapping
        self.normalize = normalize
        # Sequence length
        self.len = len(seq)
        # Calculate dinuc. frequencies
        self.du = DinucUsage(seq, Overlapping=overlapping)
        if normalize:
            self.du.normalize()
        # Calculate frequencies for individual bases
        A_base = "A"
        T_base = "T"
        U_base = "U"
        G_base = "G"
        C_base = "C"
        self.bases = (A_base, T_base, G_base, C_base)
        self.base_freqs = defaultdict(int)
        for curr_base in self.bases:
            self.base_freqs[curr_base] = \
                (self.seq.count(curr_base) / float(self.len))
        # Equalize T/U -- pick the greater of the two frequencies
        # and then set them as equal
        self.base_freqs[T_base] = max((self.base_freqs["T"],
                                       self.base_freqs["U"]))
        self.base_freqs[U_base] = self.base_freqs[T_base]


    def get_dinuc_freqs_from(self, base, all_bases="ATGC"):
        """
        Get a list of dinucleotide base from the given base to
        each other base.

        base -> A
        base -> T
        ...
        """
        return [self.du["%s%s" %(base, possible_base)] \
                for possible_base in all_bases]

            
    def prob_score(self, subseq):
        """
        Score probability of subseq in sequence.
        """
        if len(subseq) == 0:
            return 0
        # Score first base
        total_logscore = np.log(self.base_freqs[subseq[0]])
        for prev_base, next_base in utils.iter_by_pair(subseq, 1):
            # Score current dinucleotide
            curr_dinuc = "%s%s" %(prev_base, next_base)
            # Divide by sum of all other transitions from the previous base
            curr_dinuc_freq = self.du[curr_dinuc]
            denom_dinuc_freqs = np.sum(self.get_dinuc_freqs_from(prev_base))
            total_logscore += (np.log(curr_dinuc_freq) - \
                               np.log(denom_dinuc_freqs))
        total_score = np.exp(total_logscore)
        return total_score


    def get_expected_num(self, subseq):
        """
        Calculcate the number of expected occurrences of subseq
        based on dinucleotide frequencies.
        """
        subseq_len = len(subseq)
        # Calculcate score of sequence (its probability)
        subseq_score = self.prob_score(subseq)
        # Calculate the number of possible positions for the subseq
        # to occur within the larger sequence
        num_positions = self.len - subseq_len + 1
        # Expected number is the score times the number of possible
        # positions
        exp_num = num_positions * subseq_score
        return exp_num


    def __str__(self):
        return "DinucFreqs(len=%d, seq=%s)" %(self.len,
                                              self.seq)


    def __repr__(self):
        return self.__str__()

    
        
