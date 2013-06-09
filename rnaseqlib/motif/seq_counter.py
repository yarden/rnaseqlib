##
## Cogent-based sequence counter
##
import os
import sys
import time

import cogent
from cogent.core.usage import DinucUsage

import rnaseqlib
import rnaseqlib.utils as utils


class SeqCounter:
    """
    Sequence counter.
    """
    def __init__(self, seqs):
        self.seqs = seqs


def get_dinuc_freqs(seq, normed=True, overlapping=True):
    """
    Return dinucleotide frequencies
    """
    du = DinucUsage(seq, Overlapping=overlapping)
    if normed:
        du.normalize()
    return du


def score_seq_by_dinuc(seq, du):
    """
    Score probability of 'seq' given du
    """

