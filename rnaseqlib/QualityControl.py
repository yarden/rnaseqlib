import os
import sys
import time

import rnaseqlib
import rnaseqlib.fastq_utils as fastq_utils

from collections import defaultdict


class QualityControl:
    """ 
    Quality control object. Defined for every
    RNA-Seq library.
    """
    def __init__(self, settings_info):
        self.settings_info = settings_info

    def get_exon_intergenic_ratio():
        pass

    def get_exon_intron_ratio():
        pass

    def get_seq_cycle_profile(sample):
        """
        Compute the average 'N' bases (unable to sequence)
        as a function of the position of the read.
        """
        print "Computing sequence cycle profile for: %s" %(sample)
        fastq_file = read_open_fastq(sample.reads_filename)
        fastq_entries = read_fastq(fastq_file)
        # Mapping from position in read to number of Ns
        num_n_bases = defaultdict(int)
        # Mapping from position in read to total number of
        # reads in that position
        num_reads = defaultdict(int)
        for entry in fastq_entries:
            
