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

    def get_exon_intergenic_ratio(self):
        pass

    def get_exon_intron_ratio(self):
        pass

    def get_seq_cycle_profile(self, fastq_filename,
                              first_n_seqs=None):#sample):
        """
        Compute the average 'N' bases (unable to sequence)
        as a function of the position of the read.
        """
#        print "Computing sequence cycle profile for: %s" %(sample)
#        fastq_file = read_open_fastq(sample.reads_filename)
        fastq_file = fastq_utils.read_open_fastq(fastq_filename)
        fastq_entries = fastq_utils.read_fastq(fastq_file)
        # Mapping from position in read to number of Ns
        num_n_bases = defaultdict(int)
        # Mapping from position in read to total number of
        # reads in that position
        num_reads = defaultdict(int)
        num_entries = 0
        print "Computing sequence cycle profile for: %s" %(fastq_filename)
        if first_n_seqs != None:
            print "Looking at first %d sequences only" %(first_n_seqs)
        for entry in fastq_entries:
            if first_n_seqs != None:
                # Stop at requested number of entries if asked to
                if num_entries >= first_n_seqs:
                    break
            header1, seq, header2, qual = entry
            seq_len = len(seq)
            for n in range(seq_len):
                if seq[n] == "N":
                    # Record occurrences of N
                    num_n_bases[n] += 1
                num_reads[n] += 1
            num_entries += 1
        # Compute percentage of N along each position
        percent_n = []
        for base_pos in range(max(num_reads.keys())):
            curr_percent_n = float(num_n_bases[base_pos]) / num_reads[base_pos]
            percent_n.append(curr_percent_n)
        return percent_n
