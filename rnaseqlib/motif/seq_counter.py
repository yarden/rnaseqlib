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

        
def overlap_count_with_starts(s, sub):
    """
    Count overlapping occurrences of sub in s.
    Return numbers of counts as well as vector of starting
    positions (0-based) of each.
    """
    count = start = 0
    start_positions = []
    while True:
        start = s.find(sub, start) + 1
        if start > 0:
            # Important to subtract 1 to get 0-based
            # indexing
            start_positions.append(start - 1) 
            count += 1
        else:
            return count, start_positions

        
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


    def count_dinuc_subseqs(self, seq, subseqs):
        """
        Given a sequence and a list of subsequences,
        return the observed number of occurrences of all subseqs
        in curr_seq,        
        """
        dinuc_freq_obj = dinuc_freq.DinucFreqs(seq)
        obs_counts = []
        exp_counts = []
        for subseq in subseqs:
            # Observed count of subsequence
            obs_count = overlap_count(seq, subseq)
            obs_counts.append(float(obs_count))
            # Expected count of subsequence
            exp_count = dinuc_freq_obj.get_expected_num(subseq)
            exp_counts.append(exp_count)
        return np.array(obs_counts), np.array(exp_counts)


    def count_subseqs(self, seq, subseqs):
        """
        Count occurrences of subseqs in seq.
        """
        obs_counts = \
            [float(overlap_count(seq, subseq)) \
             for subseq in subseqs]
        return np.array(obs_counts)


    def count_subseqs_with_starts(self, seq, subseqs):
        """
        Count occurrences of subseqs in seq. For each subseq,
        return a list of its starting positions in seq.
        """
        start_positions = []
        for subseq in subseqs:
            counts, starts = overlap_count_with_starts(seq, subseq)
            start_positions.append(starts)
        return start_positions


    def get_subseq_densities(self, subseqs):
        """
        Collect densities of subsequences in sequence.
        """
        entries = []
        KB = float(1000)
        for curr_seq in self.seqs:
            # Sequence name is FASTA header minus '>' prefix
            seq_name = curr_seq[0][1:]
            seq = curr_seq[1]
            seq_len = len(seq)
            # Counts of occurrences for all subsequences
            subseq_counts = self.count_subseqs(seq, subseqs)
            subseq_len = len(subseqs[0])
            # Normalize by the number of possible counts
            # per kb 
            density_denom = \
                (float(len(seq) - subseq_len + 1)) / KB
            # Vector of densities, one for each subsequence
            densities = subseq_counts / density_denom
            #mean_density = np.mean(densities)
            #median_density = np.median(densities)
            sum_density = sum(densities)
            max_density = max(densities)
            obs_counts_str = \
                ",".join(["%d" %(int(c)) for c in subseq_counts])
            densities_str = \
                ",".join([str(d) for d in densities])
            seq_len_in_kb = seq_len / KB
            entries.append([seq_name,
                            sum_density,
                            max_density,
                            obs_counts_str,
                            densities_str,
                            seq_len,
                            seq_len_in_kb])
        col_names = ["header",
                     "sum_density", "max_density",
                     "obs_counts", "densities", "seq_len",
                     "seq_len_in_kb"]
        entries = pandas.DataFrame(entries,
                                   columns=col_names).set_index("header")
        # Sort in place by mean density in descending order
        entries.sort(column=["sum_density", "max_density"],
                     ascending=False,
                     inplace=True)
        return entries

    
    def obs_over_exp_counts_dinuc(self, subseqs):
        """
        Get observed over expected ratio of counts (non-log!) of
        subsequences in all sequences.
        """
        entries = []
        t1 = time.time()
        num_seqs = 0
        for curr_seq in self.seqs:
            # Sequence name is FASTA header without leading '>'
            seq_name = curr_seq[0][1:]
            # Get observed and expected counts for all subseqs
            # in the current sequence
            obs_counts, exp_counts = self.count_subseqs(curr_seq[1], subseqs)
            # Calculate ratios
            ratios = obs_counts / exp_counts
            # All ratios
            ratios_str = ",".join(["%.3f" %(r) for r in ratios])
            # The maximum ratio
            max_ratio_indx, max_ratio = utils.max_item(ratios)
            # Get the observed counts of the kmer with highest ratio
            max_ratio_obs_count = int(obs_counts[max_ratio_indx])
            obs_counts_str = ",".join(["%d" %(int(oc)) for oc in obs_counts])
            exp_counts_str = ",".join(["%.2f" %(ec) for ec in exp_counts])
            # Collect raw counts and ratio in order:
            # sequence name, obs counts, exp counts, obs / exp ratios
            entries.append([seq_name,
                            max_ratio,
                            max_ratio_obs_count,
                            obs_counts_str,
                            exp_counts_str,
                            ratios_str])
            num_seqs += 1
            if num_seqs == 100:
                print "Quitting early!"
                print "=" * 10
                break
        t2 = time.time()
        print "Counting occurrences in %d sequences took %.2f seconds" \
              %(num_seqs, (t2 - t1))
        col_names = ["header", "max_ratio", "max_ratio_obs_count",
                     "obs_counts", "exp_counts", "ratios"]
        entries = \
            pandas.DataFrame(np.array(entries),
                             columns=col_names).set_index("header")
        # Sort in descending order
        entries.sort(column=["max_ratio"],
                     ascending=False,
                     inplace=True)
        return entries
            

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
        



