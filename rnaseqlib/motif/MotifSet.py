##
## MotifSet: represent a comparison between two samples
## for motifs or a sample with its dinucleotide shuffled self.
##
import os
import sys
import time
import operator

import rnaseqlib
import rnaseqlib.motif.kmer_utils as kmer_utils
import rnaseqlib.utils as utils


class MotifSet:
    """
    Representation of a motifs comparison between two sets
    of sequences (experimental set vs. control set.)
    """
    def __init__(self,
                 exp_seqs_fname,
                 control_seqs_fname,
                 kmer_lens,
                 output_dir):
        # Coordinates representing the experimental coordinates
        # (either a BED or a GFF file)
        self.exp_seqs_fname = exp_seqs_fname
        # Control coordinates filename (optional)
        self.control_seqs_fname = control_seqs_fname
        # Kmer lengths to consider
        self.kmer_lens = kmer_lens
        self.output_dir = output_dir
        # Define Kmer tables for each kmer length
        self.exp_kmer_tables = []
        self.control_kmer_tables = []
        

    def get_enriched_kmers(self, kmer_len, exp_counts, control_counts):
        """
        Get kmers enriched in the experimental sample compared to
        the control for a particular kmer length.
        """
        all_kmers = kmer_utils.enumerate_kmers(kmer_len)
        # Make kmers into dataframe
        exp_df = kmer_utils.counts_to_df(exp_counts, label="exp")
        control_df = kmer_utils.counts_to_df(control_counts, label="control")
        all_df = exp_df.join(control_df, how="outer")
        print "---*"*5
        print "--EXP",all_df["counts_exp"]
        print "--CTRL",all_df["counts_control"]
        all_df["exp_over_control"] = \
            all_df["counts_exp"] / all_df["counts_control"]
        # Sort by fold enrichment
        all_df.sort(columns="exp_over_control",
                    inplace=True,
                    ascending=False)
        all_df.to_csv(os.path.join(self.output_dir, "enrichment.txt"),
                      sep="\t",
                      float_format="%.3f")
        print all_df["exp_over_control"]
        print all_df[0:100]
        


    def output_enriched_kmers(self):
        for kmer_len in self.kmer_lens:
            print "Counting kmers of length %d" %(kmer_len)
            exp_kmers = \
                kmer_utils.Kmers(kmer_len,
                                 fasta_fname=self.exp_seqs_fname)
            control_kmers = \
                kmer_utils.Kmers(kmer_len,
                                 fasta_fname=self.control_seqs_fname)
            # Count kmers in exp and control sequences
            exp_counts = \
                exp_kmers.count_kmers(self.exp_seqs_fname,
                                      self.output_dir)
            control_counts = \
                control_kmers.count_kmers(self.control_seqs_fname,
                                          self.output_dir)
            # Sort the kmers by counts
            self.get_enriched_kmers(kmer_len, exp_counts, control_counts)
            # exp_hits = sorted(exp_counts.iteritems(),
            #                   key=operator.itemgetter(1),
            #                   reverse=True)
            # control_hits = sorted(control_counts.iteritems(),
            #                       key=operator.itemgetter(1),
            #                       reverse=True)
            # print "EXP HITS TOP 10"
            # print exp_hits[0:10]
            # print "CONTROL HITS TOP 10"
            # print control_hits[0:10]
            
        

