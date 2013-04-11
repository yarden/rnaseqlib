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
import rnaseqlib.motif.homer_utils as homer_utils
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
                 output_dir,
                 exp_coords_fname=None,
                 control_coords_fname=None,
                 genome="mm9"):
        # FASTAs representing the sequences for exp and control
        # conditions
        self.exp_seqs_fname = exp_seqs_fname
        self.control_seqs_fname = control_seqs_fname
        # Coordinates files for exp and control conditions
        self.exp_coords_fname = exp_coords_fname
        self.control_coords_fname = control_coords_fname
        # Kmer lengths to consider
        self.kmer_lens = kmer_lens
        self.output_dir = output_dir
        # Optional genome name
        self.genome = genome
        # Define Kmer tables for each kmer length
        self.exp_kmer_tables = []
        self.control_kmer_tables = []
        self.logger = utils.get_logger("MotifSet",
                                       os.path.join(self.output_dir, "logs"))
        

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
        all_df["exp_over_control"] = \
            all_df["counts_exp"] / all_df["counts_control"].apply(float)
        # Sort by fold enrichment
        all_df.sort(columns="exp_over_control",
                    inplace=True,
                    ascending=False)
        all_df.to_csv(os.path.join(self.output_dir, "enrichment.txt"),
                      sep="\t",
                      float_format="%.3f")
        print all_df["exp_over_control"]
        print all_df[0:100]


    def find_motifs_homer(self, output_dir,
                          homer_kmer_lens=[4,5,6,7,8]):
        """
        Find motifs with Homer.
        """
        output_dir = os.path.join(output_dir, "homer_output")
        utils.make_dir(output_dir)
        params = {"-rna": "",
                  "-len": ",".join(homer_kmer_lens)}
        # Run on exp
        homer_utils.run_homer(self.logger,
                              self.exp_coords_fname,
                              self.genome,
                              os.path.join(output_dir, "exp"),
                              params)
        # Run on control
        homer_utils.run_homer(self.logger,
                              self.control_coords_fname,
                              self.genome,
                              os.path.join(output_dir, "control"),
                              params)
        

    def output_enriched_kmers(self, num_shuffles=100):
        """
        Output kmers enriched in experimental set versus control.

        Find motifs representative of each set by finding
        motifs enriched in it compared to dinucleotide shuffles.
        """
        enriched_kmers = {}
        # Shuffled exp filename and control filename
        exp_shuffled = \
            kmer_utils.ShuffledFasta(self.exp_seqs_fname,
                                     os.path.join(self.output_dir, "exp_shuffled"))
        control_shuffled = \
            kmer_utils.ShuffledFasta(self.control_seqs_fname,
                                     os.path.join(self.output_dir, "control_shuffled"))
        for kmer_len in self.kmer_lens:
            print "Counting kmers of length %d" %(kmer_len)
            exp_dinuc_enriched_fname = \
                os.path.join(self.output_dir,
                             "exp_dinuc_enriched.%d_kmer.txt" %(kmer_len))
            exp_kmers = \
                kmer_utils.Kmers(kmer_len,
                                 fasta_fname=self.exp_seqs_fname,
                                 shuffled_fasta=exp_shuffled)
            # Count Kmers and output the ones that are enriched
            exp_dinuc_kmers = \
                exp_kmers.get_enriched_kmers(self.output_dir,
                                             num_shuffles=num_shuffles)
            exp_dinuc_results = \
                exp_kmers.output_enriched_kmers(exp_dinuc_kmers,
                                                exp_dinuc_enriched_fname)
            # Do same for control sequences
            control_dinuc_enriched_fname = \
                os.path.join(self.output_dir,
                             "control_dinuc_enriched.%d_kmer.txt" %(kmer_len))
            control_kmers = \
                kmer_utils.Kmers(kmer_len,
                                 fasta_fname=self.control_seqs_fname,
                                 shuffled_fasta=control_shuffled)
            control_dinuc_kmers = \
                control_kmers.get_enriched_kmers(self.output_dir,
                                                 num_shuffles=num_shuffles)
            control_dinuc_results = \
                control_kmers.output_enriched_kmers(control_dinuc_kmers,
                                                    control_dinuc_enriched_fname)
            # Remove extraneous columns
            cols = [c for c in exp_dinuc_results.columns \
                    if not c.startswith("shuffled_counts")]
            exp_dinuc_results = exp_dinuc_results[cols]
            control_dinuc_results = control_dinuc_results[cols]
            enriched_kmers[kmer_len] = {"exp": exp_dinuc_results,
                                        "control": control_dinuc_results}
        return enriched_kmers
        

