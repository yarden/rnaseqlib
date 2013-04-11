##
## Class for representing Bind-n-Seq odds ratio/counts data
##
import os
import sys
import time
import glob

import pandas
import numpy as np
import scipy

import rnaseqlib
import rnaseqlib.motif.meme_utils as meme_utils
import rnaseqlib.utils as utils


def load_bindnseq_or_file(odds_ratio_fname, skiprows=1):
    """
    Load up bindnseq odds ratio file as a pandas DataFrame.
    """
    odds_ratio_df = \
        pandas.read_table(odds_ratio_fname,
                          sep="\t",
                          skiprows=skiprows)
    odds_ratio_df = odds_ratio_df.rename(columns={"#": "kmer"})
    return odds_ratio_df


def load_bindnseq_counts_file(counts_fname, skiprows=1):
    """
    Load up bindnseq counts file as a pandas DataFrame.
    """
    counts_df = None
    with open(counts_fname) as counts_in:
        header = counts_in.readline().strip()
        if header.startswith("#"):
            header = header[1:]
        columns = header.split("\t")
        counts_df = pandas.read_table(counts_fname,
                                      skiprows=skiprows)
        counts_df.columns = columns
    return counts_df
        

class BindnSeq:
    def __init__(self, results_dir, output_dir,
                 label=None):
        """
        Load up results directory
        """
        self.output_dir = output_dir
        self.logger_label = label
        if self.logger_label is None:
            self.logger_label = "BindnSeq"
        self.logger = utils.get_logger(self.logger_label,
                                       self.output_dir)
        self.results_dir = results_dir
        self.label = label
        # All kmer lengths to load
        self.kmer_lens = [4, 5, 6, 7, 8, 9]
        # Odds ratios (DataFrames indexed by kmer length)
        self.odds_ratios = {}
        # Counts (DataFrames indexed by kmer length)
        self.counts = {}


    def load_results(self, results_dir):
        # Load the odds ratios
        self.load_odds_ratios(results_dir)
        # Don't load the counts for now
        #self.load_counts()


    def load_odds_ratios(self, results_dir):
        """
        Load odds ratios.
        """
        self.logger.info("Loading BindnSeq results from: %s" %(results_dir))
        odds_ratio_fnames = glob.glob(os.path.join(results_dir, "*mer_OR"))
        num_or_files = len(odds_ratio_fnames)
        if num_or_files == 0:
            self.logger.critical("No OR files to load!")
            sys.exit(1)
        # Load the odds ratio files
        for or_fname in odds_ratio_fnames:
            kmer_len = int(os.path.basename(or_fname).split("mer")[0])
            odds_ratio_df = load_bindnseq_or_file(or_fname)
            self.odds_ratios[kmer_len] = odds_ratio_df
        self.logger.info("  - Found %d OR files" %(num_or_files))
        return self.odds_ratios


    def get_conc_cols(self, df, suffix='nM'):
        """
        Return concentration columns of a DataFrame. Assume concentration
        columns end in suffix.
        """
        conc_cols = [col for col in df.columns \
                     if col.endswith(suffix)]
        return conc_cols


    def rank_enriched_kmers(self, or_df,
                            rank_col="rank",
                            method="max"):
        """
        Given a DataFrame with OR values, return a new DataFrame
        with kmers sorted according to enrichment. The enrichment
        value will be in a 'rank' column and is computed based on
        the given 'method' argument:

        - method='max': use maximum enrichment across all concentrations
          as enrichment value of a kmer
        - method='mean': use average enrichment across all concentrations
          as enrichment value of a kmer
        """
        ranked_or_df = None
        conc_cols = self.get_conc_cols(or_df)
        ranked_or_df = or_df.copy()
        if method == "max":
            # maximum enrichment method
            # rank the kmers by maximum enrichment across concentrations
            ranked_or_df[rank_col] = or_df[conc_cols].max(axis=1)
        elif method == "mean":
            # mean enrichment method
            # rank the kmers by average enrichment across concentrations
            ranked_or_df[rank_col] = or_df[conc_cols].mean(axis=1)
        else:
            self.logger.criticla("Unknown enrichment method %s" %(method))
            sys.exit(1)
        # Sort resulting DataFrame by rank in descending order
        ranked_or_df.sort(columns=[rank_col],
                          inplace=True,
                          ascending=False)
        return ranked_or_df


    def run_meme_on_enriched_kmers(self, output_dir,
                                   fold_enriched_cutoff=2,
                                   method="max",
                                   len_to_output=None):
        """
        Run MEME on all enriched kmers.
        """
        self.logger.info("Running MEME on enriched BindnSeq kmers...")
        self.logger.info("  - Output dir: %s" %(output_dir))
        self.logger.info("  - Fold enrichment cutoff: %.1f" %(fold_enriched_cutoff))
        self.logger.info("  - Enrichment method: %s" %(method))
        # Make directory for all the kmer sequences to be
        # processed by MEME
        self.seqs_dir = os.path.join(output_dir, "seqs")
        utils.make_dir(self.seqs_dir)
        # Output all enriched kmers to file
        if len_to_output is None:
            len_to_output = "all"
        self.seqs_fname = \
            os.path.join(self.seqs_dir,
                         "enriched_kmers.cutoff_%.1f.method_%s.%s_kmers.fasta" \
                         %(fold_enriched_cutoff, method, str(len_to_output)))
        self.logger.info("Outputting sequences as FASTA to: %s" %(self.seqs_fname))
        seqs_out = open(self.seqs_fname, "w")
        for kmer_len in [4,5,6]:#self.kmer_lens:
            if len_to_output != "all":
                if len_to_output != kmer_len:
                    print "Skipping %d" %(kmer_len)
                    continue
            odds_ratios = self.odds_ratios[kmer_len]
            # Rank the odds ratios
            ranked_ratios = self.rank_enriched_kmers(odds_ratios)
            # Select only the kmers that meet the cutoff
            enriched_ratios = \
                ranked_ratios[ranked_ratios["rank"] >= fold_enriched_cutoff]
            # Write those to file
            for kmer in enriched_ratios["kmer"].values:
                header = ">%s\n" %(kmer)
                seq = "%s\n" %(kmer)
                seqs_out.write(header)
                seqs_out.write(seq)
        seqs_out.close()
        # Run MEME on FASTA file with kmers
        output_dir = os.path.join(output_dir, "meme_output")
        utils.make_dir(output_dir)
        self.logger.info("Running MEME on enriched BindnSeq kmers...")
        self.logger.info("  - MEME output dir: %s" %(output_dir))
        if len(glob.glob(os.path.join(output_dir, "*"))) >= 1:
            self.logger.info("MEME output exists. Skipping...")
            return
        meme_utils.run_meme(self.logger, self.seqs_fname, output_dir)
        

    def __str__(self):
        return "BindnSeq(input_dir=%s)" %(self.results_dir)


    def __repr__(self):
        return self.__str__()
        
        
