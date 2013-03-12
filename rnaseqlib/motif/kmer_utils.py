##
## Utilities for dealing with kmers
##
import os
import sys
import time

import numpy as np
import pandas

import operator

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils
import rnaseqlib.fasta_utils as fasta_utils
import rnaseqlib.motif.altschulEriksonDinuclShuffle as dinuc_shuffle
import rnaseqlib.motif.jellyfish_utils as jf_utils

from collections import defaultdict

import khmer


class Kmers:
    """
    Representation of a set of kmers a FASTA file.
    """
    def __init__(self, kmer_len, fasta_fname=None):
        self.kmer_len = kmer_len
        # FASTA file attached to kmers
        self.fasta_fname = fasta_fname
        # Shuffled FASTA filenames
        self.shuffled_fasta_fnames = []
        # Mapping from kmers to their counts
        # in the FASTA file
        self.kmer_counts = defaultdict(int)
        # Mapping from kmers to their counts for each
        # dinucleotide shuffled version of the FASTA file
        self.shuffled_kmer_counts = []
        # Get the set of kmers
        self.kmers = enumerate_kmers(kmer_len)
        # Counts file for jf
        self.jf_counts_fname = None


    def get_dinuc_shuffled_fasta(self, output_dir, num_shuffles=100):
        """
        Get dinucleotide shuffled versions of the FASTA file.

        Output FASTA files to output directory.
        """
        utils.make_dir(output_dir)
        print "Shuffling FASTA %d times into: %s" %(num_shuffles,
                                                    output_dir)
        t1 = time.time()
        shuffled_fnames = []
        for shuffle_num in range(num_shuffles):
            shuffled_basename = os.path.basename(self.fasta_fname)
            # Remove FASTA extension
            shuffled_basename = shuffled_basename.rsplit(".", 1)[0]
            # Record that it's a shuffle in the filename
            shuffled_basename = "%s.shuffle_%d.fa" %(shuffled_basename,
                                                     shuffle_num)
            shuffled_fname = os.path.join(output_dir, shuffled_basename)
            if not os.path.isfile(shuffled_fname):
                output_dinuc_shuffled_fasta(self.fasta_fname,
                                            shuffled_fname)
            shuffled_fnames.append(shuffled_fname)
        t2 = time.time()
        print "Shuffling took %.2f seconds" %(t2 - t1)
        self.shuffled_fasta_fnames = shuffled_fnames
        return self.shuffled_fasta_fnames


    def output_enriched_kmers(self, result_counts, output_fname):
        kmer_entries = []
        shuffled_counts = result_counts["shuffled_counts"]
        for kmer in self.kmers:
            # kmer name and its counts in the target FASTA file
            kmer_entry = {"kmer": kmer,
                          "counts": result_counts["counts"][kmer]}
            # Include counts from each of the shuffled runs
            shuffled_cols = []
            shuffled_vals = []
            for shuffle_num in range(len(shuffled_counts)):
                shuffled_count = shuffled_counts[shuffle_num][kmer]
                # Record as one based
                shuffled_col = "shuffled_counts_%d" %(shuffle_num + 1)
                shuffled_cols.append(shuffled_col)
                shuffled_vals.append(shuffled_count)
                kmer_entry.update({shuffled_col: shuffled_count})
            # Take mean of shuffled values
            kmer_entry.update({"shuffled_mean": np.mean(shuffled_vals)})
            # Take ratio of target count to mean of shuffled values
            target_to_shuffled_ratio = \
                kmer_entry["counts"] / kmer_entry["shuffled_mean"]
            kmer_entry.update({"target_to_shuffled":
                               target_to_shuffled_ratio})
            kmer_entries.append(kmer_entry)
        kmer_df = pandas.DataFrame(kmer_entries)
        print "Outputting kmer counts to: %s" %(output_fname)
        column_order = \
            ["kmer", "counts", "shuffled_mean", "target_to_shuffled"] + \
            shuffled_cols
        # Sort kmers by target_to_shuffled ratio
        kmer_df.sort(columns=["target_to_shuffled"],
                     ascending=False,
                     inplace=True)
        kmer_df.to_csv(output_fname,
                       sep="\t",
                       index=False,
                       cols=column_order)
            

    def get_enriched_kmers(self, output_dir, num_shuffles=1):
        """
        Get kmers that are enriched in FASTA relative to its nucleotide
        shuffled version.
        """
        # Create directory for shuffled FASTA files
        shuffled_output_dir = os.path.join(output_dir, "shuffled_fasta")
        utils.make_dir(shuffled_output_dir)
        # Make the shuffled FASTA files
        shuffled_fasta_fnames = \
            self.get_dinuc_shuffled_fasta(shuffled_output_dir,
                                          num_shuffles=num_shuffles)
        # First count the kmers in the target FASTA file
        kmer_counts = self.count_kmers(self.fasta_fname, output_dir)
        all_shuffled_counts = []
        for shuffled_fname in shuffled_fasta_fnames:
            shuffled_counts = self.count_kmers(shuffled_fname,
                                               output_dir)
            all_shuffled_counts.append(shuffled_counts)
        result_counts = {"counts": kmer_counts,
                         "shuffled_counts": all_shuffled_counts}
        return result_counts
        

    def count_kmers_dinuc_shuffled(self, output_dir):
        """
        Count Kmers in dinucleotide shuffled FASTA files.
        """
        self.shuffled_kmer_counts = []
        if len(self.shuffled_fasta_fnames):
            raise Exception, "No shuffled FASTA files to count kmers in."
        for shuffled_fname in self.shuffled_fasta_fnames:
            # Count kmers in each dinucleotide shuffled
            # FASTA file and keep these as a list
            # of Kmer counts
            kmer_counts = self.count_kmers(shuffled_fname, output_dir)
            shuffled_kmer_counts.append(kmer_counts)
        return shuffled_kmer_counts


    def count_kmers(self, fasta_fname, output_dir, method="jf"):
        # Load up fasta file
        kmer_counts = None
        if method == "jf":
            self.jf_counts_fname = self.count_kmers_jellyfish(fasta_fname,
                                                              output_dir)
            # Parse jf results
            kmer_counts = \
                jf_utils.jf_counts_to_dict(self.jf_counts_fname)
        else:
            raise Exception, "Do not support %s" %(method)
        return kmer_counts


    def output_kmer_counts(self, kmer_counts, counts_fname):
        """
        Output a set of kmer counts to a file.
        """
        pass

            
    def get_counts(self, kmer, as_df=False):
        """
        Get counts for kmer. 
        """
        if len(kmer) != self.kmer_len:
            print "WARNING: %s is not of indexed length %d" \
                  %(kmer, self.kmer_len)
            return 0
        counts = self.kmer_counts[kmer]
        return counts


    def get_sorted_kmers_by_count(self, reverse=True):
        """
        Get kmers in sorted ordering accoding to counts
        """
        sorted_counts = sorted(self.kmer_counts.iteritems(),
                               key=operator.itemgetter(1),
                               reverse=reverse)
        return sorted_counts
    
        
    ##
    ## Methods for counting kmers.
    ## 
    def count_kmers_naive(self):
        pass


    def count_kmers_jellyfish(self, fasta_fname, output_dir):
        """
        Count kmers with jellyfish.
        """
        jf_counts_fname = jf_utils.jf_count_kmers(fasta_fname,
                                                  self.kmer_len,
                                                  output_dir)
        return jf_counts_fname
        

    def count_kmers_khmer(self):
        """
        Count kmers with khmer.
        """
        # Reset counts
        self.kmer_counts = defaultdict(int)
        # Load FASTA file into ktable
        ktable = load_fastx_into_ktable(self.fasta_fname, self.kmer_len)
        for n in range(ktable.n_entries()):
            # Current kmer
            kmer = ktable.reverse_hash(n)
            # Kmer counts
            counts = ktable.get(n)
            self.kmer_counts[kmer] = counts


def counts_to_df(counts, label=None):
    """
    Convert a dictionary of counts to a pandas DataFrame.
    If label is given, make the 'counts' field into 'label_counts'.
    """
    entries = []
    for kmer in counts:
        entry = {"kmer": kmer}
        if label is not None:
            entry["counts"] = counts[kmer]
        else:
            entry["counts_%s" %(label)] = counts[kmer]
        entries.append(entry)
    df = pandas.DataFrame(df)
    return df
        


def enumerate_kmers(kmer_len):
    """
    Return all kmers as strings.
    """
    ktable = khmer.new_ktable(kmer_len)
    kmers = [ktable.reverse_hash(n) \
             for n in range(ktable.n_entries())]
    return kmers
    

def get_dinuc_shuffles(seq, num_shuffles=1):
    """
    Given a sequence, return num_shuffles-many
    dinucleotide matched shuffles of it.
    """
    shuffles = [dinuc_shuffle.dinuclShuffle(seq, DNA=True) \
                for n in range(num_shuffles)]
    return shuffles


def output_dinuc_shuffled_fasta(fasta_fname, shuffled_fasta_fname,
                                num_shuffles=1):
    """
    Given a FASTA file, output a dinucleotide shuffled version of it.
    """
    fasta_out = fastx_utils.write_open_fastx(shuffled_fasta_fname)
    for fastx_entry in fastx_utils.get_fastx_entries(fasta_fname):
        fastx_name, fastx_seq = fastx_entry
        shuffled_recs = []
        for shuffle_num in range(num_shuffles):
            shuffled_seq = \
                get_dinuc_shuffles(fastx_seq)[0]
            shuffled_rec = (fastx_name, shuffled_seq)
            shuffled_recs.append(shuffled_rec)
        fasta_utils.write_fasta(fasta_out, shuffled_recs)
    fasta_out.close()
        

def get_enriched_kmers_rel_to_bg(fasta_fname,
                                 control_fasta_fname):
    """
    Get enriched kmers relative to background set.
    """
    pass


def load_fastx_into_ktable(fastx_fname, kmer_len):
    t1 = time.time()
    ktable = khmer.new_ktable(kmer_len)
    # Load up the FASTA into ktable
    for fastx_entry in fastx_utils.get_fastx_entries(fastx_fname):
        fastx_name, fastx_seq = fastx_entry
        # Skip very short sequences
        if len(fastx_seq) < kmer_len:
            continue
        ktable.consume(fastx_seq)
    t2 = time.time()
    print "Loading up of seqs into ktable took %.2f seconds." %(t2 - t1)
    return ktable


def output_dinuc_enriched_kmers(logger,
                                fasta_fname,
                                output_dir,
                                kmer_lens,
                                num_shuffles=100):
    """
    Output enriched kmers in a FASTA file relative to
    a dinucleotide shuffled version of it.
    """
    logger.info("Output dinucleotide enriched Kmers..")
    logger.info("  - Input FASTA: %s" %(fasta_fname))
    logger.info("  - Output dir: %s" %(output_dir))
    utils.make_dir(output_dir)
    for kmer_len in kmer_lens:
        kmers = Kmers(kmer_len, fasta_fname=fasta_fname)
        # Get the enriched kmers
        results = kmers.get_enriched_kmers(output_dir,
                                           num_shuffles=num_shuffles)
        output_basename = \
            "%s.%d_kmers.counts" %(os.path.basename(fasta_fname),
                                   kmer_len)
        enrichment_fname = os.path.join(output_dir, output_basename)
        # Output enrichment result
        logger.info("Outputting enriched Kmers to: %s" %(enrichment_fname))
        kmers.output_enriched_kmers(results, enrichment_fname)


if __name__ == "__main__":
    pass
    #kmers = Kmers(8, fasta_fname=f)
    # Get enriched kmers
    #results = kmers.get_enriched_kmers("./", num_shuffles=20)
    #kmers.output_enriched_kmers(results, "./results.counts")

    
        


