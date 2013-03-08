##
## Utilities for dealing with kmers
##
import os
import sys
import time

import operator

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils
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
        # Mapping from kmers to fasta files to their counts
        self.kmer_counts = defaultdict(int)
        # Get the set of kmers
        self.kmers = enumerate_kmers(kmer_len)
        # Counts file for jf
        self.jf_counts_fname = None


    def count_kmers(self, output_dir, method="jf"):
        if self.fasta_fname is None:
            return
        # Load up fasta file
        if method == "jf":
            self.jf_counts_fname = self.count_kmers_jellyfish(output_dir)
            # Parse jf results
            self.kmer_counts = \
                jf_utils.jf_counts_to_dict(self.jf_counts_fname)
        else:
            raise Exception, "Do not support %s" %(method)
        
            
    def get_counts(self, kmer):
        """
        Get counts for kmer
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


    def count_kmers_jellyfish(self, output_dir):
        """
        Count kmers with jellyfish.
        """
        jf_counts_fname = jf_utils.jf_count_kmers(self.fasta_fname,
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



def enumerate_kmers(kmer_len):
    """
    Return all kmers as strings.
    """
    ktable = khmer.new_ktable(kmer_len)
    kmers = [ktable.reverse_hash(n) \
             for n in range(ktable.n_entries())]
    return kmers
    

def get_dinuc_shuffles(seq, num_shuffles):
    """
    Given a sequence, return num_shuffles-many
    dinucleotide matched shuffles of it.
    """
    shuffles = [dinuc_shuffle.dinucShuffle(seq) \
                for n in range(num_shuffles)]
    return shuffles


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
        if len(fastx_seq) < kmer_len: continue
        ktable.consume(fastx_seq)
    t2 = time.time()
    print "Loading up of seqs into ktable took %.2f seconds" %(t2 - t1)
    return ktable


def get_enriched_kmers(fasta_fname,
                       kmer_len=4,
                       dinuc_matched=True,
                       num_shuffles=10,
                       rna=True):
    """
    Given a FASTA file, generate dinucleotide shuffles.
    """
    ktable = load_fastx_into_ktable(fasta_fname, kmer_len)
    print ktable
    # Count kmers in a set of sequence
    for n in range(0, ktable.n_entries()):
        print ktable.reverse_hash(n), ktable.get(n)
#        time.sleep(1)


if __name__ == "__main__":
    #f = os.path.expanduser("~/jaen/gff-events/mm9/with_introns/sanitized/seqs/erbin.fa")
    f = os.path.expanduser("~/jaen/gff-events/mm9/with_introns/sanitized/seqs/SE_shortest_noAceView.mm9.with_introns.event_seqs.flank_intronic_-260_-10_10_260.fa")
    #get_enriched_kmers(f, kmer_len=4)
    kmers = Kmers(5, fasta_fname=f)
    #kmers.count_kmers("./")
    #print kmers.get_counts("AAAA")
    #print kmers.get_sorted_kmers_by_count()

    import rnaseqlib.gff.gffutils_helpers as gff_helpers
    gff_helpers.get_gff_event_seqs(["chr13:104623578:104625108:-@chr13:104620268:104620411:-@chr13:104617987:104618103:-"],
                                   f,
                                   "./events.fa")
    
        


