##
## Motif analysis for alternative splicing events
##
import os
import sys
import time
import glob

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.motif.MotifSet as motif_set

import rnaseqlib.gff.gffutils_helpers as gff_helpers

####
#### TODO: make analog of this function for motif enrichment in
####       CLIP clusters
####
#### method 1: take CLIP clusters in given gene, look for enrichment of
#### of a given kmer in that cluster by comparing its number of occurrences
#### to the number of occurrences of dinucleotide shuffled motif
####
#### method 2: use shuffleBed!  -inc the gene coordinates, and -exc the current
#### cluster coordinates per gene
####


###
def compare_events_motifs(exp_event_ids, control_event_ids,
                          all_events_seqs_fname,
                          output_dir,
                          kmer_lens=[5]):
#                          kmer_lens=[4,5,6,7,8]):
    """
    Compare the motifs in two sets of events.
    For both the experimental set and the control set, determine the motifs
    that are dinucleotide shuffled enriched.  Then compare the enriched motifs
    in either sample.

      - exp_event_ids: set of experimental event IDs
      - control_event_ids: set of control event IDs
      - all_events_seqs_fname: sequence file (FASTA) with sequences
        for all the events. Note that 'exp_event_ids' and 'control_event_ids'
        must be within this file (they are meant to be non-overlapping
        subsets of the same set of events.)
    """
    num_exp = len(exp_event_ids)
    num_control = len(control_event_ids)
    print "Comparing motifs between experimental and control set of events.."
    print "  - Exp events: %d" %(num_exp)
    print "  - Control events: %d" %(num_control)
    print "  - Output dir: %s" %(output_dir)
    utils.make_dir(output_dir)
    # Get the sequences of all the events
    seqs_dir = os.path.join(output_dir, "event_seqs")
    utils.make_dir(seqs_dir)
    exp_seqs_fname = os.path.join(seqs_dir, "exp_events.fa")
    print "Outputting exp events seqs to: %s" %(exp_seqs_fname)
    gff_helpers.output_gff_event_seqs(exp_event_ids,
                                      all_events_seqs_fname,
                                      exp_seqs_fname)
    control_seqs_fname = os.path.join(seqs_dir, "control_events.fa")
    ###
    ### TODO: MODIFY THIS TO OUTPUT SEPARATE SEQS FOR EXONS/INTRONS
    ###
    gff_helpers.output_gff_event_seqs(control_event_ids,
                                      all_events_seqs_fname,
                                      control_seqs_fname)
    # Output kmer counts for each event
    counts_dir = os.path.join(output_dir, "event_counts")
    output_event_enriched_kmers(exp_seqs_fname,
                                control_seqs_fname,
                                kmer_lens,
                                output_dir)
    # Compile counts together... make a two columns
    # format
    #
    # Ratio should be log(Experimental/Control) = log(experimental) - log(control)
    #
    # Sort the file by this value
    #
    # kmer  | counts_in_exp  dishuffle_counts_in_exp  counts_in_control  dishuffle_counts_in_control  ratio_counts_exp/ratio_counts_ctrl
    # ----------------------------------------------------------------------------------------------------------------------------------
    # kmer1 |
    # kmer2 |
    # ...
    # kmerN |
    
        

def output_event_enriched_kmers(exp_fasta_fname,
                                control_fasta_fname,
                                kmer_lens,
                                output_dir):
    """
    Output kmers enriched in events.
    """
    event_seqs = {"exp": exp_fasta_fname,
                  "control": control_fasta_fname}
    print "Outputting events kmer counts.."
    # Motifs comparison between experimental and control
    # sets of sequences
    motif_comp = motif_set.MotifSet(event_seqs["exp"],
                                    event_seqs["control"],
                                    kmer_lens,
                                    output_dir)
    motif_comp.output_enriched_kmers()


if __name__ == "__main__":
    pass
