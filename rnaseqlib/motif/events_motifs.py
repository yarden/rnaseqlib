##
## Motif analysis for alternative splicing events
##
import os
import sys
import time
import glob

import scipy
import scipy.stats

import pybedtools
import pandas

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.motif.MotifSet as motif_set
import rnaseqlib.pandas_utils as pandas_utils
import rnaseqlib.fastx_utils as fastx_utils
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

class EventSeqs:
    """
    Representation of a set of events, sequences of their various
    regions and those region coordinates.
    """
    def __init__(self, event_ids, label, input_seqs_fname,
                 remove_repeats=False,
                 entry_types=None,
                 output_dir=None):
        self.event_ids = event_ids
        self.label = label
        self.entry_types = entry_types
        self.output_dir = output_dir
        self.input_seqs_fname = input_seqs_fname
        # Whether to remove repeats or not from sequences
        self.remove_repeats = remove_repeats
        utils.make_dir(output_dir)
        # Sequence filenames for each entry type
        self.seqs_fnames = {}
        # BED filenames for each entry type
        self.bed_fnames = {}
        # Total length of sequences
        self.total_lens = {}
        self.output_event_seqs_and_coords()


    def output_event_seqs_and_coords(self):
        """
        Output the sequences for the events (FASTA) and
        their coordinates (BED format.)

        - entry_types: mapping from entry type (exon, intron) to the
          part subtypes that should be returned (e.g. se, dn, up, up_intron).

        Return the FASTA filename for its events and a BED filename.
        """
        # Output the FASTA file for the events by entry type
        if self.entry_types is None:
            return
        for entry_type in self.entry_types:
            # Output the entry type
            entry_seqs_fname = \
                os.path.join(self.output_dir,
                             "%s.%s.fa" %(self.label,
                                          entry_type))
            self.seqs_fnames[entry_type] = entry_seqs_fname
            entry_bed_fname = \
                "%s.bed" %(entry_seqs_fname.rsplit(".", 1)[0])
            self.bed_fnames[entry_type] = entry_bed_fname
            # Output FASTA for this entry type
            print "Passing %d ids" %(len(self.event_ids))
            gff_helpers.output_gff_event_seqs(self.event_ids,
                                              self.input_seqs_fname,
                                              entry_seqs_fname,
                                              entry_types=[entry_type],
                                              remove_repeats=self.remove_repeats)
            # Convert FASTA to BED coordinates
            self.total_lens[entry_type] = \
                output_bed_coords_from_fasta(entry_seqs_fname, entry_bed_fname)
            for part_type in self.entry_types[entry_type]:
                part_seqs_fname = \
                    os.path.join(self.output_dir,
                                 "%s.%s.fa" %(self.label,
                                              part_type))
                self.seqs_fnames[part_type] = part_seqs_fname
                part_bed_fname = \
                    "%s.bed" %(part_seqs_fname.rsplit(".", 1)[0])
                self.bed_fnames[part_type] = part_bed_fname
                # Output FASTA for this part type
                gff_helpers.output_gff_event_seqs(self.event_ids,
                                                  self.input_seqs_fname,
                                                  part_seqs_fname,
                                                  entry_types=[entry_type],
                                                  suffixes=[part_type],
                                                  remove_repeats=self.remove_repeats)
                # Output BED for this part type
                self.total_lens[part_type] = \
                    output_bed_coords_from_fasta(part_seqs_fname, part_bed_fname)


def parse_gff_coords(coords):
    fields = coords.split(":")
    strand = None
    chrom = fields[0]
    if len(fields) == 3:
        strand = fields[2]
    start, end = map(int, fields[1].split("-"))
    if start > end:
        # Reverse coordinates of start > end
        start, end = end, start
    return chrom, start, end, strand
    

def output_bed_coords_from_fasta(fasta_fname, bed_fname):
    """
    Output event coordinates from a FASTA file into a BED
    format.

    Assumes FASTA entry is of the form:

      >part_id:coords:entry_type
    """
    print "Converting FASTA %s to BED %s" %(fasta_fname,
                                            bed_fname)
    total_len = 0
    with open(bed_fname, "w") as bed_out:
        for fasta_entry in fastx_utils.get_fastx_entries(fasta_fname):
            fasta_name, fasta_seq = fasta_entry
            # Assume FASTA entry coordinates are in GFF format.
            # Convert them to BED
            if ";" not in fasta_name:
                raise Exception, "Malformed FASTA entry name: %s" %(fasta_name)
            gff_coords = fasta_name.split(";")[1]
            chrom, start, end, strand = parse_gff_coords(gff_coords)
            # Convert start to BED by subtracting one
            start = start - 1
            bed_entry = \
                pybedtools.create_interval_from_list(map(str, [chrom, start, end,
                                                               gff_coords, "1",
                                                               strand]))
            bed_out.write("%s" %(str(bed_entry)))
            # Accumulate total length of FASTA seqs
            total_len += len(fasta_seq)
    return total_len
    

def compare_events_motifs(exp_event_ids, control_event_ids,
                          all_events_seqs_fname,
                          output_dir,
                          kmer_lens=[4,5,6]):
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

    Also run MEME and Homer on event motifs
    """
    num_exp = len(exp_event_ids)
    num_control = len(control_event_ids)
    print "Comparing motifs between experimental and control set of events.."
    print "  - Exp events: %d" %(num_exp)
    print "  - Control events: %d" %(num_control)
    print "  - Output dir: %s" %(output_dir)
    utils.make_dir(output_dir)
    entry_types = {"exon": ["se"],
                   "intron": ["up_intron",
                              "dn_intron"]}
    remove_repeats = False
    # Experimental event set
    exp_events = EventSeqs(exp_event_ids,
                           "exp",
                           all_events_seqs_fname,
                           entry_types=entry_types,
                           output_dir=os.path.join(output_dir,
                                                   "exp"),
                           remove_repeats=remove_repeats)
    print "EXP EVENT SEQ LENS: "
    print exp_events.total_lens
    # Control event set
    control_events = EventSeqs(control_event_ids,
                               "control",
                               all_events_seqs_fname,
                               entry_types=entry_types,
                               output_dir=os.path.join(output_dir,
                                                       "control"),
                               remove_repeats=remove_repeats)
    print "CONTROL EVENT SEQ LENS: "
    print control_events.total_lens
    # Output sequences and BED for control
    for entry_type in entry_types:
        print "Processing entries of type %s" %(entry_type)
        for part_type in entry_types[entry_type]:
            curr_output_dir = os.path.join(output_dir, part_type)
            print "Running motif comparison between %s parts" %(part_type)
            print "  - Output dir: %s" %(curr_output_dir)
            # Exp seqs and coords
            exp_seqs_fname = exp_events.seqs_fnames[part_type]
            exp_coords_fname = exp_events.bed_fnames[part_type]
            # Control seqs and coords
            control_seqs_fname = control_events.seqs_fnames[part_type]
            control_coords_fname = control_events.bed_fnames[part_type]
            motif_comp = motif_set.MotifSet(exp_seqs_fname,
                                            control_seqs_fname,
                                            kmer_lens,
                                            curr_output_dir,
                                            exp_coords_fname=exp_coords_fname,
                                            control_coords_fname=control_coords_fname)
            # Get dinucleotide enriched kmers
            enriched_kmers = motif_comp.output_enriched_kmers()
            # Find differentially enriched kmers
            output_differential_kmers(enriched_kmers, curr_output_dir)
            # Run Homer
            motif_comp.find_motifs_homer(curr_output_dir)
            # Run MEME
            # ...
            

        
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


def output_differential_kmers(kmers, output_dir):
    """
    Get kmers that are enriched in exp over control.
    """
    for kmer_len in kmers:
        # Rank kmers in exp and control groups, and compare their ranks
        exp_kmers = kmers[kmer_len]["exp"].set_index("kmer")
        control_kmers = kmers[kmer_len]["control"].set_index("kmer")
        kmers_df = pandas.merge(exp_kmers, control_kmers,
                                left_index=True,
                                right_index=True,
                                suffixes=["_exp", "_control"],
                                how="outer")
        # Index results by kmer
        exp_ranks = kmers_df["target_to_shuffled_exp"].rank()
        control_ranks = kmers_df["target_to_shuffled_control"].rank()
        # Correlate the ranks
        rank_r = scipy.stats.spearmanr(exp_ranks, control_ranks)
        print "R: ", rank_r
        # Get kmers that are non-zero
        #nonzero_kmers = pandas_utils.select_df_rows(kmers_df, lambda x: x > 0,
        #                               columns=["target_to_shuffled_exp",
        #                                        "target_to_shuffled_control"],
        #                               how="all")[0]
        # Compare the fold changes
        output_fname = os.path.join(output_dir,
                                    "diff_enrichment.%d_kmer.txt" %(kmer_len))
        kmers_df["enrichment_ratio"] = \
            kmers_df["target_to_shuffled_exp"] / kmers_df["target_to_shuffled_control"]
        kmers_df["counts_ratio"] = \
            kmers_df["counts_exp"].apply(float) / kmers_df["counts_control"].apply(float)
        kmers_df.sort(columns="enrichment_ratio",
                      inplace=True,
                      ascending=False)
        kmers_df.to_csv(output_fname, sep="\t", index=True)
    


if __name__ == "__main__":
    pass
