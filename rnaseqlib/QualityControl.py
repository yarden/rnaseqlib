import os
import sys
import time

import rnaseqlib
import rnaseqlib.fastq_utils as fastq_utils

import pysam

from collections import defaultdict


class QualityControl:
    """ 
    Quality control object. Defined for every
    pipeline.
    """
    def __init__(self, settings_info):
        self.settings_info = settings_info
        # QC output dir
        self.qc_outdir = self.pipeline_outdirs["qc"]
        # Number of ribosomal reads per sample
        self.num_ribo = defaultdict(int)
        # Number of mitochondrial reads per sample
        self.num_mito = defaultdict(int)
        # Number of intronic reads per sample
        self.num_intronic = defaultdict(int)
        # Number of intergenic reads per sample
        self.num_intergenic = defaultdict(int)

    def get_exon_intergenic_ratio(self, sample):
        pass

    def get_exon_intron_ratio(self, sample):
        pass

    def get_num_ribo(self, sample,
                     chr_ribo="chrRibo"):
        """
        Compute the number of ribosomal mapping reads per
        sample.

        - chr_ribo denotes the name of the ribosome containing
          chromosome.
        """
        # Call bedtools to intersect
        bamfile = pysam.Samfile(bam_filename, 'rb')
        # Retrieve all reads on the ribo chromosome
        ribo_reads = bamfile.fetch(reference=chr_ribo,
                                   start=None,
                                   end=None)
        num_ribo = 0
        # Count reads (fetch returns an iterator)
        for r in ribo_reads:
            num_ribo += 1
        return num_ribo
        

    def get_seq_cycle_profile(self, fastq_filename,
                              first_n_seqs=None):#sample):
        """
        Compute the average 'N' bases (unable to sequence)
        as a function of the position of the read.
        """
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
