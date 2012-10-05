import os
import sys
import time

import rnaseqlib
import rnaseqlib.fastq_utils as fastq_utils
import rnaseqlib.utils as utils

import pandas
import pysam

from collections import defaultdict

class QualityControl:
    """ 
    Quality control object. Defined for
    RNA-Seq sample.
    """
    def __init__(self, sample, pipeline):
        # Pipeline instance that the sample is attached to
        self.pipeline = pipeline
        self.sample = sample
        self.settings_info = pipeline.settings_info
        # QC results
        self.qc_results = {}
        # QC output dir
        self.qc_outdir = self.pipeline.pipeline_outdirs["qc"]
        # Number of ribosomal reads per sample
        self.num_ribo = None
        # Number of mitochondrial reads per sample
        self.num_mito = None
        # Number of intronic reads per sample
        self.num_intronic = None
        # Number of intergenic reads per sample
        self.num_intergenic = None
        # Number of reads per sample
        self.num_mapped_reads = None

    def get_num_mapped(self):
        """
        Get number of mapped reads, not counting duplicates, i.e.
        reads that have alignments in the BAM file.
        """
        return 0
    

    def get_exon_intergenic_ratio(self):
        pass
    

    def get_exon_intron_ratio(self):
        pass
    

    def get_num_ribo(self, chr_ribo="chrRibo"):
        """
        Compute the number of ribosomal mapping reads per
        sample.

        - chr_ribo denotes the name of the ribosome containing
          chromosome.
        """
        bamfile = pysam.Samfile(self.sample.bam_filename, "rb")
        # Retrieve all reads on the ribo chromosome
        ribo_reads = bamfile.fetch(reference=chr_ribo,
                                   start=None,
                                   end=None)
        num_ribo = 0
        # Count reads (fetch returns an iterator)
        for r in ribo_reads:
            num_ribo += 1
        return num_ribo
        

    def get_qc(sample):
        """
        Compile all the QC results for sample.
        """
        self.qc_results = {}
        return
    
        
    def output_qc(self):
        """
        Output QC metrics for sample.
        """
        sample_outdir = os.path.join(self.qc_outdir,
                                     self.sample.label)
        utils.make_dir(sample_outdir)
        qc_filename = os.path.join(sample_outdir,
                                   "%s.qc.txt" %(self.sample.label))
        if os.path.isfile(qc_filename):
            print "SKIPPING %s, since %s already exists..." %(self.sample.label,
                                                              qc_filename)
            return None
        # Header for QC output file for sample
        qc_headers = ["num_mapped", "num_ribo"]
        qc_entry = {"num_mapped": self.get_num_mapped(),
                    "num_ribo": self.get_num_ribo()}
        qc_df = pandas.DataFrame([qc_entry])
        # Write QC information as csv
        qc_df.to_csv(qc_filename,
                     cols=qc_headers,
                     sep="\t",
                     index=False)
        

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
