##
## Base class used to initialize the pipeline. Organizes
## and collects all the necessary files to run the pipeline
##
import os
import sys
import time
import glob

import rnaseqlib
import rnaseqlib.init as init
import rnaseqlib.init as init
from rnaseqlib.init import download_seqs, tables


class RNABase:
    """
    Collection of initialization files needed to run
    the pipeline on a given genome.
    """
    def __init__(self, genome, output_dir,
                 with_index=True):
        self.genome = genome
        self.output_dir = os.path.join(output_dir,
                                       genome)
        self.with_index = with_index
        

    def download_seqs(self):
        """
        Download all necessary sequences
        """
        print "Fetching sequences.."
        # Download genome sequence files
        download_seqs.download_genome_seq(self.genome,
                                          self.output_dir)
        # Download misc sequences
        download_seqs.download_misc_seqs(self.genome,
                                         self.output_dir)
        

    def download_tables(self):
        """
        Download all necessary tables.
        """
        print "Fetching tables.."
        # Download and process UCSC tables
        tables.download_ucsc_tables(self.genome,
                                    self.output_dir)
        tables.process_ucsc_tables(self.genome,
                                   self.output_dir)


    def build_indices(self):
        """
        Build relevant genome indices for use with
        Bowtie/Tophat.
        """
        if not self.with_index:
            print "Not building indices."
            return
        print "Building indices.."
        fasta_files = self.get_bowtie_index_fasta_files()
        num_files = len(fasta_files)
        print "Building Bowtie index from %d files" %(num_files)
        fasta_str = ",".join(fasta_files)
        print fasta_str
        bowtie_build_cmd = "bowtie-build %s %s"


    def get_bowtie_index_fasta_files(self):
        """
        Return a list of genome FASTA files
        to be included in the bowtie index.
        """
        genome_dir = os.path.join(self.output_dir, "genome")
        misc_dir = os.path.join(self.output_dir, "misc")
        if not os.path.isdir(genome_dir):
            print "Error: Cannot find genome directory %s" \
                %(genome_dir)
            sys.exit(1)
        # Get the genome sequence FASTA filenames
        genome_fasta_files = map(lambda f: os.path.join(genome_dir, f),
                                 glob.glob(os.path.join(genome_dir, "*.fa")))
        # Get the misc. sequence FASTA filenames
        misc_fasta_files = map(lambda f: os.path.join(misc_dir, f),
                               glob.glob(os.path.join(misc_dir, "*.fa")))
        fasta_files = genome_fasta_files + misc_fasta_files
        return fasta_files


    def initialize(self):
        """
        Main driver function. Initialize the pipeline
        for a given genome.
        """
        print "Initializing RNA base..."
        self.download_seqs()
        self.download_tables()
        self.build_indices()
        
        
