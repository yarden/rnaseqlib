##
## Base class used to initialize the pipeline. Organizes
## and collects all the necessary files to run the pipeline
##
import os
import sys
import time

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
        pass


    def initialize(self):
        """
        Main driver function. Initialize the pipeline
        for a given genome.
        """
        print "Initializing RNA base..."
        self.download_seqs()
        self.download_tables()
        
        
