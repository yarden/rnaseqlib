##
## Base class used to initialize the pipeline. Organizes
## and collects all the necessary files to run the pipeline
##
import os
import sys
import time
import glob
import csv

import rnaseqlib
import rnaseqlib.init as init
import rnaseqlib.utils as utils
import rnaseqlib.tables as tables
from rnaseqlib.init import download_seqs


class RNABase:
    """
    Collection of initialization files needed to run
    the pipeline on a given genome.
    """
    def __init__(self, genome, output_dir,
                 with_index=True,
                 from_dir=None,
                 init_params={}):
        self.genome = genome
        self.with_index = with_index
        self.indices_dir = None
        self.init_params = init_params
        ##
        ## Gene table names for various tasks
        ##
        self.ucsc_tables_dir = None
        # Tables to use for RPKM computation
        self.gene_table_names = ["ensGene"]#, "refSeq"]
        self.rpkm_table_names = ["ensGene",
                                 "ensGene.cds_only"]
#                                 "refSeq"]
        # Gene tables indexed by table name
        self.gene_tables = {}
        # Mapping from tables to const exons information
        self.tables_to_const_exons = {}
        self.output_dir = None
        if from_dir is None:
            self.output_dir = os.path.join(output_dir,
                                           self.genome)
        else:
            self.output_dir = from_dir
            # Load the RNABase from a given directory
            self.load_base(from_dir)
            

    def load_base(self, input_dir):
        """
        Loading RNABase from directory.
        """
        print "Loading RNA base from: %s" %(input_dir)
        if not os.path.isdir(input_dir):
            print "Error: Cannot find RNA base directory: %s" %(input_dir)
            sys.exit(1)
        self.ucsc_tables_dir = os.path.join(input_dir, "ucsc")
        self.load_rpkm_info()
        self.load_qc_info()


    def load_rpkm_info(self):
        """
        Load all information needed to compute RPKM.
        """
        self.load_const_exons_info()


    def get_const_exons_dir(self):
        const_exons_dir = os.path.join(self.ucsc_tables_dir,
                                       "exons",
                                       "const_exons")
        return const_exons_dir


    def load_gene_tables(self, tables_only=False):
        """
        Load gene information.
        """
        # Load all gene tables
        for table_name in self.gene_table_names:
            table = tables.GeneTable(self.ucsc_tables_dir,
                                     table_name,
                                     tables_only=tables_only)
            self.gene_tables[table_name] = table
        return table


    def load_const_exons_info(self):
        """
        Load constitutive exons information for all tables.
        """
        const_exons_dir = self.get_const_exons_dir()
        for table_name in self.rpkm_table_names:
            const_exons = tables.ConstExons(table_name,
                                            from_dir=const_exons_dir)
            if const_exons.found:
                self.tables_to_const_exons[table_name] = const_exons


    def load_qc_info(self):
        """
        Load all information needed to compute QC.
        """
        pass
        

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
                                   self.output_dir,
                                   init_params=self.init_params)


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
        if num_files == 0:
            print "WARNING: No FASTA files to build index from."
            return
        self.indices_dir = os.path.join(self.output_dir, "indices")
        utils.make_dir(self.indices_dir)
        ##
        ## Check if the Bowtie index is already present, if so skip
        ##
        # Check for Bowtie 1 indices
        indices = glob.glob(os.path.join(self.indices_dir,
                                         "%s*.ebwt" %(self.genome)))
        # Check for Bowtie 2 indices
        indices += glob.glob(os.path.join(self.indices_dir,
                                          "%s*.bt2"))
        if len(indices) >= 1:
            print "Found Bowtie index files in %s. Skipping index build.." \
                %(self.indices_dir)
            return
        print "Building Bowtie index from %d files" %(num_files)
        for fasta_fname in fasta_files:
            print " - %s" %(os.path.basename(fasta_fname))
        fasta_str = ",".join(fasta_files)
        # Change to indices directory
        os.chdir(self.indices_dir)
        # Use the genome as basename for the bowtie index
        bowtie_build_cmd = "bowtie-build %s %s" %(fasta_str,
                                                  self.genome)
        t1 = time.time()
        os.system(bowtie_build_cmd)
        t2 = time.time()
        print "Bowtie build took %.2f minutes" %((t2 - t1) / 60.)


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
                                 glob.glob(os.path.join(genome_dir, "chr*.fa")))
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
