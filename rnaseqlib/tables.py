##
## Download all necessary tables
##
import os
import sys
import time

import pandas

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.init as init
import rnaseqlib.genes.exons as exons
from rnaseqlib.paths import *

from rnaseqlib.init.genome_urls import *

import rnaseqlib.init.download_utils as download_utils

import numpy
from numpy import *

class GeneTable:
    """
    Parse gene table.
    """
    def __init__(self, table_dir, source):
        self.table_dir = table_dir
        self.source = source
        self.delimiter = "\t"
        self.table = None
        self.genes = None
        self.genes_list = []
        # Table indexed by gene
        self.table_by_gene = None
        # Load tables
        self.load_tables()
        

    def load_tables(self):
        """
        Load table.
        """
        if self.source == "ensGene":
            self.load_ensGene_table()
        elif self.source == "ucsc":
            self.load_ucsc_table()

            
    def load_ensGene_table(self):
        """
        Load ensGene table. Expects an 'ensGene.txt' and
        the related 'ensemblToGeneName.txt' file.

        ensGene.txt format:
        
          `bin` smallint(5) unsigned NOT NULL,
          `name` varchar(255) NOT NULL,
          `chrom` varchar(255) NOT NULL,
          `strand` char(1) NOT NULL,
          `txStart` int(10) unsigned NOT NULL,
          `txEnd` int(10) unsigned NOT NULL,
          `cdsStart` int(10) unsigned NOT NULL,
          `cdsEnd` int(10) unsigned NOT NULL,
          `exonCount` int(10) unsigned NOT NULL,
          `exonStarts` longblob NOT NULL,
          `exonEnds` longblob NOT NULL,
          `score` int(11) default NULL,
          `name2` varchar(255) NOT NULL,
          `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
          `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
          `exonFrames` longblob NOT NULL,

        ensemblToGeneName.txt format:

          `name` varchar(255) NOT NULL,
          `value` varchar(255) NOT NULL,        
        """
        self.ensGene_header = ["bin",
                               "name",
                               "chrom",
                               "strand",
                               "txStart",
                               "txEnd",
                               "cdsStart",
                               "cdsEnd",
                               "exonCount",
                               "exonStarts",
                               "exonEnds",
                               "score",
                               "name2",
                               "cdsStartStat",
                               "cdsEndStat",
                               "exonFrames"]
        self.ensemblToGeneName_header = ["name",
                                         "value"]
        ensGene_filename = os.path.join(self.table_dir,
                                        "ensGene.txt")
        if not os.path.isfile(ensGene_filename):
            print "Error: Cannot find ensGene table %s" \
                %(ensGene_filename)
            sys.exit(1)
        ensGene_name_filename = os.path.join(self.table_dir,
                                             "ensemblToGeneName.txt")
        if not os.path.isfile(ensGene_name_filename):
            print "Error: Cannot find ensemblToGeneName table %s" \
                %(ensGene_name_filename)
            sys.exit(1)
        # Load the main ensGene table
        main_table = pandas.read_table(ensGene_filename,
                                       sep=self.delimiter,
                                       names=self.ensGene_header)
        # Load ensemblToGeneName table and add this info to
        # main table
        ensGene_to_name = pandas.read_table(ensGene_name_filename,
                                            sep=self.delimiter,
                                            names=self.ensemblToGeneName_header)
        self.table = pandas.merge(main_table, ensGene_to_name)
        # Load table by gene
        self.table_by_gene = self.table.set_index("name2")#self.table#self.load_by_genes()
        # Get a genes list
        self.genes_list = self.load_ensGene_list()
        self.genes = self.load_by_genes()


    def load_by_genes(self):
        """
        Load table into gene structures.
        """
        if self.source == "ensGene":
            self.load_ensGene_by_genes()
        else:
            raise Exception, "Not implemented."


    def load_ensGene_by_genes(self):
        """
        Load ensGene table as genes.
        """
        num_genes = len(self.genes_list)
        print "Loading %d genes.." %(num_genes)
        for gene in self.genes_list:
            
            # ...
            pass

    
    def load_ensGene_list(self):
        """
        Get (a non-redundant) list of genes present in the ensGene
        table.
        """
        seen_genes = {}
        for gene in self.table_by_gene.index:
            if gene not in seen_genes:
                self.genes_list.append(gene)
                seen_genes[gene] = True
        return self.genes_list
        

    def get_const_exons(self, base_diff=5):
        if self.source == "ensGene":
            self.get_ensGene_const_exons(base_diff)
        else:
            raise Exception, "Not implemented."


    def parse_string_int_list(self, int_list_as_str,
                              delim=","):
        str_list = int_list_as_str.split(delim)
        if int_list_as_str.endswith(delim):
            # Strip off last element if list ends
            # in the delimiter we split in
            str_list = str_list[0:-1]
        ints = map(int, str_list)
        return ints


    def exon_coords_from_transcripts(self, transcripts):
        """
        Parse exons from ensGene transcript.
        """
#        exon_starts = self.parse_string_int_list(transcripts["exonStarts"].values)
#        exon_ends = self.parse_string_int_list(transcripts["exonEnds"].values)
        if type(transcripts["exonStarts"]) == str:
            exon_starts_vals = [transcripts["exonStarts"]]
            exon_ends_vals = [transcripts["exonEnds"]]
        else:
            exon_starts_vals = transcripts["exonStarts"].values
            exon_ends_vals = transcripts["exonEnds"].values
        exon_starts = map(self.parse_string_int_list, exon_starts_vals)
        exon_ends = map(self.parse_string_int_list, exon_ends_vals)
        return exon_starts, exon_ends
        

    def get_ensGene_const_exons(self, base_diff):
        """
        Load constitutive exons from ensGene table.
        """
        const_exons = []
        for gene in self.table_by_gene.index:
            if gene != "ENSMUSG00000025902": continue
            # Get transcripts for the current gene
            transcripts = self.table[self.table["name2"] == gene]
            print "TRANSCRIPTS: "
            print transcripts
            exon_starts, exon_ends = self.exon_coords_from_transcripts(transcripts)
            print exon_starts
            print exon_ends
            raise Exception
            first_transcript = transcripts.ix[0]
            rest_transcripts = transcripts.ix[1:]
            # Get the exon coordinates in the first transcript
            first_starts, first_ends = self.exon_coords_from_trans(first_transcript)
            # For each exon in the first transcript, see if it
            # is constitutive wrt to other transcripts
            for exon_start, exon_end in zip(first_starts, first_ends):
                print "Checking if: ", exon_start, exon_end
            # for trans in rest_transcripts:
            #     # Compute difference with other exons
            #     curr_starts, curr_ends = self.exon_coords_from_trans(trans)
            #     start_diffs = abs(exon_start - curr_starts)
            #     end_diffs = abs(exon_end - curr_ends)
            #     if start_diffs 
        const_exons = pandas.DataFrame(const_exons)
        return const_exons
        

    def load_ucsc_table(self):
        pass
        

#import misopy
#import misopy.gff_utils as gff_utils
#import misopy.exon_utils as exon_utils

# Labels of UCSC tables to download
UCSC_TABLE_LABELS = ["knownGene.txt.gz",
                     "kgXref.txt.gz",
                     "knownToEnsembl.txt.gz",
                     "knownAlt.txt.gz",
                     "knownIsoforms.txt.gz",
                     # Ensembl-tables
                     "ensGene.txt.gz",
                     "ensemblToGeneName.txt.gz",
                     "ensGtp.txt.gz",
                     # Refseq-tables
                     "refGene.txt.gz"]

def get_ucsc_database(genome):
    return "%s/%s/database" %(UCSC_GOLDENPATH,
                              genome)


def get_ucsc_tables_urls(genome):
    """
    Return a list of all UCSC tables URLs to download
    for a particular genome.  Format:

    [[table1_label, table1_url],
     [table2_label, table2_url],
     ...]
    """
    ucsc_database = get_ucsc_database(genome)
    table_labels = []
    for table_label in UCSC_TABLE_LABELS:
        table_url = "%s/%s" %(ucsc_database, table_label)
        table_labels.append([table_label, table_url])
    return table_labels
    

def download_ucsc_tables(genome,
                         output_dir):
    """
    Download all relevant UCSC tables for a given genome.
    """
    tables_outdir = os.path.join(output_dir, "ucsc")
    utils.make_dir(tables_outdir)
    print "Download UCSC tables..."
    print "  - Output dir: %s" %(tables_outdir)
    ucsc_tables = get_ucsc_tables_urls(genome)
    for table_label, table_url in ucsc_tables:
        print "Downloading %s" %(table_label)
        # If the table exists in uncompressed form, don't download it
        table_filename = os.path.join(tables_outdir, table_label)
        unzipped_table_fname = table_filename[0:-3]
        if os.path.isfile(unzipped_table_fname):
            print "Got %s already. Skipping download.." \
                %(unzipped_table_fname)
            continue
        # Download table
        download_utils.download_url(table_url,
                                    tables_outdir)
        # Uncompress table
        utils.gunzip_file(table_filename, tables_outdir)
        

def process_ucsc_tables(genome, output_dir):
    """
    Process UCSC tables and reformat them as needed.
    """
    tables_outdir = os.path.join(output_dir, "ucsc")
    # Convert the UCSC knownGene format to GTF
    convert_knowngene_to_gtf(tables_outdir)
    # Convert the various Ensembl tables to GFF3 format
    convert_tables_to_gff(tables_outdir)
    
    ##
    ## Load tables into gene table object
    ##
    ensGene_table = GeneTable(tables_outdir, "ensGene")
    ensGene_table.get_const_exons()
    # Compute constitutive exons and output them as files
    #output_const_exons(tables_outdir)

    
def output_const_exons(tables_dir,
                       tables_to_process=["knownGene",
                                          "ensGene",
                                          "refGene"]):
    """
    Output constitutive exons 
    """
    print "Getting constitutive exons for all tables.."
    const_exons_outdir = os.path.join(tables_dir, "const_exons")
    utils.make_dir(const_exons_outdir)
    for table in tables_to_process:
        print "Getting constitutive exons for: %s" %(table)
        gff_filename = os.path.join(tables_dir, "%s.gff3" %(table))
        if not os.path.isfile(gff_filename):
            print "Error: Cannot find %s" %(gff_filename)
            sys.exit(1)
        exons_output_filename = os.path.join(tables_dir,
                                             "%s.const_exons.gff3" %(table))
        if os.path.isfile(exons_output_filename):
            print "  - File %s exists, skipping.." %(exons_output_filename)
        # Output constitutive exons
        ##
        ## TODO: Re-write this.  Write code from scratch to get
        ## 'approximately' constitutive exons in a faster way
        ## in 'genes/exons.py'. Map transcripts to exons
        ## then compute overlap
        ## 
        tables.get_const_exons(gff_filename, exons_output_filename)
        # cds_output_filename = os.path.join(tables_dir,
        #                                    "%s.cds_exons.gff3" %(table))
        # if os.path.isfile(cds_output_filename):
        #     print "  - File %s exists, skipping.." %(cds_output_filename)
        # # Output CDS constitutive exons
        # get_const_exons_by_gene(output_filename,
        #                         const_exons_outdir,
        #                         cds_only=True,
        #                         output_filename=cds_output_filename)
        


def convert_tables_to_gff(tables_outdir):
    """
    Convert various UCSC tables to GFF3 using
    foo.
    """
    print "Converting tables to GFF3 format.."
    # Use Biotoolbox script for UCSC to GFF3 conversion
    ucsc2gff = os.path.join(SCRIPTS_DIR, "ucsc_table2gff3.pl")
    # Convert knownGene, Ensembl and RefSeq to GFF3
    tables_to_convert = ["knownGene.txt",
                         "ensGene.txt",
                         "refGene.txt"]
    t1 = time.time()
    for table in tables_to_convert:
        print "  - Converting %s to GFF" %(table)
        table_filename = os.path.join(tables_outdir,
                                      table)
        output_filename = os.path.join(tables_outdir,
                                       "%s.gff3" %(table.replace(".txt", "")))
        if os.path.isfile(output_filename):
            print "  - Found %s. Skipping conversion..." \
                %(output_filename)
            continue
        os.chdir(tables_outdir)
        table_to_gff_cmd = "%s --table %s " %(ucsc2gff,
                                              table)
        os.system(table_to_gff_cmd)
    t2 = time.time()
    print "Conversion took %.2f minutes." %((t2 - t1)/60.)
    
    
def convert_knowngene_to_gtf(tables_outdir):
    """
    Convert UCSC to knowngenes from genePred
    format to GTF.
    """
    knowngene_filename = os.path.join(tables_outdir,
                                      "knownGene.txt")
    knowngene_gtf_filename = os.path.join(tables_outdir,
                                          "knownGene.gtf")
    knowngene_gff_filename = os.path.join(tables_outdir,
                                          "knownGene.gff")
    print "Converting knownGene format to GTF..."
    if not os.path.isfile(knowngene_filename):
        print "Error: Cannot find %s" %(knowngene_filename)
        sys.exit(1)
    convert_cmd = "cat %s | cut -f1-10 | genePredToGtf file stdin %s -source=knownGene" \
        %(knowngene_filename,
          knowngene_gtf_filename)
    if not os.path.isfile(knowngene_gtf_filename):
        os.system(convert_cmd)
    gtf2gff_cmd = "gtf2gff3.pl"
    convert_gff_cmd = "%s %s > %s" %(gtf2gff_cmd,
                                     knowngene_gtf_filename,
                                     knowngene_gff_filename)
    if not os.path.isfile(knowngene_gff_filename):
        os.system(convert_gff_cmd)
    return knowngene_gtf_filename, knowngene_gff_filename
    

# import os
# import sys
# import time

# from collections import defaultdict

# import pandas

# import yklib

# class GeneTable:
#     """
#     Representing a gene table.
#     """
#     def __init__(self, settings_info):
#         """
#         Take a settings info object and load gene tables from it.
#         """
#         self.genes = {}
#         self.settings_info = settings_info
#         self.gene_table_sources = ["ensembl_genes",
# #                                   "ucsc_genes",
#                                    "refseq_genes"]
#         self.index_keys = {"ensembl_genes":
#                            "ensemblToGeneName.name"}
#         # Values to consider NA
#         self.na_values = ["n/a"]
#         # Fields needed in all tables
#         self.universal_fields = ["txStart",
#                                  "txEnd",
#                                  "chrom",
#                                  "strand",
#                                  "desc",
#                                  "geneSymbol"]                                 
#         # Fields needed to intersect events with genes
#         self.needed_fields = defaultdict(list)
#         for source in self.gene_table_sources:
#             self.needed_fields[source].extend(self.universal_fields)
#         # Ensembl specific fields
#         self.needed_fields["ensembl_genes"].extend(["name2",
#                                                     "ensemblToGeneName.value",
#                                                     "ensemblToGeneName.name"])
#         self.load_gene_table()
#         # Fields for getting relevant information from each
#         # table
#         self.table_fields = defaultdict(dict)
#         # Find the specific name for each desired field in the table
#         self.set_fields()
#         # Index the tables by their unique keys
#         self.index_gene_table()

        
#     def set_fields(self):
#         """
#         Record the relevant field names for each
#         table.
#         """
#         needed_fields = self.needed_fields
#         for table_source in self.genes:
#             for field_key in needed_fields[table_source]:
#                 # See which field matches the key
#                 found_key = False
#                 for table_field in self.genes[table_source].columns:
#                     if field_key in str(table_field):
#                         print "Setting: %s" %(field_key)
#                         self.table_fields[table_source][field_key] = table_field
#                         found_key = True
#                 if not found_key:
#                     print "Cannot find %s in: " %(field_key), self.genes[table_source].columns
            
        
#     def load_gene_table(self, tables_to_load=[]):
#         """
#         Load gene table.
#         """
#         if tables_to_load == []:
#             # By default, load all tables
#             tables_to_load = self.gene_table_sources
#         for table_source in tables_to_load:
#             print "Loading gene table: %s" %(table_source)
#             if len(self.settings_info["settings"]) == 0:
#                 print "Error: could not read from settings object: ", self.settings_info
#                 sys.exit(1)
#             if table_source in self.settings_info["settings"]:
#                 print "Loading: %s" %(table_source)
#                 table_filename = self.settings_info["settings"][table_source]
#                 table_filename = os.path.abspath(os.path.expanduser(table_filename))
#                 if not os.path.isfile(table_filename):
#                     print "Error: %s gene table filename does not exist." \
#                         %(table_filename)
#                     sys.exit(1)
#                 # Load gene table as dataframe
#                 self.genes[table_source] = pandas.read_table(table_filename,
#                                                              na_values=self.na_values)

                
#     def index_gene_table(self, table_sources=['ensembl_genes']):
#         if table_sources == None:
#             table_sources = self.gene_table_sources
#         for table_source in table_sources:
#             if table_source not in self.index_keys:
#                 print "Not indexing %s" %(table_source)
#                 continue
#             index_key = self.table_fields[table_source][self.index_keys[table_source]]
#             print "Indexing %s table by %s" %(table_source,
#                                               index_key)
#             self.genes[table_source] = self.genes[table_source].set_index(index_key)
        
    
                
#     def get_event_gene_info(self, event):
#         """
#         Get gene information from current gene tables.
#         """
#         # Get gene information for each event
#         for table_source in self.genes:
#             # Get relevant keys
#             txStart_key = self.table_fields[table_source]["txStart"]
#             txEnd_key = self.table_fields[table_source]["txEnd"]
    
