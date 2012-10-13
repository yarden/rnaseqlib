##
## Download all necessary tables
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.init as init

from rnaseqlib.init.genome_urls import *

import rnaseqlib.init.download_utils as download_utils

# Labels of UCSC tables to download
UCSC_TABLE_LABELS = ["knownGene.txt.gz",
                     "kgXref.txt.gz",
                     "knownToEnsembl.txt.gz",
                     "knownAlt.txt.gz",
                     "knownIsoforms.txt.gz"]


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
        # Download table
        table_filename = download_utils.download_url(table_url,
                                                     tables_outdir)
        # Uncompress table
        utils.gunzip_file(table_filename, tables_outdir)
        

def process_ucsc_tables(tables_outdir):
    """
    Process UCSC tables and reformat them as needed.
    """
    # Convert the UCSC knownGene format to GTF
    convert_knowngene_to_gtf(tables_outdir)
        

def convert_knowngene_to_gtf(tables_outdir):
    """
    Convert UCSC to knowngenes from genePred
    format to GTF.
    """
    # UCSC utilities genePredToGtf 
    genePredToGtf = "genePredToGtf"
    if which(geneToPred) is None:
        print "Error: need genePredToGtf from UCSC utilities " \
              "to be on path."
        sys.exit(1)
    knowngene_filename = os.path.join(tables_outdir,
                                      "knownGene.txt")
    if not os.path.isfile(knowngene_filename):
        print "Error: Cannot find %s" %(knowngene_filename)
        sys.exit(1)
    convert_cmd = "%s "

    

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
    
