##
## Representation of gene tables
##
## Utilities for representing, download and processing
## the tables.
##
import os
import sys
import time
import csv

import pandas

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.init as init
import rnaseqlib.genes.exons as exons
import rnaseqlib.genes.GeneModel as GeneModel

from rnaseqlib.paths import *
from rnaseqlib.init.genome_urls import *
import rnaseqlib.init.download_utils as download_utils

from collections import defaultdict

import numpy
from numpy import *

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

class GeneTable:
    """
    Parse gene table.
    """
    def __init__(self, table_dir, source):
        self.table_dir = table_dir
        self.source = source
        self.delimiter = "\t"
        self.table = None
        self.genes = []
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
                                       names=self.ensGene_header,
                                       converters={"exonStarts":
                                                   self.parse_string_int_list,
                                                   "exonEnds":
                                                   self.parse_string_int_list})
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
        genes = None
        if self.source == "ensGene":
            genes = self.load_ensGene_by_genes()
        else:
            raise Exception, "Not implemented."
        return genes
        

    def load_ensGene_by_genes(self):
        print "Loading Ensembl table into genes.."
        ensGene_filename = os.path.join(self.table_dir,
                                        "ensGene.txt")
        t1 = time.time()
        table_in = csv.DictReader(open(ensGene_filename),
                                  delimiter="\t",
                                  fieldnames=self.ensGene_header)
        self.genes_table = defaultdict(list)
        for gene in table_in:
            self.genes_table[gene["name2"]].append(gene)
        # Convert genes into gene models
        for gene_id, gene_entries in self.genes_table.iteritems():
            gene = self.ensGene_entries_to_gene(gene_id, gene_entries)
            self.genes.append(gene)
        t2 = time.time()
        print "Loading took %.2f secs" %(t2 - t1)
        return self.genes



    def ensGene_entries_to_gene(self, gene_id, gene_entries):
        """
        Convert a set of ensGene entries related to a gene
        into a GeneModel. 'gene' is the gene ID and 'gene_entries'
        is a list of DictReader entries.
        """
        all_transcripts = []
        for entry in gene_entries:
            chrom = entry["chrom"]
            strand = entry["strand"]
            exon_starts = self.parse_string_int_list(entry["exonStarts"])
            exon_ends = self.parse_string_int_list(entry["exonEnds"])
            exon_coords = zip(*[exon_starts, exon_ends])
            parts = [GeneModel.Part(exon[0], exon[1], chrom, strand) \
                     for exon in exon_coords]
            transcript_id = entry["name"]
            transcript = GeneModel.Transcript(parts, chrom, strand,
                                              label=transcript_id)
            all_transcripts.append(transcript)
        gene_model = GeneModel.Gene(all_transcripts, chrom, strand,
                                    label=gene_id)
        return gene_model



    def load_ensGene_by_genes_pandas(self):
        """
        Load ensGene table as genes.
        """
        num_genes = len(self.genes_list)
        print "Loading %d genes.." %(num_genes)
        t1 = time.time()
        for gene in self.genes_list:
#            if gene != "ENSMUSG00000025902": continue
            # Get all transcripts related to gene
#            transcripts = self.table[self.table["name2"] == gene]
            transcripts = self.table_by_gene.ix["ENSMUSG00000025902"]
            pass
            # Index by transcripts
#            transcripts = transcripts.set_index("name")
#            exon_starts = transcripts["exonStarts"]
            # for trans_name in transcripts.index:
            #     curr_trans = transcripts.ix[trans_name]
            #     # Get chrom and strand
            #     chrom = curr_trans["chrom"]
            #     strand = curr_trans["strand"]
            #     exon_starts = curr_trans["exonStarts"]
            #     exon_ends = curr_trans["exonEnds"]
            #     exon_coords = zip(*[exon_starts, exon_ends])
            #     parts = [GeneModel.Part(exon[0], exon[1],
            #                             chrom, strand) \
            #              for exon in exon_coords]
            #     transcript = GeneModel.Transcript(parts, chrom, strand,
            #                                       label=trans_name)
            #     all_transcripts.append(transcript)
            #self.genes.append(transcripts)
        t2 = time.time()
        print "Loading took: %.2f secs" %(t2 - t1)
        raise Exception
            

    
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
        

    def get_const_exons(self, base_diff=6):
        for gene in self.genes:
            const_exons = gene.get_const_exons(base_diff=base_diff)
            


    def parse_string_int_list(self, int_list_as_str,
                              delim=","):
        str_list = int_list_as_str.split(delim)
        if int_list_as_str.endswith(delim):
            # Strip off last element if list ends
            # in the delimiter we split in
            str_list = str_list[0:-1]
        ints = tuple(map(int, str_list))
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
        

    def load_ucsc_table(self):
        raise Exception, "Not implemented."

    

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
    
