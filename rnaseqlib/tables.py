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

import itertools
import operator

import pandas

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.init as init
import rnaseqlib.genes.exons as exons
import rnaseqlib.genes.GeneModel as GeneModel

from rnaseqlib.paths import *
from rnaseqlib.init.genome_urls import *
import rnaseqlib.init.download_utils as download_utils
import rnaseqlib.mapping.bedtools_utils as bedtools_utils

import misopy
import misopy.gff_utils as gff_utils

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
                     "refGene.txt.gz",
                     # tRNA tables
                     "tRNAs.txt.gz"]

class GeneTable:
    """
    Parse gene table.
    """
    def __init__(self, table_dir, source,
                 tables_only=False):
        self.table_dir = table_dir
        self.exons_dir = os.path.join(self.table_dir, "exons")
        self.const_exons_dir = os.path.join(self.exons_dir,
                                            "const_exons")
        self.introns_dir = os.path.join(self.table_dir, "introns")
        self.utrs_dir = os.path.join(self.table_dir, "utrs")
        self.source = source
        self.delimiter = "\t"
        self.table = None
        self.genes = {}
        self.genes_list = []
        self.na_val = "NA"
        # Mapping from transcripts to gene names/symbols
        self.trans_to_names = defaultdict(lambda: self.na_val)
        # Mapping from genes to gene names/symbols
        self.genes_to_names = defaultdict(lambda: self.na_val)
        # Mapping from genes to descriptions
        self.genes_to_desc = defaultdict(lambda: self.na_val)        
        # Table indexed by gene
        self.table_by_gene = None
        # UCSC known to Ensembl
        self.known_to_ensembl = defaultdict(lambda: self.na_val)
        # kgXref table
        self.kgXref_table = None
        # Initialize output directories
        self.init_dirs()
        # Load tables
        self.load_tables(tables_only=tables_only)
        

    def init_dirs(self):
        """
        Make sure directories exist.
        """
        utils.make_dir(self.exons_dir)
        utils.make_dir(self.const_exons_dir)
        utils.make_dir(self.introns_dir)
        utils.make_dir(self.utrs_dir)


    def load_kgXref_table(self):
        """
        Load kgXref table mapping UCSC transcript names to
        Ensembl transcript names.

        kgXref.txt format:

          `kgID` varchar(40) NOT NULL,
          `mRNA` varchar(40) default NULL,
          `spID` varchar(40) default NULL,
          `spDisplayID` varchar(40) default NULL,
          `geneSymbol` varchar(40) default NULL,
          `refseq` varchar(40) default NULL,
          `protAcc` varchar(40) default NULL,
          `description` longblob NOT NULL,        
        """
        self.kgXref_header = ["kgID",
                              "mRNA",
                              "spID",
                              "spDisplayID",
                              "geneSymbol",
                              "refseq",
                              "protAcc",
                              "description"]
        kgXref_filename = os.path.join(self.table_dir, "kgXref.txt")
        if not os.path.isfile(kgXref_filename):
            print "Error: Cannot find kgXref table %s" \
                %(kgXref_filename)
            sys.exit(1)
        self.kgXref_table = pandas.read_table(kgXref_filename,
                                              sep="\t",
                                              names=self.kgXref_header)
            

    def load_tables(self, tables_only=False):
        """
        Load table.
        """
        # Load kgXref for all tables
        self.load_kgXref_table()
        if self.source == "ensGene":
            self.load_ensGene_table(tables_only=tables_only)
        elif self.source == "knownGene":
            self.load_knownGene_table(tables_only=tables_only)
        elif self.source == "refSeq":
            self.load_refSeq_table(tables_only=tables_only)

        
    def load_knownGene_table(self, tables_only=False):
        raise Exception, "Not implemented."


    def load_refSeq_table(self, tables_only=False):
        raise Exception, "Not implemented."

    
    def get_ensembl_to_refseq(self):
        """
        Return mapping from Ensembl gene IDs
        to RefSeq.
        """
        # Load the combined table with kgXref
        combined_filename = os.path.join(self.table_dir,
                                         "ensGene.kgXref.combined.txt")
        if not os.path.isfile(combined_filename):
            print "Error: cannot get Ensembl to RefSeq ID mapping since " \
                  "%s does not exist." %(combined_filename)
            return None
        # Ensembl -> RefSeq mapping
        ensembl_to_refseq = defaultdict(list)
        # RefSeq -> Ensembl mapping
        refseq_to_ensembl = defaultdict(list)
        ensGene_table = csv.DictReader(open(combined_filename),
                                       delimiter="\t")
        for entry in ensGene_table:
            ensembl_id = entry["name2"]
            refseq_id = entry["refseq"]
            # If one of the ids is not available, skip the entry
            if (refseq_id == self.na_val) or (ensembl_id == self.na_val):
                continue
            # Add the entry mapping Ensembl -> Refseq and
            # RefSeq -> Ensembl only if it doesn't already exist
            if refseq_id not in ensembl_to_refseq[ensembl_id]:
                ensembl_to_refseq[ensembl_id].append(refseq_id)
            if ensembl_id not in refseq_to_ensembl[refseq_id]:
                refseq_to_ensembl[refseq_id].append(ensembl_id)
        return ensembl_to_refseq, refseq_to_ensembl

            
    def load_ensGene_table(self, tables_only=False):
        """
        Load ensGene table. Expects an 'ensGene.txt'

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

        if tables_only is True, do not parse table into
        genes but only load tables.
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
        self.knownToEnsembl_header = ["knownGene_name",
                                      "name"]
        ensGene_filename = os.path.join(self.table_dir,
                                        "ensGene.txt")
        if not os.path.isfile(ensGene_filename):
            print "Error: Cannot find ensGene table %s" \
                %(ensGene_filename)
            sys.exit(1)
        known_to_ensembl_filename = os.path.join(self.table_dir,
                                                 "knownToEnsembl.txt")
        if not os.path.isfile(known_to_ensembl_filename):
            print "Error: Cannot find knownToEnsembl table %s" \
                %(known_to_ensembl_filename)
            sys.exit(1)
        # Load the main ensGene table
        main_table = pandas.read_table(ensGene_filename,
                                       sep=self.delimiter,
                                       names=self.ensGene_header)
        self.table = main_table
#                                       converters={"exonStarts":
#                                                   self.parse_string_int_list,
#                                                   "exonEnds":
#                                                   self.parse_string_int_list})
        self.ensemblToGeneName_header = ["name",
                                         "value"]
        ensGene_name_filename = os.path.join(self.table_dir,
                                             "ensemblToGeneName.txt")
        self.ensGene_to_name_avail = True
        # Column name to use for geneSymbol
        self.gene_symbol_field = "geneSymbol"
        if not os.path.isfile(ensGene_name_filename):
            print "WARNING: Cannot find ensemblToGeneName table %s" \
                %(ensGene_name_filename)
            # Signal that ensGene to name table is unavailable
            self.ensGene_to_name_avail = False
        # If available, load ensemblToGeneName table and add this info to
        # main table
        if self.ensGene_to_name_avail:
            ensGene_to_names = pandas.read_table(ensGene_name_filename,
                                                 sep=self.delimiter,
                                                 names=self.ensemblToGeneName_header)
            # Merge names into table
            self.table = pandas.merge(main_table, ensGene_to_names,
                                      how="left")
            self.gene_symbol_field = "value"
        known_to_ensembl = pandas.read_table(known_to_ensembl_filename,
                                             sep=self.delimiter,
                                             names=self.knownToEnsembl_header)
        # Add mapping from Ensembl to gene names
        self.ensembl_to_known = known_to_ensembl.set_index("name")
        self.raw_table = self.table
        # Add mapping from Ensembl transcripts to UCSC transcripts
        self.table = pandas.merge(self.raw_table, known_to_ensembl,
                                  # try left index
                                  how="left")
#                                  how="outer")
        # Output combined table
        self.output_ensGene_combined(self.table,
                                     "ensGene.combined")
        # Bring information from kgXref
        # Note that ensGene table keys are used only in the join,
        # to avoid introducing into the table entries that have
        # kgXref info and a UCSC transcript name but *do not*
        # have an Ensembl transcript ID
        self.table = pandas.merge(self.table, self.kgXref_table,
                                  # use ensGene table keys 
                                  how="left",
                                  left_index=True,
                                  left_on=["knownGene_name"],
                                  right_on=["kgID"])
        # Output combined table with kgXref
        self.output_ensGene_combined(self.table,
                                     "ensGene.kgXref.combined")
        self.table_by_trans = self.table.set_index("name")
        # Get mapping from transcripts to genes
        self.trans_to_genes = self.table.set_index("name")
        ##
        ## Index table by gene and load a list of genes
        ## Also compute mapping from gene to symbol
        ## and gene to description
        ##
        self.table_by_gene = defaultdict(list)
        seen_genes = {}
        for idx, series in self.table.iterrows():
            gene_id = series["name2"]
            gene_info = series.to_dict()
            trans = gene_info["name"]
            self.table_by_gene[gene_id].append(gene_info)
            if gene_id in seen_genes:
                continue
            else:
                self.genes_list.append(gene_id)
                seen_genes[gene_id] = True
                # Record mapping from gene to name via transcript
                self.genes_to_names[gene_id] = gene_info[self.gene_symbol_field]
                # Get Ensembl transcript's UCSC name and from that get
                # the gene description
                self.genes_to_desc[gene_id] = gene_info["description"]
        # Parse table into actual gene objects if asked
        if not tables_only:
            self.genes = self.get_genes()

        
    def output_ensGene_combined(self, table, basename):
        """
        Output combined ensGene table.
        """
        combined_filename = os.path.join(self.table_dir,
                                         "%s.txt" %(basename))
        print "Outputting combined ensGene table..."
        print "  - Output file: %s" %(combined_filename)
        if os.path.isfile(combined_filename):
            print "Found combined file: skipping..."
            return
        table.to_csv(combined_filename,
                     sep="\t",
                     na_rep=self.na_val,
                     index=False)
            

    def load_introns(self):
        """
        Load introns.
        """
        if self.source == "ensGene":
            self.load_ensGene_introns()
        else:
            print "WARNING: Loading of introns not implemented yet " \
                  "for %s" %(self.source)


    def load_ensGene_introns(self):
        """
        Load ensGene introns.

        For each transcript, compute the coordinates of introns between the exons.
        """
        pass
        

    def get_genes(self):
        """
        Load table into gene structures.

        Return a generator.
        """
        genes = None
        if self.source == "ensGene":
            genes = self.get_ensGene_by_genes()
        else:
            raise Exception, "Not implemented."
        return genes
        

    def load_genes_by_id(self):
        genes_by_id = {}
        genes = self.get_genes()
        for gene in genes:
            genes_by_id[gene.label] = gene
        return genes_by_id
            

    def get_ensGene_by_genes(self):
        print "Loading Ensembl table into genes..."
        # Main ensGene table
        ensGene_filename = os.path.join(self.table_dir,
                                        "ensGene.txt")
        t1 = time.time()
        table = dictread_groupby_col(open(ensGene_filename),
                                     "name2",
                                     delimiter="\t",
                                     fieldnames=self.ensGene_header)
        for gene_id, gene_entries in table.iteritems():
            gene_symbol = self.na_val
            all_transcripts = []
            entries_copy = []
            for entry in gene_entries:
                entries_copy.append(entry)
                chrom = entry["chrom"]
                strand = entry["strand"]
                # Convert start coordinates into 1-based coordinates
                exon_starts = (int(start) + 1 for start in entry["exonStarts"].rstrip(",").split(","))
                exon_ends = (int(end) for end in entry["exonEnds"].rstrip(",").split(","))
                exon_coords = zip(exon_starts, exon_ends)
                # Convert cds coordinates into 1-based as well
                cds_start = int(entry["cdsStart"]) + 1
                cds_end = int(entry["cdsEnd"])
                transcript_id = entry["name"]
                parts = [GeneModel.Part(exon[0], exon[1], chrom, strand,
                                        parent=transcript_id) \
                         for exon in exon_coords]
                transcript = GeneModel.Transcript(parts, chrom, strand,
                                                  label=transcript_id,
                                                  cds_start=cds_start,
                                                  cds_end=cds_end,
                                                  parent=gene_id)
                gene_symbol = self.trans_to_names[transcript_id]
                all_transcripts.append(transcript)
            gene_model = GeneModel.Gene(all_transcripts, chrom, strand,
                                        label=gene_id,
                                        gene_symbol=gene_symbol)
            self.genes[gene_id] = gene_model
        t2 = time.time()
        print "Loading took %.2f secs" %(t2 - t1)
        return self.genes


    def load_ensGene_name_table(self, delimiter="\t"):
        """
        Load mapping from genes to names.
        """
        # Table mapping ensembl gene IDs to gene symbols
        ensembl_to_name_filename = os.path.join(self.table_dir,
                                                "ensemblToGeneName.txt")
        name_table = open(ensembl_to_name_filename, "r")
        for line in name_table:
            fields = line.strip().split(delimiter)
            key_header, val_header = self.ensemblToGeneName_header
            self.trans_to_names[key_header] = val_header
        name_table.close()
        return self.trans_to_names

        
    def output_exons_as_gff(self,
                            base_diff=6,
                            const_only=False,
                            cds_only=False):
        """
        Output constitutive exons for all genes as GFF.

        - const_only: if True, output only constitutive exons
        - cds_only: if True, output CDS only exons 
        """
        exons_type = "exons"
        exons_outdir = self.exons_dir
        if const_only:
            exons_type = "const_exons"
            exons_outdir = self.const_exons_dir
        if cds_only:
            exons_basename = "%s.cds_only.%s.gff" %(self.source, exons_type)
        else:
            exons_basename = "%s.%s.gff" %(self.source, exons_type)
        gff_output_filename = os.path.join(exons_outdir, exons_basename)
        print "Outputting exons..."
        print "  - Exons type: %s" %(exons_type)
        print "  - Output file: %s" %(gff_output_filename)
        print "  - CDS only: %s" %(cds_only)
        if os.path.isfile(gff_output_filename):
            print "%s exists. Skipping.." %(gff_output_filename)
            return
        # Output a map from genes to constitutive exons
        # for convenience
        genes_to_exons_fname = os.path.join(exons_outdir,
                                            exons_basename.replace(".gff",
                                                                   ".to_genes.txt"))
        gff_out = gff_utils.Writer(open(gff_output_filename, "w"))
        rec_type = "exon"
        genes_to_exons = []
        genes_to_exons_header = ["gene_id", "exons"]
        for gene_id, gene in self.genes.iteritems():
            if const_only:
                # Get only constitutive exons
                exons = gene.compute_const_exons(base_diff=base_diff,
                                                 cds_only=cds_only)
            elif cds_only:
                # Get all CDS exons
                exons = gene.cds_parts
            else:
                # Get all exons
                exons = gene.parts
            exon_labels = [e.label for e in exons]
            if len(exon_labels) == 0:
                exon_labels = self.na_val
            else:
                exon_labels = ",".join(exon_labels)
            entry = {"gene_id": gene_id,
                     "exons": exon_labels}
            genes_to_exons.append(entry)
            # Output constitutive exons to GFF file
            GeneModel.output_parts_as_gff(gff_out,
                                          exons,
                                          gene.chrom,
                                          gene.strand,
                                          source=self.source,
                                          rec_type=rec_type,
                                          gene_id=gene_id)
        genes_to_exons = pandas.DataFrame(genes_to_exons)
        genes_to_exons.to_csv(genes_to_exons_fname,
                              cols=genes_to_exons_header,
                              index=False,
                              sep="\t")


    def output_exons_as_bed(self):
        """
        Output the table's exons as BED.
        
        Only implemented for ensGene.txt; probably not
        necessary to work out for other tables since this is
        only used for aggregate statistics.
        """
        if self.source != "ensGene":
            return
        output_filename = os.path.join(self.exons_dir,
                                       "%s.exons.bed" %(self.source))
        print "Outputting exons..."
        if os.path.isfile(output_filename):
            print "  - Found %s. Skipping..." %(output_filename)
            return output_filename
        exons_file = open(output_filename, "w")
        for idx, series in self.raw_table.iterrows():
            gene_info = series.to_dict()
            gene_id = gene_info["name2"]
            chrom = gene_info["chrom"]
            strand = gene_info["strand"]
            exonStarts = gene_info["exonStarts"]
            exonEnds = gene_info["exonEnds"]
            # Keep 0-based start of ensGene table since
            # this will be outputted as a BED
            exon_starts = (int(start) for start in exonStarts.rstrip(",").split(","))
            exon_ends = (int(end) for end in exonEnds.rstrip(",").split(","))
            exon_coords = zip(exon_starts, exon_ends)
            # Output as BED: encode gene ID
            bedtools_utils.output_intervals_as_bed(exons_file,
                                                   chrom, exon_coords, strand,
                                                   name=gene_id)
        exons_file.close()
        return output_filename


    def output_merged_exons(self):
        """
        Output the table's merged exons as a (sorted) BED file.

        Used to determine the exonic content of a sample. Relies on
        sortBed and mergeBed to do the heavy lifting.

        Only implemented for ensGene.txt; probably not
        necessary to work out for other tables since this is
        only used for aggregate statistics.
        """
        if self.source != "ensGene":
            return
        print "Outputting merged exons..."
        exons_filename = os.path.join(self.exons_dir,
                                      "ensGene.exons.bed")
        if not os.path.isfile(exons_filename):
            print "Error: Could not find exons filename %s" %(exons_filename)
            print "Did a previous step fail?"
            sys.exit(1)
        output_filename = os.path.join(self.exons_dir,
                                       "ensGene.merged_exons.bed")
        # Merge the exons and output them as a new BED file
        bedtools_utils.merge_bed(exons_filename, output_filename)


    def load_merged_exons_by_gene(self):
        """
        Load ensGene exons from a BED file, indexed by gene.
        """
        self.merged_exons_header = ["chrom",
                                    "start",
                                    "end",
                                    "name",
                                    "strand"]
        merged_exons_filename = os.path.join(self.exons_dir,
                                             "ensGene.merged_exons.bed")
        merged_exons_file = open(merged_exons_filename)
        ensGene_bed = csv.DictReader(merged_exons_file,
                                     fieldnames=self.merged_exons_header,
                                     delimiter="\t")
        merged_exons_by_gene = defaultdict(list)
        for entry in ensGene_bed:
            # Gene names might appear multiple times due in mergeBed
            # output in which case they are semi-colon delimited
            gene_id = entry["name"].split(";")[0]
            merged_exons_by_gene[gene_id].append(entry)
        return merged_exons_by_gene


    def output_introns(self, min_intron_size=50):
        """
        Given the merged exons, compute the intronic coordinates.
        
        Only implemented for ensGene.txt; probably not
        necessary to work out for other tables since this is
        only used for aggregate statistics.

        Exclude intronic content that is less than 'min_intron_size'.
        """
        if self.source != "ensGene":
            return
        output_filename = os.path.join(self.introns_dir,
                                       "ensGene.introns.bed")
        print "Outputting introns..."
        if os.path.isfile(output_filename):
            print "  - Found %s. Skipping..." %(output_filename)
            return
        print " - Output file: %s" %(output_filename)
        introns_file = open(output_filename, "w")
        # Load ensGene exons
        merged_exons_by_gene = self.load_merged_exons_by_gene()
        for gene_id, merged_exons in merged_exons_by_gene.iteritems():
            chrom = merged_exons[0]["chrom"]
            strand = merged_exons[0]["strand"]
            # For each gene, get its list of introns and serialize them
            exon_coords = [(int(exon["start"]), int(exon["end"])) \
                           for exon in merged_exons]
            intron_coords = []
            for first_exon, second_exon in zip(exon_coords, exon_coords[1::1]):
                # Intron start coordinate is the coordinate right after
                # the end of the first exon, intron end coordinate is the
                # coordinate just before the beginning of the second exon
                intron_start = first_exon[1] + 1
                intron_end = second_exon[0] - 1
                if intron_start >= intron_end:
                    continue
                # Filter on intron size (in 0-based coordinates)
                intron_size = intron_end - intron_start
                if intron_size < min_intron_size:
                    continue
                intron_coords.append((intron_start, intron_end))
            bedtools_utils.output_intervals_as_bed(introns_file,
                                                   chrom,
                                                   intron_coords,
                                                   strand,
                                                   name=gene_id)
        introns_file.close()
                                                   

    def parse_string_int_list(self, int_list_as_str,
                              delim=","):
        """
        Parse a comma-separated string list into a
        list of integers.
        """
        str_list = int_list_as_str.split(delim)
        # Strip off last element if list ends
        # in the delimiter we split in
        ints = tuple(map(int, str_list[0:-1]))
        return ints


class ConstExons:
    """
    A table storing a set of constitutive exons.

    Consists of a GFF filename specifying the exons
    and a text file mapping genes to their constitutive exons.
    """
    def __init__(self, table_name,
                 from_dir=None):
        self.table_name = table_name
        self.from_dir = from_dir
        self.gff_filename = None
        self.na_val = "NA"
        self.genes_to_exons_filename = None
        self.found = False
        # A list of genes to exons mapping
        self.genes_to_exons = []
        self.exon_lens = defaultdict(int)
        if from_dir is not None:
            self.load_const_exons()


    def load_const_exons(self):
        """
        Load constitutive exons from a table name.
        """
        print "Loading constitutive exons for %s from dir: %s" \
            %(self.table_name,
              self.from_dir)
        self.gff_filename = os.path.join(self.from_dir,
                                         "%s.const_exons.gff" %(self.table_name))
        if not os.path.isfile(self.gff_filename):
            print "WARNING: Cannot find constitutive exons GFF for %s" \
                %(self.table_name)
            return
        # Look for the mapping from genes to exons
        self.genes_to_exons_filename = os.path.join(self.from_dir,
                                                    "%s.const_exons.to_genes.txt" \
                                                    %(self.table_name))
        if not os.path.isfile(self.genes_to_exons_filename):
            print "WARNING: Cannot find mapping from genes to constitutive " \
                "exons for %s" %(table_to_get)
            return
        # Load genes to exons mapping
        self.load_genes_to_exons()
        self.found = True


    def load_genes_to_exons(self):
        """
        Return a dictionary mapping each table name to its file
        that specifies a mapping from genes to constitutive
        exons.
        """
        table_file = open(self.genes_to_exons_filename, "r")
        table_in = csv.DictReader(table_file,
                                  delimiter="\t")
        for entry in table_in:
            self.genes_to_exons.append(entry)
            # Compute the length of each set of exons
            if entry["exons"] == self.na_val:
                continue
            exons = entry["exons"].split(",")
            exon_coords = map(lambda e: e.split(":")[1].split("-"), exons)
            exon_lens = map(lambda coords: int(coords[1]) - int(coords[0]) + 1,
                            exon_coords)
            for exon, exon_len in itertools.izip(exons, exon_lens):
                self.exon_lens[exon] = exon_len
                

    def __repr__(self):
        return "ConstExons(table=%s, gff=%s, genes_to_exons=%d entries)" \
            %(self.table_name,
              self.gff_filename,
              len(self.genes_to_exons))
        
        


##
## Related table utilities
##
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
        download_status = download_utils.download_url(table_url,
                                                      tables_outdir)
        if download_status is None:
            print "Failed to get %s, skipping.." %(table_label)
            continue
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
    ## Process misc. tables
    ##
    # tRNA table
    process_tRNA_table(tables_outdir)
    # snoRNA table
    # ...
    # mitoRNA table
    # ...
    ##
    ## Process gene tables
    ##
    table_names = ["ensGene"]#, "refGene"]
    for table_name in table_names:
        table = GeneTable(tables_outdir, table_name)
        # Output the table's exons as GFF
        table.output_exons_as_gff()
        # Output the table's CDS-only exons as GFF
        table.output_exons_as_gff(cds_only=True)        
        # Output the table's exons as BED
        table.output_exons_as_bed()
        # Output the table's merged exons
        table.output_merged_exons()
        # Output the table's constitutive exons
        table.output_exons_as_gff(const_only=True)
        # Output the table's CDS-only constitutive exons
        table.output_exons_as_gff(const_only=True,
                                  cds_only=True)
        # Output introns
        table.output_introns()


def process_tRNA_table(tables_outdir,
                       delimiter="\t"):
    """
    Process tRNA table from UCSC.
    
        CREATE TABLE `tRNAs` (
          `bin` smallint(5) unsigned NOT NULL,
          `chrom` varchar(255) NOT NULL,
          `chromStart` int(10) unsigned NOT NULL,
          `chromEnd` int(10) unsigned NOT NULL,
          `name` varchar(255) NOT NULL,
          `score` int(10) unsigned NOT NULL,
          `strand` char(1) NOT NULL,
          `aa` varchar(255) NOT NULL,
          `ac` varchar(255) NOT NULL,
          `intron` varchar(255) NOT NULL,
          `trnaScore` float NOT NULL,
          `genomeUrl` varchar(255) NOT NULL,
          `trnaUrl` varchar(255) NOT NULL,
          PRIMARY KEY  (`chrom`,`chromStart`,`chromEnd`),
          KEY `binChrom` (`chrom`,`bin`)
        ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
        SET character_set_client = @saved_cs_client;
    """
    # Process tRNA table
    tRNA_filename = os.path.join(tables_outdir, "tRNAs.txt")
    tRNA_header = ["bin",
                   "chrom",
                   "chromStart",
                   "chromEnd",
                   "name",
                   "score",
                   "strand",
                   "aa",
                   "ac",
                   "intron",
                   "trnaScore",
                   "genomeUrl",
                   "trnaUrl"]
    if not os.path.isfile(tRNA_filename):
        print "WARNING: Could not find tRNA table %s" \
            %(tRNA_filename)
        return
    tRNA_outdir = os.path.join(tables_outdir, "tRNAs")
    utils.make_dir(tRNA_outdir)
    tRNA_bed_filename = os.path.join(tRNA_outdir,
                                     "tRNAs.bed")
    print "Writing tRNA BED.."
    print "  - Output file: %s" %(tRNA_bed_filename)
    if os.path.isfile(tRNA_bed_filename):
        print "Found %s. Skipping.." %(tRNA_bed_filename)
        return
    # Output tRNA table as BED
    tRNA_table = csv.DictReader(open(tRNA_filename, "r"),
                                delimiter=delimiter,
                                fieldnames=tRNA_header)
    tRNA_bed = open(tRNA_bed_filename, "w")
    for entry in tRNA_table:
        bed_fields = [entry["chrom"],
                      entry["chromStart"],
                      entry["chromEnd"],
                      entry["name"],
                      entry["trnaScore"],
                      entry["strand"]]
        bed_line = "%s\n" %("\t".join(bed_fields))
        tRNA_bed.write(bed_line)
    tRNA_bed.close()
    

def convert_tables_to_gff(tables_outdir):
    """
    Convert various UCSC tables to GFF3 using
    foo.
    """
    print "Converting tables to GFF3 format.."
    # Use Biotoolbox script for UCSC to GFF3 conversion
    ucsc2gff = "ucsc_table2gff3.pl"
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


##
## Random table utilities
##
def dictread_groupby_col(file_in, col,
                         delimiter="\t",
                         fieldnames=None):
    """
    Read entries via DictReader grouped by column 'col'.

    Returns a defaultdict of lists.
    """
    table_in = csv.DictReader(file_in,
                              delimiter=delimiter,
                              fieldnames=fieldnames)
    table = defaultdict(list)
    for row in table_in:
        table[row[col]].append(row)
    return table
                          
