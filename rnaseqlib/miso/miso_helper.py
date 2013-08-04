##
## MISO helper script
##
## Alternative interface to MISOWrap. Meant for when you don't want to
## define a settings file.
##
##
## Filter MISO comparisons
##
import os
import sys
import time

import argh
from argcomplete.completers import EnvironCompleter
from argh import arg

import rnaseqlib
import rnaseqlib.utils as utils

import pybedtools


DEFAULT_FILTERS = \
    {
    "SE":
     # Default counts filters
     # SE
     {"atleast_inc": 10,
      "atleast_exc": 1,
      "se_atleast_sum": 20},
    "AFE":
     {"atleast_inc": 10,
      "atleast_exc": 10,
      "afe_atleast_sum": 10},
    "ALE":
     {"atleast_inc": 10,
      "atleast_exc": 10,
      "atleast_sum": 20},
    "TandemUTR":
     {"atleast_inc": 5,
      "atleast_exc": 0,
      "atleast_sum": 10},
    "RI":
     {"atleast_inc": 10,
      "atleast_exc": 1,
      "atleast_sum": 20}
    }


def get_events_to_genes_from_gff(fname):
    """
    Create dictionary mapping event IDs to gene information
    from the given GFF file.
    """
    gff_entries = pybedtools.BedTool(fname)
    gene_entries = gff_entries.filter(lambda x: x.fields[2] == "gene")
    events_to_genes = {}
    for gene in gene_entries:
        # Parse Ensembl gene, RefSeq and gene symbols
        attrs = gene.attrs
        event_id = attrs["ID"]
        events_to_genes[event_id] = {}
        if "ensg_id" in attrs:
            events_to_genes[event_id] = attrs["ensg_id"]
        if "gsymbol" in attrs:
            events_to_genes[event_id] = attrs["gsymbol"]
        if "refseq_id" in attrs:
            events_to_genes[event_id] = attrs["refseq_id"]
    return events_to_genes

@arg("dirname", help="MISO directory containing sample comparisons.")
@arg("output-dir", help="Output directory.")
@arg("--get-genes-from-gff",
     help="If passed in a GFF3 file, incorporate gene IDs/symbols from " \
          "\'gene\' entries of GFF if they are present.")
@arg("--dry-run", help="Dry run. Do not execute commands.")
def combine_comparisons(fname, output_dir,
                        get_genes_from_gff=None,
                        gene_id_cols=["ensg_id", "gsymbol"],
                        dry_run=False):
    """
    Given a directory containing MISO comparisons, output
    a combined file pooling information from all the comparisosn.

    - get_genes_from_gff: if a GFF filename, incorporate
    gene information (determined by the attributes in 'gene_id_cols')
    into the combined file.
    """
    pass



@arg("fname", help="MISO comparisons file (*.miso_bf) to filter.")
@arg("output-dir", help="Output directory.")
@arg("event-type", help="Event type to filter.")
# Read count filters
@arg("--atleast-inc", help="TandemUTR atleast inc counts.")
@arg("--atleast-exc", help="TandemUTR atleast exc counts.")
@arg("--atleast-sum", help="TandemUTR atleast sum counts.")
# Gene table-related arugments
# GFF to read gene symbols from
@arg("--get-genes-from-gff",
     help="If passed in a GFF3 file, incorporate gene IDs/symbols from " \
          "\'gene\' entries of GFF if they are present.")
@arg("--gene-table",
     help="If given an ensGene combined table (which contains kgXref) " \
          "add gene description from table into the file.")
@arg("--dry-run", help="Dry run. Do not execute commands.")
def filter_comparisons(fname, output_dir,
                       event_type=None,
                       atleast_inc=None,
                       atleast_exc=None,
                       atleast_sum=None,
                       gene_table=None,
                       gene_id_cols=["ensg_id", "gsymbol"],
                       dry_run=False):
    """
    Filter a MISO comparison file (*.miso_bf)
    Annotate a GFF file with useful information. For now, add annotation
    of gene IDs based on an input GFF annotation of genes.

    Computes the most inclusive transcription start/end coordinates
    fonr each gene, and then uses pybedtools to intersect (in strand-specific 
    manner) with the input annotation.
    """
    pass
        

def main():
    argh.dispatch_commands([
        filter_comparisons,
        combine_comparisons
    ])
    

if __name__ == "__main__":
    main()
