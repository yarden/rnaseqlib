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
    fname = utils.pathify(fname)
    output_dir = utils.pathify(output_dir)
    print "Filtering MISO comparisons file..."
    print "  - MISO comparisons: %s" %(fname)
    print "  - Event type: %s" %(event_type)
    if event_type is not None:
        output_dir = os.path.join(output_dir, event_type)
    utils.make_dir(output_dir)
    print " - Output dir: %s" %(output_dir)
    if "UTR" in event_type:
        def_atleast_inc = tandemutr_atleast_inc
        def_atleast_exc = tandemutr_atleast_exc
        def_atleast_sum = tandemutr_atleast_sum
    elif "SE" in event_type:
        def_atleast_inc = se_atleast_inc
        def_atleast_exc = se_atleast_exc
        def_atleast_sum = se_atleast_sum
    elif "AFE" in event_type:
        def_atleast_inc = afe_atleast_inc
        def_atleast_exc = afe_atleast_exc
        def_atleast_sum = afe_atleast_sum
    elif "ALE" in event_type:
        def_atleast_inc = ale_atleast_inc
        def_atleast_exc = ale_atleast_exc
        def_atleast_sum = ale_atleast_sum
    elif "RI" in event_type:
        def_atleast_inc = ri_atleast_inc
        def_atleast_exc = ri_atleast_exc
        def_atleast_sum = ri_atleast_sum
    else:
        def_atleast_inc = 0
        def_atleast_exc = 0
        def_atleast_sum = 0
    # If read count filters are not given, use the default
    if atleast_inc is None:
        atleast_inc = def_atleast_inc
    if atleast_exc is None:
        atleast_exc = def_atleast_exc
    if atleast_sum is None:
        atleast_sum = def_atleast_sum
    # Filter the events file
    if not os.path.isfile(fname):
        print "Error: Cannot find MISO comparisons file %s" %(fname)
        sys.exit(1)
    if not fname.endswith(".miso_bf"):
        print "Warning: MISO comparisons file %s does not end in " \
              ".miso_bf.  Are you sure it is a comparisons file?" \
              %(fname)
    # Filter comparisons
    # ...
    filtered_df = None
            comparison_counts = \
                self.load_comparisons_counts_from_df(comparisons_df[event_type])
            # Get counts for each read class for sample 1 and sample 2
            comparison_counts = \
                miso_utils.get_counts_by_class("sample1_counts_int",
                                               "sample1",
                                               comparison_counts)
            comparison_counts = \
                miso_utils.get_counts_by_class("sample2_counts_int",
                                               "sample2",
                                               comparison_counts)
            filtered_df = comparison_counts
            # Filter exclusion reads
            # Only apply this to events other than TandemUTRs!
            if "TandemUTR" in event_type:
                atleast_exc = 0
                atleast_const = 5
            # Filter inclusion reads
            filtered_df = \
                filtered_df[filtered_df["sample1_inc_counts"] \
                            | filtered_df["sample2_inc_counts"] \
                            >= atleast_inc]
            # Filter exclusion reads
            filtered_df = \
                filtered_df[filtered_df["sample1_exc_counts"] \
                            | filtered_df["sample2_exc_counts"] \
                            >= atleast_exc]
            # Filter the sum of inclusion and exclusion reads
            sample1_sum = \
                filtered_df["sample1_inc_counts"] + \
                filtered_df["sample1_exc_counts"]
            sample2_sum = \
                filtered_df["sample2_inc_counts"] + \
                filtered_df["sample2_exc_counts"]
            filtered_df = \
                filtered_df[sample1_sum | sample2_sum >= atleast_sum]
            # Filter constitutive reads
            filtered_df = \
                filtered_df[filtered_df["sample1_const_counts"] \
                            | filtered_df["sample2_const_counts"] \
                            >= atleast_const]
            self.filtered_events[event_type] = filtered_df
    
    if not dry_run:
        # Call filtered comparisons here
        pass
    # Add gene information
    if get_genes_from_gff is not None:
        gene_table_fname = utils.pathify(get_genes_from_gff)
        print "Adding gene information from %s" %(gene_table_fname)
        if not os.path.isfile(gene_table_fname):
            print "Error: GFF file %s not found." %(gene_table_fname)
            sys.exit(1)
        events_to_genes = get_events_to_genes(gene_table_fname)
        
        

def main():
    argh.dispatch_commands([
        filter_comparisons,
        combine_comparisons
    ])
    

if __name__ == "__main__":
    main()
