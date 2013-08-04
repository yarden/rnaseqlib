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
import glob

import argh
from argcomplete.completers import EnvironCompleter
from argh import arg

import rnaseqlib
import rnaseqlib.miso.miso_utils as miso_utils
import rnaseqlib.pandas_utils as pandas_utils
import rnaseqlib.utils as utils

import pybedtools

NA_VAL = "NA"

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
@arg("event-type", help="Event type (label that will be used to make combined file.")
@arg("output-dir", help="Output directory.")
@arg("--genes-from-gff",
     help="If passed in a GFF3 file, incorporate gene IDs/symbols from " \
          "\'gene\' entries of GFF if they are present.")
@arg("--logger-name", help="Name for logging file.")
@arg("--dry-run", help="Dry run. Do not execute commands.")
def combine_comparisons(dirname, event_type, output_dir,
                        genes_from_gff=None,
                        logger_name="misowrap_combine",
                        dry_run=False):
    """
    Given a directory containing MISO comparisons, output
    a combined file pooling information from all the comparisosn.

    - dirname: directory containing MISO comparisons to process

    - genes_from_gff: if a GFF filename, incorporate
    gene information (determined by the attributes in 'gene_id_cols')
    into the combined file.
    """
    gene_id_cols = ["ensg_id", "gsymbol"]
    print "Merging comparisons for %s" %(event_type)
    dirname = utils.pathify(dirname)
    output_dir = utils.pathify(output_dir)
    logs_dir = os.path.join(output_dir, "logs")
    if not os.path.isdir(dirname):
        logger.error("Cannot find directory %s" %(dirname))
        sys.exit(1)
    logger = utils.get_logger(logger_name, logs_dir)
    glob_path = os.path.join(dirname, "*_vs_*")
    comp_dirs = glob.glob(glob_path)
    if len(comp_dirs) == 0:
        logger.info("No comparison directories in %s. Quitting." %(dirname))
        return
    comparison_dfs = []
    for comp_num, comp_dirname in enumerate(comp_dirs):
        comparison_name = os.path.basename(comp_dirname)
        if (comp_num % 50 == 0) and comp_num > 0:
            logger.info("Loaded %d comparisons" %(comp_num))
        sample1, sample2 = comparison_name.split("_vs_")
        bf_data = miso_utils.load_miso_bf_file(dirname,
                                               comparison_name,
                                               substitute_labels=True)
        if bf_data is None:
            logger.warning("Could not find comparison %s" \
                           %(comparison_name))
            continue
        comparison_dfs.append(bf_data)
    # Merge the comparison dfs together
    combined_df = pandas_utils.combine_dfs(comparison_dfs)
    output_dir = os.path.join(output_dir, "combined_comparisons")
    utils.make_dir(output_dir)
    output_filename = os.path.join(output_dir,
                                   "%s.miso_bf" %(event_type))
    # Combine gene information if asked
    if genes_from_gff is not None:
        genes_gff_fname = genes_from_gff
        combined_df = \
          add_gene_info_to_events_df(logger, combined_df, genes_gff_fname)
    if not dry_run:
        logger.info("Outputting to: %s" %(output_filename))
        combined_df.to_csv(output_filename,
                           float_format="%.4f",
                           sep="\t",
                           na_rep=NA_VAL,
                           index=True,
                           index_label="event_name")


def add_gene_info_to_events_df(logger, df, genes_gff_fname):
    """
    Add gene information to events DataFrame.
    Assume the DataFrame is indexed by the event_name field.
    """
    logger.info("Importing genes from %s" %(genes_gff_fname))
    if not os.path.isfile(genes_gff_fname):
        logger.error("Cannot find genes GFF %s" %(genes_gff_fname))
        sys.exit(1)
    # Load mapping from events in GFF to gene info
    events_to_genes = get_events_to_genes_from_gff(genes_gff_fname)
    # Add the events to genes to the combined DataFrame
    # Assume the combined dataframe is indexed by event name
    if df.index[0] == 0:
        logger.error("Combined dataframe not indexed by event name.")
        sys.exit(1)
    gene_symbols = []
    ensembl_ids = []
    for event_name in df.index:
        gene_symbols.append(events_to_genes[event_name]["gsymbol"])
        ensembl_ids.append(events_to_genes[event_name]["ensg_id"])
    # Add gene symbol / Ensembl IDs to 
    df["ensembl_id"] = ensembl_ids
    df["gene_symbol"] = gene_symbols
    return df
    

@arg("fname", help="MISO comparisons file (*.miso_bf) to filter.")
@arg("output-dir", help="Output directory.")
@arg("event-type", help="Event type to filter.")
# Read count filters
@arg("--atleast-inc", help="TandemUTR atleast inc counts.")
@arg("--atleast-exc", help="TandemUTR atleast exc counts.")
@arg("--atleast-sum", help="TandemUTR atleast sum counts.")
def filter_comparisons(fname, output_dir, event_type,
                       atleast_inc=None,
                       atleast_exc=None,
                       atleast_sum=None,
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
