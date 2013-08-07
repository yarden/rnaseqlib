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

import pandas

from collections import defaultdict, OrderedDict

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


def get_events_to_genes_from_gff(fname,
                                 key_names=["ID",
                                            "human_gff_id",
                                            "mouse_gff_id"],
                                 gene_id_cols=["ensg_id",
                                               "gsymbol",
                                               "refseq_id"]):
    """
    Create dictionary mapping event IDs to gene information
    from the given GFF file.

    Uses as a key each attribute name in 'key_names'. 
    """
    gff_entries = pybedtools.BedTool(fname)
    gene_entries = gff_entries.filter(lambda x: x.fields[2] == "gene")
    events_to_genes = defaultdict(dict)
    for gene in gene_entries:
        # Parse Ensembl gene, RefSeq and gene symbols
        attrs = gene.attrs
        for key_name in key_names:
            if key_name not in attrs:
                continue
            event_id = attrs[key_name]
            for gene_col in gene_id_cols:
                if gene_col in attrs:
                    events_to_genes[event_id][gene_col] = attrs[gene_col]
        for key_name in key_names:
            if key_name not in attrs: continue
            key_value = attrs[key_name]
            events_to_genes[key_value][key_name] = key_value
    return events_to_genes


def add_gene_info_to_df(df, genes_gff_fname,
                        key_names=["ID",
                                   "mouse_gff_id",
                                   "human_gff_id"],
                        gene_id_cols=["ensg_id", "gsymbol"]):
    """
    Given dataframe of comparisons add gene information
    from GFF of genes in 'genes_gff_gff'.

    Assumes 'df' is indexed by event_name. Searches
    the attributes field of all 'gene' entries in
    'genes_from_gff' to see if they have

    'key_names' are the attributes of each 'gene' entry
    that should be used as keys.
    """
    events_to_genes = \
      get_events_to_genes_from_gff(genes_gff_fname,
                                   key_names=key_names,
                                   gene_id_cols=gene_id_cols)
    # Make genes information dataframe
    gene_info_df = []
    combined_df = None
    for event_name in df.index:
        entry = {"event_name": event_name}
        possible_event_names = \
          [event_name,
           convert_event_name_to_gff_format(event_name),
           convert_event_name_to_gff_format(event_name,
                                            reverse_up_and_down=True)]
        event_found = False
        for curr_event_name in possible_event_names:
            if curr_event_name in events_to_genes:
                # Try colonifying and reversing up and down exons
                event_name = curr_event_name
                event_found = True
                break
        if not event_found:
            print "Possible events were: ", possible_event_names
            raise Exception, "Could not find event %s gene info." %(event_name)
        event_info = events_to_genes[event_name]
        for gene_col in event_info:
            entry[gene_col] = event_info[gene_col]
        # Add remaining keys' information 
        for key_name in key_names:
            if key_name in event_info:
                entry[key_name] = event_info[key_name]
            else:
                entry[key_name] = NA_VAL
        gene_info_df.append(entry)
    # Merge gene information into DataFrame
    gene_info_df = pandas.DataFrame(gene_info_df).set_index("event_name")
    combined_df = \
      pandas.merge(df, gene_info_df,
                   how="left",
                   left_index=True,
                   right_index=True)
    print combined_df.columns
    for k in key_names:
        assert (k in combined_df.columns), "Could not find key %s in df." %(k)
    return combined_df
    

@arg("dirname", help="MISO directory containing sample comparisons.")
@arg("event-type", help="Event type (label that will be used to make combined file.")
@arg("output-dir", help="Output directory.")
@arg("--logger-name", help="Name for logging file.")
@arg("--dry-run", help="Dry run. Do not execute commands.")
def combine_comparisons(dirname, event_type, output_dir,
                        logger_name="misowrap_combine",
                        dry_run=False):
    """
    Given a directory containing MISO comparisons, output
    a combined file pooling information from all the comparisons.
    
    - dirname: directory containing MISO comparisons to process
    """
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
    output_dir = os.path.join(output_dir, "combined_comparisons")
    utils.make_dir(output_dir)
    output_filename = os.path.join(output_dir,
                                   "%s.miso_bf" %(event_type))
    if os.path.isfile(output_filename):
        logger.info("Found combined comparison file %s already, skipping" \
                    %(output_filename))
        return output_filename
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
    if not dry_run:
        logger.info("Outputting to: %s" %(output_filename))
        combined_df.to_csv(output_filename,
                           float_format="%.4f",
                           sep="\t",
                           na_rep=NA_VAL,
                           index=True,
                           index_label="event_name")
    return output_filename


def convert_event_name_to_gff_format(event_name, reverse_up_and_down=False):
    """
    Given an event name like:

    'chr10:100146958-100147064:-@chr10:100148111-100148265:-@chr10:100150355-100150511:-'

    Convert it to GFF formatting using ':' instad of '-', as in:

    'chr10:100146958:100147064:-@chr10:100148111:100148265:-@chr10:100150355:100150511:-'

    Note that for minus event strands like above, the flanking exons are also
    flipped.
    """
    strand = event_name.split(":")[-1]
    def colonify(exon_name):
        chrom, coords, strand = exon_name.split(":")
        return ":".join([chrom, coords.replace("-", ":"), strand])
    if "@" not in event_name:
        return None
    # Replace ':' with '-'
    new_event_name = event_name
    fields = event_name.split("@")
    if len(fields) != 3:
        return None
    up_exon, se_exon, dn_exon = fields
    # Switch upstream and downstream exons if asked
    if reverse_up_and_down:
        up_exon, dn_exon = dn_exon, up_exon
    new_event_name = "@".join(map(colonify, [up_exon, se_exon, dn_exon]))
    return new_event_name


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
