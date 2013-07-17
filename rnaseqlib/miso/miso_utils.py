##
## Utilities for processing MISO output
##
import os
import sys
import time
import glob
import re

import pandas

import numpy as np

import rnaseqlib
import rnaseqlib.utils as utils

from collections import defaultdict


def parse_miso_counts(counts_str):
    """
    Parse two-isoform MISO counts.
    """
    counts = defaultdict(int)
    fields = re.findall("(\(.{3}\):\d+)", counts_str)
    for field in fields:
        read_class, num_reads = field.split(":")
        counts[read_class] = num_reads
    # Canonical ordering
    counts_vector = map(int, [counts["(1,0)"],
                              counts["(0,1)"],
                              counts["(1,1)"],
                              counts["(0,0)"]])
    return np.array(counts_vector,
                    dtype=np.int64)
    

def load_comparisons_counts_from_df(df,
                                    counts_labels=["sample1_counts",
                                                   "sample2_counts"]):
    """
    Return sample1 and sample2 counts from comparisons
    MISO file.
    """
    # Don't process empty dfs
    if df.empty:
        return
    # Get list of counts for each sample
    col1, col2 = counts_labels[0], counts_labels[1]
    sample1_col = "%s_int" %(col1)
    sample2_col = "%s_int" %(col2)
    df[sample1_col] = df[col1].apply(parse_miso_counts)
    df[sample2_col] = df[col2].apply(parse_miso_counts)
    return df


def get_counts_by_class(col_label, df_col, df):
    """
    Return counts for each MISO read class.
    """
    df["%s_inc_counts" %(df_col)] = \
        np.array(map(lambda x: x[0], df[col_label].values))
    df["%s_exc_counts" %(df_col)] = \
        np.array(map(lambda x: x[1], df[col_label].values))
    df["%s_const_counts" %(df_col)] = \
        np.array(map(lambda x: x[2], df[col_label].values))
    df["%s_neither_counts" %(df_col)] = \
        np.array(map(lambda x: x[3], df[col_label].values))
    return df


def read_pe_params(insert_len_filename):
    """
    Get paired-end parameters from .insert_len file.
    """
    insert_len_filename = utils.pathify(insert_len_filename)
    if not os.path.isfile(insert_len_filename):
        print "Error: %s not a file." %(insert_len_filename)
        sys.exit(1)
    insert_file = open(insert_len_filename, "r")
    fields = insert_file.readline()[1:].strip().split(",")
    pe_params = {}
    for field in fields:
        k, v = field.split("=")
        pe_params[k] = float(v)
    insert_file.close()
    return pe_params


def load_miso_bf_file(comparisons_dir, comparison_name,
                      substitute_labels=False):
    """
    Load MISO information for a comparison name.
    """
    sample_comparison_dir = os.path.join(comparisons_dir,
                                         comparison_name)
    bf_filename = get_bf_filename(sample_comparison_dir)
    if bf_filename is None or (not os.path.isfile(bf_filename)):
        return None
    miso_bf_data = pandas.read_table(bf_filename,
                                     sep="\t",
                                     index_col="event_name")
    if substitute_labels:
        # If asked, replace sample1 and sample2
        # with names of the samples
        sample1_label, sample2_label = comparison_name.split("_vs_")
        columns = []
        for c in miso_bf_data.columns:
            new_col = c.replace("sample1", sample1_label)
            new_col = new_col.replace("sample2", sample2_label)
            if (c == "bayes_factor") or (c == "diff"):
                new_col = "%s_%s" %(c, comparison_name)
            columns.append(new_col)
        miso_bf_data.columns = columns
    return miso_bf_data
    

def get_event_types_dirs(settings_info):
    """
    Return event types.
    """
    miso_events_dir = \
      utils.pathify(settings_info["settings"]["miso_events_dir"])
    event_types_dirs = [os.path.join(miso_events_dir, dirname) \
                        for dirname in os.listdir(miso_events_dir)]
    return event_types_dirs


def get_summary_filename(sample_dir):
    """
    Get summary filename from directory.
    """
    summary_dir = os.path.join(sample_dir, "summary")
    if not os.path.isdir(summary_dir):
        raise Exception, "%s not a summary dir." %(summary_dir)
    summary_files = glob.glob(os.path.join(summary_dir,
                                           "*.miso_summary"))
    summary_files = [os.path.join(summary_dir, fname) \
                     for fname in summary_files]
    if len(summary_files) > 1:
        raise Exception, "Warning: more than 1 summary file for %s" \
            %(summary_dir)
    return summary_files[0]

    
def get_bf_filename(pairwise_comparison_dir):
    """
    Return a Bayes factor filename from a
    pairwise comparisons directory.
    """
    pairwise_comparison_dir = utils.pathify(pairwise_comparison_dir)
    if not os.path.isdir(pairwise_comparison_dir):
        print "WARNING: Could not find %s" %(pairwise_comparison_dir)
        return None
    bf_dir = os.path.join(pairwise_comparison_dir, "bayes-factors")
    if not os.path.isdir(bf_dir):
        # Attempt current directory without "bayes-factor"
        # inner directory
        bf_dir = pairwise_comparison_dir
    bf_filename = glob.glob(os.path.join(bf_dir,
                                         "*.miso_bf"))
    if len(bf_filename) > 1:
        print "Error: Multiple BF filenames in %s" %(bf_dir)
        return None
    bf_filename = bf_filename[0]
    return bf_filename
    

def get_comparisons_dirs(comparisons_dir):
    """
    Get all comparisons directories.
    """
    comparisons_dirs = []
    comparisons_dir = os.path.abspath(os.path.expanduser(comparisons_dir))
    candidate_dirs = glob.glob(os.path.join(comparisons_dir, "*_vs_*"))
    for dirname in candidate_dirs:
        if not os.path.isdir(dirname):
            continue
        comparisons_dirs.append(dirname)
    return comparisons_dirs




