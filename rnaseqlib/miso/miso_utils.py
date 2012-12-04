##
## Utilities for processing MISO output
##
import os
import sys
import time

import pandas

import rnaseqlib
import rnaseqlib.utils as utils

def read_pe_params(insert_len_filename):
    """
    Get paired-end parameters from .insert_len file.
    """
    insert_len_filename = os.path.abspath(os.path.expanduser(insert_len_filename))
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


def load_miso_bf_file(comparisons_dir, comparison_name):
    """
    Load MISO information for a comparison name.
    """
    sample_comparison_dir = os.path.join(comparisons_dir, comparison_name)
    bf_filename = get_bf_filename(sample_comparison_dir)
    if bf_filename is None or (not os.path.isfile(bf_filename)):
        return None
    miso_bf_data = pandas.read_table(bf_filename, sep="\t")
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
    pairwise_comparison_dir = os.path.abspath(os.path.expanduser(pairwise_comparison_dir))
    bf_dir = os.path.join(pairwise_comparison_dir, "bayes-factors")
    if not os.path.isdir(bf_dir):
        print "WARNING: Could not get BF dir %s" %(bf_dir)
        return None
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




