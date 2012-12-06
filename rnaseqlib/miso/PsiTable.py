##
## PsiTable: class for representing MISO results
## from a set of samples and parsing their results
##
import os
import re
import sys
import time
import glob

import misopy
import pandas

import numpy as np
from collections import defaultdict

import rnaseqlib
import rnaseqlib.tables as tables
import rnaseqlib.miso.miso_utils as miso_utils

from miso_utils import \
     get_summary_filename, \
     get_comparisons_dirs, \
     get_bf_filename


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
    

class PsiTable:
    """
    Representation of Psi values from a set of samples.

    Takes as input a MISOWrap object.
    """
    def __init__(self, misowrap_obj,
                 verbose=True):
        self.delimiter = "\t"
        self.misowrap_obj = misowrap_obj
        self.settings_info = self.misowrap_obj.settings_info
        self.sample_labels = self.misowrap_obj.sample_labels
        # Where MISO output for samples is
        self.miso_outdir = self.misowrap_obj.miso_outdir
        # Where the MISO sample comparisons are
        self.comparisons_dir = self.misowrap_obj.comparisons_dir
        self.filtered_events = {}
        # Load gene tables
        self.gene_table = None
        self.load_gene_table()
        self.verbose = verbose
        # Event types to process
        self.event_types = self.misowrap_obj.event_types
        # Summaries dataframe
        self.summaries_df = None
        # Comparisons dataframe
        self.comparisons_df = None
        ##
        ## Load annotation of events, like a map
        ## events to genes.
        ##
        self.events_to_genes = {}
        self.load_events_to_genes()
        ## Load the MISO output
        ##
        # Load summaries
        self.load_summaries(self.miso_outdir)
        # Load comparisons
        self.load_comparisons(self.comparisons_dir)
        # Filter events based on coverage
        self.filter_coverage_events()


    def load_events_to_genes(self):
        """
        Load mapping from events to genes.
        """


    def load_gene_table(self):
        """
        Try to load a gene table if one is given.
        """
        # Load gene table based on settings info
        self.gene_table = tables.GeneTable(self.misowrap_obj.tables_dir,
                                           "ensGene",
                                           tables_only=True)

        
    def add_event_genes(self):
        """
        Add the gene information to each event.
        """
        pass
    

    def load_summaries(self, miso_samples_dir):
        """
        Load MISO summary files.
        """
        print "Loading summary files.."
        summaries_dict = defaultdict(dict)
        for sample in self.sample_labels:
            for event_type in self.event_types:
                sample_name, label = sample
                sample_dir = os.path.join(miso_samples_dir,
                                          sample_name,
                                          event_type)
                if not os.path.isdir(sample_dir):
                    print "WARNING: Skipping %s..." \
                        %(sample_dir)
                    continue
                summary_filename = get_summary_filename(sample_dir)
                if not os.path.isfile(summary_filename):
                    print "WARNING: %s not a summary file" \
                        %(summary_filename)
                    continue
                summary_df = pandas.read_table(summary_filename,
                                               sep=self.delimiter)
                summaries_dict[event_type][sample_name] = summary_df
        self.summaries_df = pandas.DataFrame(summaries_dict)


    def filter_only_two_isoform(self, df,
                                col_label="sample1_posterior_mean"):
        """
        Given DataFrame, keep only entries with Psi values
        for two isoforms.
        """
        where_scalar = map(lambda x: "," not in str(x),
                           df[col_label])
        df_filtered = df.ix[where_scalar]
        # Ensure strings are cast to float
        df_filtered = self.set_datatypes(df_filtered)
        return df_filtered

    
    def set_datatypes(self, df,
                      float_cols=["sample1_posterior_mean",
                                  "sample2_posterior_mean",
                                  "sample1_ci_low",
                                  "sample1_ci_high",
                                  "sample2_ci_low",
                                  "sample2_ci_high",
                                  "bayes_factor",
                                  "diff"]):
        for col in float_cols:
            df[col] = map(float, df[col])
        return df


    def load_comparisons(self, comparisons_dir,
                         only_two_isoform=True):
        """
        Load MISO comparisons files.
        """
        print "Loading comparisons.."
        comparisons_dict = defaultdict(list)
        dataframe_dict = {}
        for event_type in self.event_types:
            event_comparisons_dir = os.path.join(comparisons_dir,
                                                 event_type)
            comparisons_dirnames = get_comparisons_dirs(event_comparisons_dir)
            comparison_labels = []
            if len(comparisons_dirnames) == 0:
                print "WARNING: No comparisons for event type %s in %s" \
                    %(event_type, event_comparisons_dir)
                continue
            for curr_comp_dir in comparisons_dirnames:
                comparison_label = os.path.basename(curr_comp_dir)
                comparison_labels.append(comparison_label)
                bf_filename = get_bf_filename(curr_comp_dir)
                if not os.path.isfile(bf_filename):
                    raise Exception, "No BF file: %s" %(bf_filename)
                curr_df = pandas.read_table(bf_filename,
                                            sep=self.delimiter,
                                            index_col=[0])
                # Keep only events for which we have two isoforms
                if only_two_isoform:
                    if event_type == "TandemUTR_3pseq":
                        curr_df = self.filter_only_two_isoform(curr_df)
                comparisons_dict[event_type].append(curr_df)
            # Concatenate all the DataFrames for each comparison together
            dataframe_dict[event_type] = \
                pandas.concat(comparisons_dict[event_type],
                              keys=comparison_labels)
        self.comparisons_df = dataframe_dict


    def load_comparisons_counts_from_df(self, df,
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


    def filter_coverage_events(self, comparisons_df=None,
                               atleast_inc=1,
                               atleast_exc=5,
                               atleast_sum=20,
                               atleast_const=1):
        """
        Filter events for coverage.
        """
        print "filter_coverage_events::Filtering..."
        if comparisons_df == None:
            comparisons_df = self.comparisons_df
        if len(comparisons_df.keys()) == 0:
            print "Not filtering - no comparisons found."
            return
        for event_type in self.event_types:
            if event_type not in comparisons_df:
                continue
            ##
            ## Load read count filters from the settings
            ##
            if event_type in self.misowrap_obj.event_filters:
                event_filters = self.misowrap_obj.event_filters[event_tyoe]
                if "atleast_inc" in event_filters:
                    atleast_inc = event_filters["atleast_inc"]
                if "atleast_exc" in event_filters:
                    atleast_exc = event_filters["atleast_exc"]
                if "atleast_sum" in event_filters:
                    atleast_sum = event_filters["atleast_sum"]
                if "atleast_const" in event_filters["atleast_const"]:
                    atleast_const = event_filters["atleast_const"]
            print "Filtering event type: %s" %(event_type)
            comparison_counts = \
                self.load_comparisons_counts_from_df(comparisons_df[event_type])
            # Get counts for each read class for sample 1 and sample 2
            comparison_counts = self.get_counts_by_class("sample1_counts_int",
                                                         "sample1",
                                                         comparison_counts)
            comparison_counts = self.get_counts_by_class("sample2_counts_int",
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

        
    def get_counts_by_class(self, col_label, df_col, df):
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


    def output_filtered_comparisons(self, output_dir=None,
                                    sort_column="bayes_factor",
                                    columns_to_write=[#"event_name",
                                                      "sample1_posterior_mean",
                                                      "sample1_ci_low",
                                                      "sample1_ci_high",
                                                      "sample2_posterior_mean",
                                                      "sample2_ci_low",
                                                      "sample2_ci_high",
                                                      "diff",
                                                      "bayes_factor",
                                                      "isoforms",
                                                      "sample1_counts",
                                                      "sample1_assigned_counts",
                                                      "sample2_counts",
                                                      "sample2_assigned_counts",
                                                      "chrom",
                                                      "strand",
                                                      "mRNA_starts",
                                                      "mRNA_ends"],
                                    with_gene_table=True):
        """
        Output filtered comparisons table to.
        """
        if output_dir == None:
            output_dir = self.output_dir
        # Output each file by event type
        output_dir = os.path.join(output_dir, "test_filtered_events")
        print "output_filtered_comparisons::writing to dir: %s" %(output_dir)
        for event_type, filtered_df in self.filtered_events.iteritems():
            curr_output_dir = os.path.join(output_dir, event_type)
            print "Event type: %s" %(event_type)
            # View by comparison
            comparison_labels = \
                utils.unique_list(filtered_df.index.get_level_values(0))
            print "Outputting %d comparisons" %(len(comparison_labels))
            for label in comparison_labels:
                print "Comparison: %s" %(label)
                comparison_output_dir = os.path.join(curr_output_dir,
                                                     label)
                if not os.path.isdir(comparison_output_dir):
                    os.makedirs(comparison_output_dir)
                output_filename = os.path.join(comparison_output_dir,
                                               "%s.%s.filtered.miso_bf" \
                                               %(label,
                                                 event_type))
                print "Outputting to: %s" %(output_filename)
                curr_df = filtered_df.ix[label].sort_index(by=sort_column,
                                                           ascending=False)
                curr_df.to_csv(output_filename,
                               sep=self.delimiter,
                               cols=columns_to_write)
#                print filtered_df.ix[label].to_csv(output_filename)
            
            
    
    def get_differential_events(self):
        """
        Get differential events.
        """
        pass


def main():
    pass


if __name__ == '__main__':
    main()
