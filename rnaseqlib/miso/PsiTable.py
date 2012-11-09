##
## Psi Table
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

import yklib
import yklib.utils as utils
import yklib.GeneTable as gt
from yklib.misowrap.miso import \
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
    return np.array(counts_vector, dtype=np.int64)
    

class PsiTable:
    """
    Class for representing Psi values.
    """
    def __init__(self, sample_labels,
                 miso_samples_dir,
                 comparisons_dir=None,
                 event_types=["A3SS",
                              "A5SS",
                              "AFE",
                              "ALE",
                              "MXE",
                              "RI",
                              "SE",
                              "SE_noAceView",
                              "TandemUTR",
                              "TandemUTR_3pseq",
                              # shortest, noAceView events
                              "A3SS_shortest_noAceView",
                              "A5SS_shortest_noAceView",
                              "MXE_shortest_noAceView",
                              "SE_shortest_noAceView"],
                              verbose=False,
                 settings_info=None,
                 gene_table=None):
        self.filtered_events = {}
        self.settings_info = settings_info
        # Load gene tables
        self.gene_table = gene_table
        self.load_gene_table()
        self.csv_output_dir = miso_samples_dir
        self.sample_labels = sample_labels
        # Directory where MISO output for samples are 
        self.miso_samples_dir = os.path.abspath(os.path.expanduser(miso_samples_dir))
        if comparisons_dir == None:
            self.comparisons_dir = os.path.join(self.miso_samples_dir,
                                                "comparisons")
        else:
            self.comparisons_dir = comparisons_dir
        self.verbose = verbose
        self.event_types = event_types
        # Summaries dataframe
        self.summaries_df = None
        # Comparisons dataframe
        self.comparisons_df = None
        # Load summaries
        self.load_summaries(self.miso_samples_dir)
        # Load comparisons
        self.load_comparisons(self.comparisons_dir)
        self.filter_coverage_events()
        # Add gene information to comparisons
        self.add_event_genes()


    def load_gene_table(self):
        """
        Try to load a gene table if one is given.
        """
        if self.settings_info == None:
            return
        if self.gene_table != None:
            # If given a table already, use it
            return
        # Load gene table based on settings info
        print "PsiTable::loading gene table..."
        self.gene_table = gt.GeneTable(self.settings_info)

        
    def add_event_genes(self, use_tables=["ensembl_genes"]):
        """
        Add the gene information to each event.
        """
        return
        # Add gene information to comparisons
        for event_type in self.event_types:
            curr_df = self.filtered_events[event_type]
            for table_source in self.gene_table.genes:
                gene_names = []
                gene_symbs = []
                if table_source not in use_tables:
                    continue
                gene_table = self.gene_table.genes[table_source]
                gene_keys = self.gene_table.table_fields[table_source]
                for event_name, event_info in curr_df.T.iteritems():
                    comparison_label, event_id = event_name
                    same_chrom_events = gene_table[event_info["chrom"] == gene_table[gene_keys["chrom"]]]
                    event_start = int(event_info["mRNA_starts"].split(",")[0])
                    event_end = int(event_info["mRNA_ends"].split(",")[0])
                    event_genes = same_chrom_events[same_chrom_events[gene_keys["strand"]] == event_info["strand"]]
                    # Check if either start or end of event lies within txStart/txEnd
                    #event_genes = event_genes[(event_start >= event_genes[gene_keys["txStart"]]) &\
                    #                           (event_end <= event_genes[gene_keys["txEnd"]])]
                    
                    first_candidates = event_genes[(event_start >= event_genes[gene_keys["txStart"]]) &\
                                                   (event_start <= event_genes[gene_keys["txEnd"]])]
                    if len(first_candidates) == 0:
                        first_candidates = event_genes[(event_end >= event_genes[gene_keys["txStart"]]) &\
                                                       (event_end <= event_genes[gene_keys["txEnd"]])]
                    event_genes = first_candidates
                    gene_name_key = gene_keys["ensGene.name2"]
                    candidate_genes = utils.unique_list(event_genes[gene_name_key])
                    num_candidates = len(candidate_genes)
                    # Event gene information
                    if num_candidates == 0:
                        print "WARNING: no gene for event %s" %(event_id)
                        event_gene_name = "NA"
                        event_gene_symbol = "NA"
                    elif num_candidates >= 1:
                        if num_candidates > 1:
                            print "WARNING: there are %d candidate genes for event %s" \
                                %(num_candidates,
                                  event_id)
                        # If there's one only gene or more, take first
                        event_gene_name = event_genes[gene_name_key].values[0]
                        event_gene_symbol = event_genes[gene_keys["geneSymbol"]].values[0]
                    gene_names.append(event_gene_name)
                    gene_symbs.append(event_gene_symbol)

                self.filtered_events[event_type][gene_keys["ensGene.name2"]] = gene_names
                self.filtered_events[event_type][gene_keys["geneSymbol"]] = gene_symbs


    def load_summaries(self, miso_samples_dir,
                       delimiter='\t'):
        """
        Load MISO summary files.
        """
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
                if self.verbose:
                    print "Processing event type %s" %(event_type)
                    print "Loading %s summary" %(sample_name)
                    print "Loading: %s" %(summary_filename)
                summary_df = pandas.read_table(summary_filename,
                                               sep=delimiter)
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
                         delimiter='\t',
                         only_two_isoform=True):
        """
        Load comparisons files.
        """
        comparisons_dict = defaultdict(list)
        dataframe_dict = {}
        for event_type in self.event_types:
            if self.verbose:
                print "Processing event type %s" %(event_type)
            event_comparisons_dir = os.path.join(comparisons_dir,
                                                 event_type)
            comparisons_dirnames = get_comparisons_dirs(event_comparisons_dir)
            comparison_labels = []
            if len(comparisons_dirnames) == 0:
                print "WARNING: No comparisons for event type %s" \
                    %(event_type)
                raise Exception
            for curr_comp_dir in comparisons_dirnames:
                comparison_label = os.path.basename(curr_comp_dir)
                comparison_labels.append(comparison_label)
                bf_filename = get_bf_filename(curr_comp_dir)
                if not os.path.isfile(bf_filename):
                    raise Exception, "No BF file: %s" %(bf_filename)
                if self.verbose:
                    print "Loading comparisons dir: %s" %(curr_comp_dir)
                    print "  - Type: %s" %(event_type)
                    print "  - Comparison: %s" %(comparison_label)
                curr_df = pandas.read_table(bf_filename,
                                            sep=delimiter,
                                            index_col=[0])
                # Keep only events for which we have two isoforms
                if only_two_isoform:
                    if event_type == "TandemUTR_3pseq":
                        curr_df = self.filter_only_two_isoform(curr_df)
                comparisons_dict[event_type].append(curr_df)
            # Concatenate all the DataFrames for each comparison together
            dataframe_dict[event_type] = pandas.concat(comparisons_dict[event_type],
                                                       keys=comparison_labels)
            # Add gene information
#        self.comparisons_df = pandas.DataFrame(dataframe_dict)
        self.comparisons_df = dataframe_dict



    def load_comparisons_counts_from_df(self, df,
                                        counts_labels=["sample1_counts",
                                                       "sample2_counts"]):
        """
        Return sample1 and sample2 counts from comparisons
        MISO file.
        """
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
                               atleast_sum_inc_exc=20):
        """
        Filter events for coverage.
        """
        print "filter_coverage_events::Filtering..."
        if comparisons_df == None:
            comparisons_df = self.comparisons_df
        for event_type in self.event_types:
            print "Filtering event type: %s" %(event_type)
            comparison_counts = self.load_comparisons_counts_from_df(comparisons_df[event_type])
            # Get counts for each read class for sample 1 and sample 2
            comparison_counts = self.get_counts_by_class("sample1_counts_int", "sample1",
                                                         comparison_counts)
            comparison_counts = self.get_counts_by_class("sample2_counts_int", "sample2",
                                                         comparison_counts)
            filtered_df = comparison_counts
            # Filter exclusion reads
            # Only apply this to events other than TandemUTRs!
            if "TandemUTR" in event_type:
                # For tandem UTRs, apply a filter on the (1,1) reads instead
                filtered_df = filtered_df[filtered_df["sample1_const_counts"] | filtered_df["sample2_const_counts"] >= atleast_exc]
            elif "AFE" == event_type:
                # Use more aggressive filtering for AFEs
                AFE_atleast_inc = 10
                AFE_atleast_exc = 10
                # Filter inclusion reads
                filtered_df = filtered_df[filtered_df["sample1_inc_counts"] | filtered_df["sample2_inc_counts"] >= AFE_atleast_inc]
                filtered_df = filtered_df[filtered_df["sample1_exc_counts"] | filtered_df["sample2_exc_counts"] >= AFE_atleast_exc]
            else:
                # Filter inclusion reads
                filtered_df = filtered_df[filtered_df["sample1_inc_counts"] | filtered_df["sample2_inc_counts"] >= atleast_inc]
                filtered_df = filtered_df[filtered_df["sample1_exc_counts"] | filtered_df["sample2_exc_counts"] >= atleast_exc]
                
            # Filter the sum of inclusion and exclusion
            sample1_sum = filtered_df["sample1_inc_counts"] + filtered_df["sample1_exc_counts"]
            sample2_sum = filtered_df["sample2_inc_counts"] + filtered_df["sample2_exc_counts"]
            filtered_df = filtered_df[sample1_sum | sample2_sum >= atleast_sum_inc_exc]
            self.filtered_events[event_type] = filtered_df

        
    def get_counts_by_class(self, col_label, df_col, df):
        """
        Return counts for each MISO read class.
        """
        df["%s_inc_counts" %(df_col)] = np.array(map(lambda x: x[0], df[col_label].values))
        df["%s_exc_counts" %(df_col)] = np.array(map(lambda x: x[1], df[col_label].values))
        df["%s_const_counts" %(df_col)] = np.array(map(lambda x: x[2], df[col_label].values))
        df["%s_neither_counts" %(df_col)] = np.array(map(lambda x: x[3], df[col_label].values))
        return df


    def output_filtered_comparisons(self, output_dir=None,
                                    sort_column="bayes_factor",
                                    delimiter='\t',
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
            output_dir = self.csv_output_dir
        # Output each file by event type
        output_dir = os.path.join(output_dir, "filtered_events")
        print "output_filtered_comparisons::writing to dir: %s" %(output_dir)
        for event_type, filtered_df in self.filtered_events.iteritems():
            curr_output_dir = os.path.join(output_dir, event_type)
            print "Event type: %s" %(event_type)
            # View by comparison
            comparison_labels = utils.unique_list(filtered_df.index.get_level_values(0))
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
                curr_df.to_csv(output_filename, sep=delimiter,
                               cols=columns_to_write)
#                print filtered_df.ix[label].to_csv(output_filename)
            
            
    
    def get_differential_events(self):
        """
        Get differential events.
        """
        pass


def main():
    sample_labels = [["KH2_NoDox_A", "KH2 -Dox (A)"],
                     ["MSI1_NoDox_A", "MSI1 -Dox (A)"],
                     ["MSI1_NoDox_B", "MSI1 -Dox (B)"],
                     ["MSI1_DOX_A", "MSI1 +DOX (A)"]]
    miso_samples_dir = os.path.expanduser("~/jaen/Musashi/rna-seq/miso_output/")
    psi_table = PsiTable(sample_labels,
                         miso_samples_dir)


if __name__ == '__main__':
    main()
