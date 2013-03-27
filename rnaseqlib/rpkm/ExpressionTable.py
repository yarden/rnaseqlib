##
## Expression table
##
import os
import time
import glob
import pandas
pandas.set_option('use_inf_as_null', True)

import scipy
import numpy as np

from scipy import *

import misopy

from scipy.stats.stats import zscore

import rnaseqlib
import rnaseqlib.ribo.ribo_utils as ribo_utils


def compute_fold_changes_table(table,
                               pairs_to_compare,
                               c=0.01,
                               logged=False):
    """
    Return matrix of pairs to compare fold changes.

    - logged: if True, assume values are logged, so compute
      fold change by subtraction.
    """
    fold_changes = []
    for pair in pairs_to_compare:
        label1, label2 = pair
        sample1 = table[label1].values + c
        sample2 = table[label2].values + c
        if logged:
#            pair_fc = sample1 - sample2
            pair_fc = 2**(sample1) / 2**(sample2)
        else:
            pair_fc = sample1 / sample2
        fold_changes.append(list(pair_fc))
    fold_changes = np.array(fold_changes).T
    return fold_changes


def get_ribo_rpkm_table(rpkm_filename,
                        ribo_sample_pairs,
                        rna_expr_table=None,
                        rna_to_ribo_samples={}):
    """
    Get ribo RPKM table with fold changes.
    """
    expr_table = RiboExpressionTable(from_file=rpkm_filename,
                                     comparison_pairs=ribo_sample_pairs,
                                     rna_expr_table=rna_expr_table,
                                     rna_to_ribo_samples=rna_to_ribo_samples)
    # Add fold changes
    expr_table.add_fold_changes(ribo_sample_pairs)
    return expr_table



class RiboExpressionTable:
    """
    Class for representing expression values
    from Ribo-Seq.
    """
    def __init__(self, from_file=None,
                 label=None,
                 delimiter="\t",
                 index_col="gene_id",
                 comparison_pairs=None,
                 na_val="NA",
                 ribo_comparison_pairs=[],
                 rna_expr_table=None,
                 rna_to_ribo_samples={},
                 rna_comparison_pairs=[]):
        self.from_file = from_file
        self.label = label
        self.table = None
        self.indexed_table = None
        self.delimiter = delimiter
        self.ribo_comparison_pairs = ribo_comparison_pairs
        self.rna_expr_table = rna_expr_table
        self.rna_to_ribo_samples = rna_to_ribo_samples
        self.rna_comparison_pairs = rna_comparison_pairs
        # Column to index by
        self.index_col = index_col
        # Samples to compare
        self.comparison_pairs = comparison_pairs
        # NA value to use
        self.na_val = na_val

        if self.from_file is not None:
            self.load_table(from_file)

            
    def load_table(self, table_filename):
        print "Loading table from: %s" %(table_filename)
        self.table = pandas.read_table(table_filename,
                                       sep=self.delimiter)
        # Index the table
        self.indexed_table = self.table.set_index(self.index_col)
        # Add RNA expression to table if available
        if self.rna_expr_table is not None:
            self.add_rna_expr_table(self.rna_expr_table)
            # Compute TEs
            self.add_te(self.rna_to_ribo_samples)


    def add_rna_expr_table(self, rna_expr_table):
        """
        Add RNA expression table to ribo table.
        """
        # Merge the tables using the left index
        merged_rpkm_table = pandas.merge(self.table,
                                         rna_expr_table.data,
                                         left_on="gene_id",
                                         right_on="#Gene",
                                         # Merge on left index, the index of
                                         # the Ribo RPKM
                                         how="left")
        self.table = merged_rpkm_table
        return self.table
    

    def get_counts_columns(self):
        """
        Return the counts columns.
        """
        counts_cols = []
        for col in self.table.columns:
            if col.startswith("counts_"):
                counts_cols.append(col)
        return counts_cols


    def get_rpkm_columns(self):
        rpkm_cols = []
        for col in self.table.columns:
            if col.startswith("rpkm_"):
                rpkm_cols.append(col)
        return rpkm_cols


    def add_te(self, rna_to_ribo_samples):
        """
        Add TE calculations to the given RPKM table
        that contains merged RPKM and Ribo values.

        - rna_to_ribo_samples: mapping from RNA sample names
          to ribo sample names
        Move me later to RiboExpressionTable class.
        """
        for rna_sample, ribo_sample in rna_to_ribo_samples.iteritems():
            te_label = "TE_%s_%s" %(ribo_sample, rna_sample)
            te = ribo_utils.compute_te(self.table[ribo_sample],
                                       self.table[rna_sample])
            # Add TE values to table
            self.table[te_label] = te
        # Normalize TE
        self.add_normalized_te()
        return self.table

            
    def add_normalized_te(self, normed_prefix="norm"):
        """
        z-score normalize the TE values.

        Creates new columns corresponding to normed TEs,
        beginning with 'normed_prefix'.

        Normed TEs are first logged (base 2) and then z-score
        normalized.
        """
        print "Normalizing TE..."
        te_cols = [c for c in self.table.columns \
                   if c.startswith("TE_")]
        for col in te_cols:
            normed_col = "%s_%s" %(normed_prefix, col)
            self.table[normed_col] = zscore(self.table[col].apply(log2).dropna())
            

    def filter_table(self, table,
                     thresholds={'rpkm': {'mode': 'any',
                                          'cutoff': 0},
                                 'counts': {'mode': 'all',
                                            'cutoff': 5}}):
        """
        Filter table.

        Return filtered table.
        """
        ##
        ## TODO: Add option to only filter on a subset of
        ## samples
        ##
        counts_cols = self.get_counts_columns()
        rpkm_cols = self.get_rpkm_columns()
        filtered_data = None
        # Apply RPKM cutoff: slice only relevant rows
        # out of the RPKM index
        rpkm_cutoff = thresholds["rpkm"]["cutoff"]
        filtered_rpkm_index = self.table.ix[:, rpkm_cols].values > rpkm_cutoff
        filtered_rpkm = None
        if thresholds["rpkm"]["mode"] == "any":
            filtered_rpkm = self.table[filtered_rpkm_index.any(1)]
        elif thresholds["rpkm"]["mode"] == "all":
            filtered_rpkm = self.table[filtered_rpkm_index.all(1)]
        # Apply counts cutoff
        # Get dataframe containing only counts column for each sample
        counts_cutoff = thresholds["counts"]["cutoff"]
        filtered_counts_index = filtered_rpkm.ix[:, counts_cols].values > counts_cutoff
        if thresholds["counts"]["mode"] == "any":
            filtered_data = filtered_rpkm[filtered_counts_index.any(1)]
        elif thresholds["counts"]["mode"] == "all":
            filtered_data = filtered_rpkm[filtered_counts_index.all(1)]
        return filtered_data
        


    def get_rpkm_sample_names(self, header, prefix="rpkm_"):
        """
        Return any sample beginning with RPKM prefix.
        """
        sample_names = []
        for col in header:
            if col.startswith(prefix):
                sample_names.append(col)
        return sample_names


    def get_expr_matrix(self, samples, table):
        """
        Get expression as matrix for the given samples.
        """
        expr_matrix = []
        for sample in samples:
            expr_matrix.append(table[sample].values)
        expr_matrix = np.array(expr_matrix)
        return expr_matrix
        

    def add_fold_changes(self, sample_pairs,
                         logged=False):
        """
        Add fold changes to the table in the form of
        X_vs_Y fields, where X, Y are samples.

        Takes a list of sample pairs.
        """
        if self.table is None:
            raise Exception, "No data to compute foldchanges."
        fold_changes = compute_fold_changes_table(self.table,
                                                  sample_pairs,
                                                  logged=logged)
        for pair_num, pair in enumerate(sample_pairs):
            label = "%s_vs_%s" %(pair[0],
                                 pair[1])
            fc_column = fold_changes[:, pair_num]
            self.table[label] = fc_column
        return self.table


    # def loess_normalize(self, sample_pairs):
    #     """
    #     Apply loess normalization to filtered table.
    #     """
    #     ##
    #     ## Add fold changes with loess normalization
    #     ##
    #     print "Loess normalizing..."
    #     normalized_table = loess_normalize_table(self.filtered_table, sample_pairs)
    #     if normalized_table is not None:
    #         # Normalization was successful
    #         self.filtered_table = normalized_table
    #     return normalized_table

        
    def __repr__(self):
        return "RiboExpressionTable(from_file=%s)" %(self.from_file)

    
class ExpressionTable:
    """
    Class for representing gene expression values.
    """
    def __init__(self, label=None,
                 header_fields=None,
                 from_file=None,
                 delimiter='\t',
                 counts_dir=None,
                 index_col=None):
        self.label = label
        self.header_fields = header_fields
        self.from_file = from_file
        self.delimiter = delimiter
        self.index_col = index_col
        self.counts_dir = counts_dir
        self.counts_panel = {}
        self.indexed_data = None
        # Label to use when plotting
        self.plot_label = label

        # data frame
        self.data = None
        
        if self.from_file != None:
            self.load_rpkm_table(from_file)

            # Load counts information if given
            if self.counts_dir != None:
                self.load_counts_dir(counts_dir)
        self.index_data()


    def index_data(self):
        if (self.data is not None) and (self.index_col is not None):
            self.indexed_data = self.data.set_index(self.index_col)
            
                
    def add_fold_changes(self, sample_pairs):
        """
        Add fold changes to the table in the form of
        X_vs_Y fields, where X, Y are samples.

        Takes a list of sample pairs.
        """
        if self.data is None:
            raise Exception, "No data to compute foldchanges."
        fold_changes = compute_fold_changes_table(self.data,
                                                  sample_pairs)
        for pair_num, pair in enumerate(sample_pairs):
            label = "%s_vs_%s" %(pair[0],
                                 pair[1])
            fc_column = fold_changes[:, pair_num]
            self.data[label] = fc_column
        return self.data

        
    def filter_rpkm_table(self,
                          sample_labels,
                          thresholds={'rpkm': [{'mode': 'any',
                                                'cutoff': 0.3}],
                                      'counts': [{'mode': 'any',
                                                  'cutoff': 5}]}):
        """
        Filter RNA RPKM table.

        - thresholds: Describes a set of filters to apply to
          RPKM values and or count values.
        """
        self.index_data()
        column_names = list(self.data.columns.values)
        # Get columns corresponding to selected labels
        selected_cols = [column_names.index(label) \
                         for label in sample_labels]
        filtered_data = None
        ##
        ## Apply RPKM value filters
        ##
        for rpkm_filter in thresholds["rpkm"]:
            # Slice only relevant columns and apply RPKM cutoff
            rpkm_cutoff = rpkm_filter["cutoff"]
            filtered_rpkm_index = \
                self.indexed_data.ix[:, selected_cols].values > rpkm_cutoff
            if rpkm_filter["mode"] == "any":
                filtered_data = self.indexed_data[filtered_rpkm_index.any(1)]
            elif rpkm_filter["mode"] == "all":
                filtered_data = self.indexed_data[filtered_rpkm_index.all(1)]
        ##
        ## Apply counts filters
        ##
        # Get dataframe containing only counts column for each sample
        counts_by_samples = self.counts_panel.minor_xs("counts")
        for counts_filter in thresholds["counts"]:
            counts_threshold = counts_filter["cutoff"]
            # Apply threshold disjunctively in any sample
            counts_met = None
            if counts_filter["mode"] == "any":
                counts_met = \
                    reduce(lambda x, y: x | y,
                           [col >= counts_threshold \
                            for _, col in counts_by_samples.iteritems()])
            else:
                raise Exception, "Not implemented mode %s" \
                    %(counts_filter["mode"])
            filtered_data = filtered_data.reindex(counts_met.index).dropna()
        return filtered_data

        
    # def loess_normalize(self, sample_pairs):
    #     """
    #     Apply loess normalization to filtered table.
    #     """
    #     ##
    #     ## Add fold changes with loess normalization
    #     ##
    #     print "Loess normalizing..."
    #     normalized_table = loess_normalize_table(self.filtered_data, sample_pairs)
    #     if normalized_table is not None:
    #         # Normalization was successful
    #         self.filtered_data = normalized_table
    #     return normalized_table
        

    def load_rpkm_table(self, rpkm_filename, delimiter=None):
        """
        Load RPKM table.
        """
        if delimiter == None:
            delimiter = self.delimiter
        print "Loading RPKM table from: %s" %(rpkm_filename)
        # Load table
        self.data = pandas.read_table(rpkm_filename,
                                      sep=delimiter)

    def load_counts_dir(self, counts_dir,
                        counts_index=['#Gene',
                                      'counts',
                                      'length',
                                      'rpkm']):
        """
        Load RPKM counts directories.
        """
        counts_dir = os.path.abspath(os.path.expanduser(counts_dir))
        
        if not os.path.isdir(counts_dir):
            raise Exception, "%s not a counts dir." \
                  %(counts_dir)

        counts = {}
        for fname in glob.glob(os.path.join(counts_dir, "*.rpkm")):
            table_name = os.path.basename(fname).split(".")[0]
            counts_filename = os.path.join(counts_dir, fname)
            print "Loading: %s" %(counts_filename)
            print "  - Sample name: %s" %(table_name)
            counts[table_name] = pandas.read_table(counts_filename,
                                                   sep=self.delimiter,
                                                   names=counts_index)
            counts[table_name] = counts[table_name].set_index("#Gene")
            
        # Cast to DataFrame
        self.counts_panel = pandas.Panel(counts)
