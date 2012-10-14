##
## Expression table
##
import os
import time
import glob

import misopy
import pandas

import numpy as np

import yklib

def compute_fold_changes_table(table,
                               pairs_to_compare,
                               c=0.01):
    """
    Return matrix of pairs to compare fold changes.
    """
    fold_changes = []
    for pair in pairs_to_compare:
        label1, label2 = pair
        sample1 = table[label1].values + c
        sample2 = table[label2].values + c
        pair_fc = sample1 / sample2
        fold_changes.append(list(pair_fc))
    fold_changes = np.array(fold_changes).T
    return fold_changes

    
class ExpressionTable:
    """
    Class for representing gene expression values.
    """
    def __init__(self, label=None, header_fields=None,
                 from_file=None, delimiter='\t',
                 counts_dir=None, index_col=None):
        self.label = label
        self.header_fields = header_fields
        self.from_file = from_file
        self.delimiter = delimiter
        self.index_col = index_col
        self.counts_dir = counts_dir
        self.counts_panel = {}

        # Label to use when plotting
        self.plot_label = label

        # data frame
        self.data = None
        
        if self.from_file != None:
            self.load_rpkm_table(from_file)

            # Load counts information if given
            if self.counts_dir != None:
                self.load_counts_dir(counts_dir)

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
                          thresholds={'rpkm': {'mode': 'any',
                                               'cutoff': 0.3},
                                      'counts': {'mode': 'any',
                                                 'cutoff': 5}}):
        column_names = list(self.data.columns.values)
        # Get columns corresponding to selected labels
        selected_cols = [column_names.index(label) \
                         for label in sample_labels]

        filtered_data = None
        
        # Slice only relevant columns and apply RPKM cutoff
        rpkm_cutoff = thresholds['rpkm']['cutoff']
        filtered_rpkm_index = self.data.ix[:, selected_cols].values > rpkm_cutoff
        if thresholds['rpkm']['mode'] == 'any':
            filtered_data = self.data[filtered_rpkm_index.any(1)]
        elif thresholds['rpkm']['mode'] == 'all':
            filtered_data = self.data[filtered_rpkm_index.all(1)]

        # Apply counts cutoff
        # Get dataframe containing only counts column for each sample
        counts_by_samples = self.counts_panel.minor_xs("counts")
        counts_threshold = thresholds["counts"]["cutoff"]

        # Apply threshold disjunctively in any sample
        counts_met = None
        if thresholds["counts"]["mode"] == "any":
            counts_met = reduce(lambda x, y: x | y,
                                [col >= counts_threshold \
                                 for _, col in counts_by_samples.iteritems()])
        filtered_data = filtered_data.reindex(counts_met.index).dropna()
        return filtered_data
        
                        
    # def filter_rpkm_table(self, rpkm_cutoff=0.3,
    #                       cutoff_mode='any',
    #                       min_read_count=None,
    #                       sample_rpkm_files=None,
    #                       entry_id_header='#Gene'):
    #     """
    #     Filter RPKM table.
    #     """
    #     filtered_rpkm_table = []
    #     genes_to_sample_counts = {}
        
    #     if min_read_count != None:
    #         if sample_rpkm_files == None:
    #             raise Exception, "Sample RPKM files must be given to " \
    #                   "use read count cutoff."
    #         genes_to_sample_counts = load_rpkm_read_counts(sample_rpkm_files,
    #                                                        column_labels_to_use)

    #     rpkm_table = self.table_dictlist

    #     for gene in rpkm_table:
    #         # Convert fields
    #         gene = evalDict(gene)

    #         # Filter based on read counts
    #         counts_profile = None
    #         if min_read_count != None:
    #             gene_id = gene[entry_id_header]
    #             counts_profile = array(genes_to_sample_counts[gene_id].values())

    #             if cutoff_mode == 'any':
    #                 if not (any(counts_profile >= min_read_count)):
    #                     continue
    #             elif cutoff_mode == 'all':
    #                 if not (all(counts_profile >= min_read_count)):
    #                     continue

    #         # Filter based on expression
    #         expression_profile = array([gene[label] \
    #                                     for label in column_labels_to_use])

    #         # enforce all samples to have RPKM >= 0.3
    #         if any(expression_profile < rpkm_cutoff):
    #             continue

    #         if cutoff_mode == 'any':
    #             if any(expression_profile >= rpkm_cutoff):
    #                 filtered_rpkm_table.append(gene)
    #         elif cutoff_mode == 'all':
    #             if all(expression_profile >= rpkm_cutoff):
    #                 filtered_rpkm_table.append(gene)

    #     print "Filtering resulted in total of %d genes." \
    #         %(len(filtered_rpkm_table))
    #     self.filtered_table_dictlist = filtered_rpkm_table
        

    def load_rpkm_table(self, rpkm_filename, delimiter=None,
                        index_col=None):
        """
        Load RPKM table.
        """
        if delimiter == None:
            delimiter = self.delimiter

        if index_col == None:
            index_col = self.index_col
            
        print "Loading RPKM table from: %s" %(rpkm_filename)

        # Load table
        self.data = pandas.read_table(rpkm_filename,
                                      sep=delimiter)
        if index_col is not None:
            self.data = self.data.set_index(index_col)

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