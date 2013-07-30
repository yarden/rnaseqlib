##
## Pandas related utilities
##
import os
import sys
import time

from numpy import *
import numpy as np

import pandas

import rnaseqlib
import rnaseqlib.utils as utils



def get_abs_cov_df(df):
    """
    Return the absolute coefficient of variation of the rows
    in the given DataFrame with respect to columns.
    """
    # Compute absolute coefficient of variation for each entry in dataframe
    abs_cov_df = df.std(axis=1).div(df.mean(axis=1)).abs()
    return abs_cov_df


def get_variable_rows(df, top_percentile=75):
    """
    Return subset of dataframe that contains the most variable rows.
    Most variable rows are defined as rows whose absolute coefficient of
    variation, abs(sigma / mu), across columns is in the
    'top_percentile' of deviations for the dataframe.
    """
    abs_cov_df = get_abs_cov_df(df)
    # Normalized df should have same number of rows as original
    num_rows = abs_cov_df.shape[0]
    if df.ndim == 1:
        orig_num_rows = df.values.shape
    else:
        orig_num_rows, _ = df.values.shape
    assert (num_rows == orig_num_rows), \
        "zscore normed df should have same number of rows as original " \
        "dataframe."
    percentile_cutoff = np.percentile(abs_cov_df.values, top_percentile)
    # Get subset of original dataframe containing only the most
    # variable columns
    var_df = df[abs_cov_df >= percentile_cutoff]
    return var_df


def zscore_norm_df(df, how):
    """
    Z-score normalize a DataFrame by either column
    or row
    """
    norm_df = df.copy()
    if how == "column":
        # Take mean across the rows (axis=0)
        mu_val = df.mean(axis=0)
        # Take standard deviation across the rows (axis=0)
        std_val = df.std(axis=0)
        # Subtract mean across columns (axis=1) and divide
        # by standard deviation
        norm_df = df.sub(mu_val, axis=1).sub(std_val, axis=1)
    elif how == "row":
        raise Exception, "Not implemented yet."
    else:
        raise Exception, "Don\'t know how to zscore normalize by: %s" %(how)


def replace_inf_with_na(df, val_to_use=np.nan):
    """
    Replace infinite values with np.nan
    """
    df = df.replace([np.inf, -np.inf], np.nan)
    return df


def common_cols(df1, df2, except_cols=[]):
    """
    Get common cols between df1/df2 except
    the ones given in except cols.
    """
    common = set(df1.columns).intersection(set(df2.columns))
    common = common - set(except_cols)
    return list(common)


# def shuffle_df(df, n, axis=0):
#     """
#     Shuffle dataframe. If axis is 0, shuffle by rows (each column
#     shuffled independently), if 1 shuffle by columns.

#     n is number of shuffles.
#     """
#     shuffled_df = df.copy()
#     for k in range(n):
#         if axis == 0:
#             for c in shuffled_df.columns:
#                 # Shuffle the current column
#                 curr_col = shuffled_df[c]
#                 new_index = np.arange(len(curr_col))
#                 print "index before: ", new_index
#                 np.random.shuffle(new_index)
#                 print "index after: ", new_index
#                 # Assign this shuffled column as new column
#                 shuffled_df[c] = curr_col[new_index]
#         else:
#             raise Exception, "Not implemented."
#     if axis == 0:
#         shuffled_df.columns = df.columns
#     return shuffled_df

def shuffle_df(df, n, axis=0):
    """
    Shuffle dataframe rows (independently) if axis == 0.
    If axis == 1, shuffle columns independently.
    """
    shuffled_df = df.copy()
    for k in range(n):
        shuffled_df.apply(np.random.shuffle, axis=axis)
    return shuffled_df

# def shuffle_df(df, n, axis=0):
#     shuffled_df = df.copy()
#     for k in range(n):
#         if axis == 0:
#             np.random.shuffle(shuffled_df.values)
#         elif axis == 1:
#             np.random.shuffle(shuffled_df.T.values)
#         else:
#             raise Exception, "Not implemented."
#     return shuffled_df


def merge_dfs(dfs_list, how="outer"):
    if len(dfs_list) == 1:
        return dfs_list[0]
    merged_df = dfs_list[0]
    for other_df in dfs_list[1:]:
        common_cols = \
            set(list(merged_df.columns)).intersection(set(list(other_df.columns)))
        common_cols = list(common_cols)
        merged_df = merged_df.merge(other_df,
                                    on=common_cols,
                                    how=how)
    return merged_df
              

def is_df_indexed(df):
    """
    Check if dataframe has the default index or not.
    """
    num_rows = len(df)
    indexed = not pandas.Index(np.arange(0, num_rows)).equals(df.index)
    return indexed


def combine_dfs(dfs_list, axis=1, combine_first=True):
    """
    Concatenate dataframe columns together, non-redundantly.

    Assumes the concatenation occurs based on the index of
    each df (each df needs to be indexed.)

    If 'combine_first' is True, first fill in all the non-missing
    values between dataframes on common columns
    """
    first_df = dfs_list[0]
    if len(dfs_list) == 1:
        return first_df
    # Check that the dataframes are indexed
    for df in dfs_list:
        assert (is_df_indexed(df) == True), \
               "Dataframe %s is not indexed! Cannot combine." \
               %(str(df.columns))
    # If asked, combine dataframes first to fill in non-missing
    # values where there are missing values. IMPORTANT: when combining
    # first, consider only columns that are *common*
    if combine_first:
        for next_df in dfs_list:
            common_cols = first_df.columns.intersection(next_df.columns)
            first_df = first_df.combine_first(next_df[common_cols])
            # Check for redundancy
            assert (len(first_df.columns) == len(set(first_df.columns))), \
                   "Redundancy found in combine_first step!"
    nonredundant_dfs = [first_df]
    cols = list(first_df.columns)
    for df in dfs_list:
        # Only consider new columns 
        new_cols = list(df.columns.diff(cols))
        if len(new_cols) == 0:
            continue
        nonredundant_dfs.append(df[new_cols])
        cols = new_cols + list(cols)
        print "Checking redunancy of combined columns"
        assert (len(cols) == len(set(cols))), \
               "Combined columns are redundant!"
    new_df = pandas.concat(nonredundant_dfs, axis=axis)
    assert (len(new_df.columns) == len(set(new_df.columns))), \
           "Combined columns are redundant at end!"
    return new_df
    

def select_df_rows(df, cond,
                   columns=None,
                   how="any"):
    """
    Select rows from DataFrame where a condition
    holds.  cond is a lambda.
    """
    if columns is None:
        columns = df.columns
    # Apply condition result and get array as result
    cond_result = df[columns].applymap(cond).values
    result = None
    if how == 'any':
        result_ind = cond_result.any(1)
    elif how == 'all':
        result_ind = cond_result.all(1)
    elif type(how) == int:
        # apply to N many columns
        result_ind = cond_result.sum(axis=1) >= how
    else:
        raise Exception, "Do not know %s" %(how)
    result = df[result_ind]
    return result, result_ind


if __name__ == "__main__":
    from numpy import *
    from scipy import *
    df = pandas.DataFrame([[0.5, 0.1, 0.7],
                           [0.6, 0.8, 0.6],
                           [0.505, 0.02, 0.9]])
    print "df: "
    print df
    # Select rows where all cols above 0.5
#    any_result = select_df_rows(df, lambda x: x >= 0.5,
#                                how='any')
#    print "any:"
#    print any_result
    # Select rows where any cols above 0.5
#    all_result = select_df_rows(df, lambda x: x >= 0.5,
#                                how='all')
#    print "all:"
#    print all_result
    # Select rows where 2 cols above 0.5
    two_result = select_df_rows(df, lambda x: x > 0.5,
                                how=2)
    print two_result
    
