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


def combine_dfs(dfs_list, **kwargs):
    """
    Combine dataframes into one. 
    """
    if len(dfs_list) == 1:
        combined_df = dfs_list[0]
        return combined_df
    combined_df = reduce(lambda x, y: x.combine_first(y), dfs_list)
#    combined_df = dfs_list[0]
#    for right_df in dfs_list[1:]:
        # Get the new columns that are not common
        #new_noncommon_cols = [c for c in right_df \
        #                      if (c not in combined_df.columns) or \
        #                         (c in kwargs["on"])]
#        new_noncommon_cols = [c for c in right_df \
#                              if (c not in combined_df.columns)]
#        print "COMBINING USING: ", new_noncommon_cols
#        combined_df = pandas.merge(combined_df,
#                                   right_df[new_noncommon_cols],
#                                   left_index=True,
#                                   right_index=True,
#                                   **kwargs)
#                                   left_index=True,
#                                   right_index=True,
#                                   suffixes=["", ""],
#                                   **kwargs)
    # combined_df = reduce(lambda left_df, right_df:
    #                      pandas.merge(left_df, right_df,
    #                                   left_index=True,
    #                                   right_index=True,
    #                                   suffixes=["_%s" %(df_labels_dict[left_df]),
    #                                             "_%s" %(df_labels_dict[right_df])],
    #                                   **kwargs),
    #                      dfs_list)
    return combined_df


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
    
