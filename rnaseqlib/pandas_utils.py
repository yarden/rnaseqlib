##
## Pandas related utilities
##
import os
import sys
import time

from numpy import *

import pandas

def select_df_rows(df, cond,
                   columns=None,
                   how='any'):
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
    
