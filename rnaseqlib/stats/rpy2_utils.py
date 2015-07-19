##
## Utilities related to Rpy2 library
##
import os
import sys
import time

import pandas
import numpy as np

from collections import OrderedDict

import rnaseqlib
import rnaseqlib.utils as utils

try:
    import rpy2
    from rpy2.robjects import r
    import rpy2.robjects as robj
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr
    rpy2.robjects.numpy2ri.activate()
    from rpy2.robjects.lib import grid
    from rpy2.robjects import r, Formula
    py2ri_orig = rpy2.robjects.conversion.py2ri
except:
    print "WARNING: Cannot import rpy2"

def series2ndarray(s):
    if isinstance(s, pandas.core.series.Series):
        return s.values
    else:
        return s

py2ri_orig = rpy2.robjects.conversion.py2ri
    
def conversion_pydataframe(obj):
    if isinstance(obj, pandas.core.frame.DataFrame):
        od = OrderedDict()
        for name, values in obj.iteritems():
            if values.dtype.kind == 'O':
                od[name] = rpy2.robjects.vectors.StrVector(series2ndarray(values))
            else:
                od[name] = rpy2.robjects.conversion.py2ri(series2ndarray(values))
        return rpy2.robjects.vectors.DataFrame(od)
    else:
        return py2ri_orig(obj)

rpy2.robjects.conversion.py2ri = conversion_pydataframe
df_to_r = conversion_pydataframe

def r_df_get_col(r_df, col_name):
    """
    Return the column named 'col_name' from the R DataFrame
    object.
    """
    col_names = [curr_col for curr_col in r_df.colnames]
    if col_name not in col_names:
        raise Exception, "No column %s found in R dataframe." %(col_name)
    # Add 1 because R is 0-based (ugh)
    col_index = col_names.index(col_name) + 1
    return r_df.rx(col_index)


def run_ma_loess(x, y, span=0.75):
    """
    Run MA-based loess normalization on X and Y. Computes

      M = log(X/Y)
      A = 0.5 * log(X*Y)

    Fits loess regression M ~ A and corrects X and Y accordingly.

    Assumes input X and Y values are non-logged.
    """
    M = np.log2(x) - np.log2(y)
    # A = average intensity 1/2(XY)
    A = 0.5 * (np.log2(x) + np.log2(y))
    # Fit M ~ A
    corrected_m, correction_factor = run_loess(A, M, span=span)
    corrected_x = 2**((2*A + corrected_m)/2.)
    corrected_y = 2**((2*A - corrected_m)/2.)
    return corrected_x, corrected_y


def where_na_like(l):
    """
    Return indices where array is NA-like
    """
    bool_index = np.array(map(lambda x: np.isinf(x) or \
                              pandas.isnull(x), l))
    return np.where(bool_index)[0]


def run_loess(x, y, span=0.75):
    """
    Predict y as function of x. Takes two numpy vectors.
    """
    # Ensure that Inf/-Inf values are substituted
    x[where_na_like(x)] = robj.NA_Real
    y[where_na_like(x)] = robj.NA_Real
    data = robj.DataFrame({"x": x, "y": y})
    loess_fit = r.loess("y ~ x", data=data, span=span,
                        family="symmetric")
    correction_factor = np.array(list(r.predict(loess_fit, x)))
    corrected_y = \
        np.array(list(y)) - correction_factor
    return corrected_y, correction_factor


def run_lowess(x, y, span=0.75):
    """
    Predict y as function of x. Takes two numpy vectors.

    Uses 'r.lowess' instead of 'r.loess'.
    """
    # Ensure that Inf/-Inf values are substituted
    x[where_na_like(x)] = robj.NA_Real
    y[where_na_like(x)] = robj.NA_Real
    data = robj.DataFrame({"x": x, "y": y})
    print "x: ", x, "y: ", y
    lowess_fit = r.lowess(data, f=span)
    print "LOWESS FIT: ", lowess_fit
    corrected_y = np.array(list(lowess_fit.y))
    return corrected_y, correction_factor



def main():
    pass
    #x = np.array([1, 0.5, 3, 4, 5, 5.5, 6, 7], dtype=np.float)
    #y = np.array([10, 25, 38, 44.5, 500, 550, 600, 705], dtype=np.float)
    #print "Results loess: ", run_loess(x, y)
    #print "Results LOWESS: ", run_lowess(x, y)
    


if __name__ == "__main__":
    main()
