##
## Statistics-related utilities
##
import scipy
import numpy
from numpy import *
from scipy import *

from scipy.stats.stats import pearsonr, spearmanr

def spearman_dist(u, v, na_vals=["NA"]):
    """
    Compuute Spearman distance for vectors u, v.

    Returns 1 - spearmanr(u, v).
    """
    matrix = [[x, y] for x, y in zip(u, v) \
              if (u not in na_vals) and (v not in na_vals)]
    matrix = array(matrix)
    spearman = scipy.stats.spearmanr(matrix[:, 0], matrix[:, 1])[0]
    return 1 - spearman


def my_pdist(X, dist_func,
             na_values=["NA"]):
    """
    (inefficient) pdist function that takes an arbitrary
    distance function (a lambda).

    X: data matrix
    dist_func: lambda that returns distance on vector
    na_values: values to consider as missing data
    """
    X = array(X, dtype=object)
    num_rows, num_cols = X.shape
    dist_matrix = []
    for col1 in range(num_cols):
        pdist_row = []
        for col2 in range(num_cols):
            pairs = array([[x, y] for x, y in zip(X[:, col1], X[:, col2]) \
                           if (x not in na_values) and (y not in na_values)])
            if len(pairs) == 0:
                continue
            dist = dist_func(pairs[:, 0],
                             pairs[:, 1])
            pdist_row.append(dist)
        dist_matrix.append(pdist_row)
    dist_matrix = array(dist_matrix)
    return dist_matrix

