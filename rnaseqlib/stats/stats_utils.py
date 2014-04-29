##
## Statistics-related utilities
##
import scipy
import numpy
import numpy as np
from numpy import *
from scipy import *

from scipy.stats.stats import pearsonr, spearmanr


def spearman_dist(u, v, na_vals=["NA", np.nan]):
    """
    Compute Spearman distance for vectors u, v.

    Returns 1 - spearmanr(u, v).
    """
    matrix = [[x, y] for x, y in zip(u, v) \
              if (x not in na_vals) and (y not in na_vals)]
    matrix = array(matrix)
    spearman = scipy.stats.spearmanr(matrix[:, 0], matrix[:, 1])[0]
    return 1 - spearman


def pearson_dist(u, v, na_vals=["NA", np.nan]):
    """
    Compute Pearson distance for vectors u, v.

    Returns 1 - spearmanr(u, v).
    """
#    matrix = [[x, y] for x, y in zip(u, v) \
#              if (x not in na_vals) and (y not in na_vals)]
#    matrix = array(matrix)
    pearson = scipy.stats.pearsonr(u, v)[0]
    return 1 - pearson


def pval_to_stars(pval, thresh=0.05):
    """
    Convert p-values to star notation (idiotic). 
    """
    stars = None
    if pval > thresh:
        stars = "NS"
    elif pval <= 10**(-4):
        stars = "****"
    elif pval <= 10**(-3):
        stars = "***"
    elif pval <= 10**(-2):
        stars = "**"
    elif pval <= thresh:
        stars = "*"
    assert (stars is not None), "Failed to get star notation."
    return stars
        


def my_pdist(X, dist_func,
             na_vals=["NA", np.nan, np.inf, -np.inf]):
    """
    (inefficient) pdist function that takes an arbitrary
    distance function (a lambda).

    X: data matrix
    dist_func: lambda that returns distance on vector
    na_vals: values to consider as missing data
    """
    X = array(X, dtype=object)
    if len(X.shape) == 1:
        num_rows = X.shape[0]
        num_cols = 1
    else:
        num_rows, num_cols = X.shape
    dist_matrix = []
    for col1 in range(num_cols):
        pdist_row = []
        for col2 in range(num_cols):
            pairs = array([[x, y] for x, y in zip(X[:, col1], X[:, col2]) \
                           if (x not in na_vals) and (y not in na_vals)])
            if len(pairs) == 0:
                continue
            dist = dist_func(pairs[:, 0],
                             pairs[:, 1])
            pdist_row.append(dist)
        dist_matrix.append(pdist_row)
    dist_matrix = array(dist_matrix)
    return dist_matrix


def leven_dist(first, second, na_vals=[]):
    """
    Levenshtein distance between two strings.
    """
    if (len(first) == 1) or len(second) == 1:
        raise Exception, "One of two vectors passed to leven dist has only 1 " \
                         "element in it."
    if len(first) > len(second):
        first, second = second, first
    if len(second) == 0:
        return len(first)
    first_length = len(first) + 1
    second_length = len(second) + 1
    distance_matrix = [[0] * second_length for x in range(first_length)]
    for i in range(first_length):
       distance_matrix[i][0] = i
    for j in range(second_length):
       distance_matrix[0][j]=j
    for i in xrange(1, first_length):
        for j in range(1, second_length):
            deletion = distance_matrix[i-1][j] + 1
            insertion = distance_matrix[i][j-1] + 1
            substitution = distance_matrix[i-1][j-1]
            if first[i-1] != second[j-1]:
                substitution += 1
            distance_matrix[i][j] = min(insertion, deletion, substitution)
    return distance_matrix[first_length-1][second_length-1]


def zscore_df(df, scale="row"):
    if scale == "row":
        # Get mean across columns (mean per row)
        df_mu = df.mean(axis=1)
        # Standard deviation across columns (sdev per row)
        df_sdev = df.std(axis=1)
        # Subtract mean across rows and divide result by sdev
        norm_df = df.sub(df_mu, axis=0).div(df_sdev, axis=0)
        return norm_df
    elif scale.startswith("col"):
        if len(df.index) == 0:
            raise Exception, "Cannot normalize by column with only one row."
        # Get mean across rows
        df_mu = df.mean(axis=0)
        # Standard deviation across columns (sdev per row)
        df_sdev = df.std(axis=0)
        # Subtract mean across rows and divide result by sdev
        norm_df = df.sub(df_mu, axis=1).div(df_sdev, axis=1)
        return norm_df
    else:
        raise Exception, "Unknown scale method %s" %(scale)
        

