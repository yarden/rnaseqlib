##
## Utilities for clustering
##
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram

from numpy import *

import rnaseqlib
import rnaseqlib.stats.stats_utils as stats_utils

def hierarchical_clust(data_matrix,
                       dist_func,
                       linkage_method,
                       normalize=False):
    """
    Wrapper for hierarchical clustering.
    """
    print "Hierarchical clustering..."
    dist_matrix = stats_utils.my_pdist(data_matrix,
                                       dist_func)
    print "Computing linkage (method = %s).." %(linkage_method)
    linkage_matrix = linkage(squareform(dist_matrix),
                             linkage_method)
    clustering = {"linkage": linkage_matrix,
                  "dist": dist_matrix}
    return clustering

#...
#...
#...
#...
