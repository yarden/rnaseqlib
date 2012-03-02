##
## Quality control metrics for RNA-Seq
##
import os
import sys
import time
import utils

from numpy import *
from scipy import *

def compute_qc_metrics(bam_filename, output_dir, settings):
    """
    Compute QC metrics.
    """
    pass


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-qc", dest="compute_qc", nargs=1, default=None,
                      help="Compute QC metrics for RNA-Seq library. Takes as input a BAM file.")
    parser.add_option("--settings", dest="settings", nargs=1, default=None,
                      help="Settings filename.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    output_dir = options.output_dir
    if output_dir != None:
        output_dir = pathify(output_dir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    

if __name__ == '__main__':
    main()
