##
## Quality control metrics for RNA-Seq
##
import os
import sys
import time
import utils

from numpy import *
from scipy import *

import rnaseqlib
from rnaseqlib.QualityControl import QualityControl

def compute_qc_metrics(settings_filename, output_dir, settings):
    """
    Compute QC metrics.
    """
    pipeline = Pipeline(settings_filename, output_dir)
    # Compute the QC metrics for each sample
    qc_obj = QualityControl()


def get_cycle_profile(fastq_filename):
    qc_obj = QualityControl(None)
    percent_n = qc_obj.get_seq_cycle_profile(fastq_filename,
                                             first_n_seqs=1000000)
    print "Percent n for %s" %(fastq_filename)

    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compute-qc", dest="compute_qc", nargs=1, default=None,
                      help="Compute QC metrics for RNA-Seq library. Takes as input a BAM file.")
    parser.add_option("--get-cycle-profile", dest="get_cycle_profile", nargs=1, default=None,
                      help="Get cycle profile for sequence. Takes as input a FASTQ file.")
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

    if options.get_cycle_profile != None:
        fastq_filename = utils.pathify(options.get_cycle_profile)
        get_cycle_profile(fastq_filename)
        
    

if __name__ == '__main__':
    main()
