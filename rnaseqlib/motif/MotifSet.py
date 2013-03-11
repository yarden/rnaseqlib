##
## MotifSet: represent a comparison between two samples
## for motifs or a sample with its dinucleotide shuffled self.
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

class MotifSet:
    def __init__(self, exp_coords_fname, control_fname=None):
        # Coordinates representing the experimental coordinates
        # (either a BED or a GFF file)
        self.exp_coords_fname = exp_coords_fname
        # Control coordinates filename (optional)
        self.control_fname = control_fname

    def foo():
        pass
