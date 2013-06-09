##
## Cogent-based sequence counter
##
import os
import sys
import time

import cogent
from cogent.core.usage import DinucUsage

import rnaseqlib
import rnaseqlib.utils as utils


class SeqCounter:
    """
    Sequence counter.
    """
    def __init__(self, seqs):
        self.seqs = seqs



