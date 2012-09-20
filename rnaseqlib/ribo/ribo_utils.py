##
## Utilities for Ribo-Seq data
##
import os
import sys
import time

import yklib
import yklib.fastq_utils

def rstrip_stretch(s, letter):
    """
    Strip (from right) consecutive stretch of letters.
    """
    stripped_s = ""
    first = True
    for l in s[::-1]:
        if first:
            if l == letter:
                continue
            else:
                stripped_s += l
                first = False
        else:
            stripped_s += l
    return stripped_s[::-1]
            

def trim_polyA(fastq_filename,
               output_dir,
               compressed=False,
               min_polyA_len=3):
    """
    Trim polyA stretches from reads
    """
    for line in fastq_utils.read_fastq(open(fastq_filename)):
        header, seq, header2, qual = line
        if not seq.endswith("A"):
            # Skip sequences that do not end with at least N
            # many As
            if seq[-min_polyA_len] != ("A" * min_polyA_len):
                continue
            stop_index = 0
            stripped_seq = rstrip_stretch(s[min_polyA_len:], "A")
            

def read_len_dist():
    """
    Compute distribution of read lengths.
    """
    pass
    
    
