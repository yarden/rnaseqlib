##
## Utilities for Ribo-Seq data
##
import os
import sys
import time

import rnaseqlib.fastq_utils

import yklib
import yklib.utils as utils

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
            

def trim_polyA_ends(fastq_filename,
                    output_dir,
                    compressed=False,
                    min_polyA_len=3):
    """
    Trim polyA ends from reads.
    """
    print "Trimming polyA trails from: %s" %(fastq_filename)
    # Strip the trailing extension
    output_basename = ".".join(os.path.basename(fastq_filename).split(".")[0:-1])
    output_basename = "%s.trimmed_polyA.fastq" %(output_basename)
    output_filename = os.path.join(output_dir, output_basename)
    utils.make_dir(output_dir)
    if os.path.isfile(output_filename):
        print "Error: %s already exists!" %(output_filename)
        sys.exit(1)
    print "  - Outputting trimmed sequences to: %s" %(output_filename)
    output_file = open(output_filename, "w")
    t1 = time.time()
    for line in fastq_utils.read_fastq(open(fastq_filename)):
        header, seq, header2, qual = line
        if not seq.endswith("A"):
            # Skip sequences that do not end with at least N
            # many As
            if seq[-min_polyA_len] != ("A" * min_polyA_len):
                continue
            # Get sequence stripped of contiguous strech of polyAs
            stripped_seq = rstrip_stretch(s[min_polyA_len:], "A")
            new_rec = (header, stripped_seq, header2, qual)
            # Write the record with trimmed sequence back out to file
            fastq_utils.write_fastq(output_file, new_rec)
    t2 = time.time()
    print "Trimming took %.2f mins." %((t2 - t1)/60.)
    output_file.close()
    return output_filename
            

def compute_read_len_dist():
    """
    Compute distribution of read lengths.
    """
    pass

if __name__ == "__main__":
    test_file = "/home/yarden/jaen/test_ribo.fastq"
    trim_polyA_ends(test_file, "/home/yarden/jaen/test_polyA/")
    
