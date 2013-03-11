##
## Utilities for working with jellyfish
##
import os
import sys
import time
import glob

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils

from collections import defaultdict

jf_path = utils.which("jellyfish")
    
def jf_merge(fname_base):
    """
    Merge all the jellyfish output files in the directory
    that they're in.

    Returns merged filename.
    """
    output_dir = os.path.dirname(fname_base)
    fname_pat = "%s_*" %(fname_base)
    jf_files = glob.glob(fname_pat)
    merged_fname = None
    if len(jf_files) == 0:
        raise Exception, "Cannot merge, no jf files with " \
                         "basename %s" %(fname_base)
    if len(jf_files) == 1:
        # No need to merge; only one jf file, so use
        # it
        merged_fname = jf_files[0]
    else:
        # There are multiple files to be merged
        merged_fname = "%s.merged.jf" %(fname_base)
        merge_cmd = "%s merge -o %s " %(fname_base)
    return merged_fname
    

def jf_count_kmers(fastx_fname, kmer_len,
                   output_dir,
                   hash_size=100000000):
    """
    Count kmers using jellyfish.
    """
    if not os.path.isfile(fastx_fname):
        print "Error: fastx file %s not found." %(fastx_fname)
        sys.exit(1)
    # Count kmers, use temporary file for db
    fastx_basename = os.path.basename(fastx_fname)
    db_fname = "%s.jf" %(os.path.join(output_dir, fastx_basename))
    output_fname = "%s_counts" %(db_fname)
    if os.path.isfile(db_fname):
        print "Overwriting %s" %(db_fname)
        os.remove(db_fname)
    if os.path.isfile(output_fname):
        print "Overwriting %s" %(output_fname)
        os.remove(output_fname)
    count_cmd = "%s count -m %d -o %s -s %d %s" \
        %(jf_path,
          kmer_len,
          db_fname,
          hash_size,
          fastx_fname)
    print "Counting kmers with jf: %s" %(count_cmd)
    ret_val = os.system(count_cmd)
    if ret_val != 0:
        raise Exception, "jellyfish count call failed."
        sys.exit(1)
    # Merge db results
    merged_fname = jf_merge(db_fname)
    # Load up kmer results
    dump_cmd = "%s dump -o %s %s" \
        %(jf_path,
          output_fname,
          merged_fname)
    ret_val = os.system(dump_cmd)
    return output_fname


def jf_counts_to_dict(jf_counts_fname):
    """
    Load jellyfish counts file (a FASTA file)
    into a dictionary mapping kmers to counts.
    """
    kmer_counts = defaultdict(int)
    if not os.path.isfile(jf_counts_fname):
        print "Error: Cannot find jf counts file %s" %(jf_counts_fname)
        sys.exit(1)
    for fastx_entry in fastx_utils.get_fastx_entries(jf_counts_fname,
                                                     fasta=True):
        kmer_count, kmer = fastx_entry
        # Remove prefix '>' from FASTA entry
        kmer_count = int(kmer_count[1:])
        kmer_counts[kmer] = kmer_count
    return kmer_counts


