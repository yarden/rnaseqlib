##
## BAM-related utilities
##
import os
import sys
import time

import pysam

import rnaseqlib
import rnaseqlib.fastx_utils as fastx_utils


##
## Utilities for extracting FASTX sequences
## from BAM files. Not intended for handling
## of paired-end reads.
##
def bam_to_rec(in_file):
    """
    Generator to convert BAM files into Biopython SeqRecords.
    """
    from Bio import SeqIO, Seq, SeqRecord
    bam_file = pysam.Samfile(in_file, "rb")
    for read in bam_file:
        seq = Seq.Seq(read.seq)
        if read.is_reverse:
            seq = seq.reverse_complement()
        rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
        yield rec


def bam_to_fastx(in_file, out_file, record_type="fasta"):
    """
    BAM to FASTX converter, based on code by Brad Chapman.

    By default converts to FASTA record.
    """
    from Bio import SeqIO, Seq, SeqRecord
    out_file = "%s.fa" %(os.path.splitext(in_file)[0])
    out_handle = fastx_utils.write_open_fastx(out_file)
    SeqIO.write(bam_to_rec(in_file), out_handle, record_type)
    out_handle.close()



