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
## Utilities for converting bam to UCSC formats like
## bigWig
## 
def bam_to_bigWig_file(bam_fname, bigWig_fname, genome):
    """
    Convert a bam file to bigWig
    """
    try:
        import pybedtools
        from pybedtools.contrib.bigwig import bam_to_bigwig
    except:
        print "Cannot convert BAM to Wig without pybedtools."
        return None
    bam_to_bigwig(bam=bam_fname,
                  genome=genome,
                  output=bigWig_fname)
    return None
    

##
## Utilities for extracting FASTX sequences
## from BAM files. Not intended for handling
## of paired-end reads.
##
def bam_to_rec(in_file,
               make_unique_recs=False):
    """
    Generator to convert BAM files into Biopython SeqRecords.
    """
    from Bio import SeqIO, Seq, SeqRecord
    bam_file = pysam.Samfile(in_file, "rb")
    rec_num = 1
    # Keep track of which read IDs have been outputted
    read_ids_outputted = {}
    for read in bam_file:
        seq = Seq.Seq(read.seq)
        if read.is_reverse:
            seq = seq.reverse_complement()
        read_name = read.qname
        if make_unique_recs:
            read_name = "%s_%d" %(read_name, rec_num)
        else:
            # If we're not asked to make the records unique,
            # then don't output the same read ID twice
            if read_name in read_ids_outputted:
                continue
            # Record that we've seen this read
            read_ids_outputted[read_name] = True
        rec = SeqRecord.SeqRecord(seq, read_name, "", "")
        rec_num += 1
        yield rec


def bam_to_fastx(logger, in_file, out_file,
                 record_type="fasta",
                 make_unique_recs=False):
    """
    BAM to FASTX converter, based on code by Brad Chapman.

    By default converts to FASTA record.

    If 'make_unique_recs' is set to True, then make each FASTA
    record unique (append a number to it) so that reads with
    multiple alignments can be considered.
    """
    logger.info("Converting %s to FASTA" %(in_file))
    from Bio import SeqIO, Seq, SeqRecord
    out_handle = fastx_utils.write_open_fastx(out_file)
    logger.info("  - Output file: %s" %(out_file))
    SeqIO.write(bam_to_rec(in_file, make_unique_recs=make_unique_recs),
                out_handle,
                record_type)
    logger.info("Finished FASTA conversion.")
    out_handle.close()



