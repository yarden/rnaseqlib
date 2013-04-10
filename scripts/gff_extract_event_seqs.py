import os
import sys
import re
import time

import numpy as np

from collections import defaultdict

import rnaseqlib
import rnaseqlib.gff.gffutils_helpers as gffutils_helpers
import rnaseqlib.utils as utils

import misopy
import misopy.gff_utils as gff_utils
import misopy.Gene as gene_utils

import pybedtools


def greeting():
    print "gff_extract_event_seqs:\nExtract sequences from GFF file " \
          "encoding alternative events.\n"
    print "Example:\n"
    print "gff_extract_event_seqs.py --input-gff events.gff3 --fi mm9.fasta " \
          "--output-dir event_seqs/\n"
    print "See --help for options.\n"


def main():
    from optparse import OptionParser
    parser = OptionParser()
      
    parser.add_option("--input-gff", dest="input_gff", default=None, nargs=1,
                      help="Fetch sequence from GFF events file. Takes as input: "
                      "GFF filename.")
    parser.add_option("--fi", dest="fasta_fname", default=None, nargs=1,
                      help="FASTA filename to fetch sequences from.")
    parser.add_option("--with-flanking-introns", dest="with_flanking_introns",
                      default=False, action="store_true",
                      help="Get sequence of flanking introns relative to skipped exon.")
    parser.add_option("--flanking-introns-coords", dest="flanking_introns_coords",
                      default=None, nargs=4,
                      help="Fetch the sequences of the flanking introns "
                      "(for SpliceGraph events). Takes as input the intervals to " 
                      "be used, which are: "
                      "(1) start position relative to 5 prime splice site of SE "
                      "(negative int), "
                      "(2) end position 5 prime splice site (negative int), "
                      "(3) start position relative to 3 prime splice site "
                      "(positive int), "
                      "(4) end position relative to 3 prime splice site. "
                      "(posiitve int). "
                      "Suggested settings are -250, -20, 20, -250.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.output_dir is None:
        greeting()
        print "Error: need --output-dir to be provided."
        sys.exit(1)
        
    output_dir = options.output_dir
    output_dir = utils.pathify(options.output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.input_gff is not None:
        if options.fasta_fname is None:
            greeting()
            print "Error: Must provide input fasta file with --fi."
            sys.exit(1)
        # Check for FASTA
        gff_filename = utils.pathify(options.input_gff)
        fasta_fname = utils.pathify(options.fasta_fname)
        flanking_introns_coords = options.flanking_introns_coords
        gffutils_helpers.fetch_seq_from_gff(gff_filename, fasta_fname, output_dir,
                                            with_flanking_introns=options.with_flanking_introns,
                                            flanking_introns_coords=options.flanking_introns_coords)
        
if __name__ == '__main__':
    main()
