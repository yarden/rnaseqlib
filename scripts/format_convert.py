##
## Format conversions for various formats.
##
import os
import sys
import time

import fileinput

import rnaseqlib
import rnaseqlib.utils as utils

def bed_to_gff_line(line,
                    source,
                    field_type):
    """
    BED line to GFF line.
    """
    fields = line.strip().split("\t")
    num_fields = len(fields)
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = "."
    attributes = "."
    if num_fields >= 6:
        field_type = fields[3]
        # Do nothing with score for now
        score = fields[4]
        strand = fields[5]
    gff_fields = [chrom,
                  # Add 1 to start
                  start + 1,
                  end]
    gff_line = "\t".join(map(str, gff_fields))
    return gff_line


def bed_to_gff(args, source, field_type):
    """
    Convert input filename from BED to GFF.
    """
    if field_type is None:
        field_type = "exon"
    for line in fileinput.input(args):
        line = line.strip()
        if source is None:
            if fileinput.isstdin():
                source = "stdin"
            else:
                source = os.path.basename(fileinput.filename())
        output_line = bed_to_gff_line(line, source, field_type)
        sys.stdout.write("%s\n" %(output_line))


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--bed2gff", dest="bed_to_gff", action="store_true", default=False,
                      help="Convert BED to GFF format. BED is 0-based start, "
                      "GFF is 1-based start. Takes a filename or \'-\' to read input "
                      "from stdin. Outputs to stdout.")
    parser.add_option("--source", dest="source", default=None,
                      help="Source for resulting GFF.")
    parser.add_option("--type", dest="type", default=None,
                      help="Type for resulting GFF.")
    (options, args) = parser.parse_args()

    if options.bed_to_gff:
        bed_to_gff(args,
                   options.source,
                   options.type)

    

if __name__ == '__main__':
    main()
