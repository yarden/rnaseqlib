##
## Sanitize a GFF file using Ryan Dale's gffutils
##
## Ensure that start < end for all coordinates. Annotate
## each entry with 'gene_id'
## 
##
import os
import sys
import time

import gffutils

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.gff.gffutils_helpers as gffutils_helpers


def greeting():
    print "gff_sanitize:\n\tSanitize a GFF database file."
    print "See --help for options."


def sanitize_gff(db_fname, output_dir, db_subdir):
    """
    Sanitize a GFF file.
    """
    db_fname = gffutils_helpers.get_default_db_fname(db_fname)
    print "Sanitizing GFF file..."
    if db_fname is None:
        print "Error: Cannot find a GFF database for %s. Please index " \
              "and try again." %(db_fname)
        sys.exit(1)
    print "  - Input GFF: %s" %(db_fname)
    print "  - Output directory: %s" %(output_dir)
    # Iterate through GFF file
    output_db_fname = \
        gffutils_helpers.get_output_db_fname(db_fname,
                                             output_dir,
                                             db_subdir=db_subdir)
    gffutils_helpers.sanitize_gff(db_fname, output_db_fname, output_dir)


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--input-gff", dest="input_gff", default=None, nargs=1,
                      help="Input GFF filename for a GFF database.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    parser.add_option("--gtf", dest="gtf", default=False,
                      action="store_true",
                      help="Output a GTF file instead of GFF.")
    parser.add_option("--db-subdir", dest="db_subdir", default="gff_db",
                      help="Name of output subdirectory containing GFF " \
                      "database. By default, creates \'gff_db\' " \
                      "subdirectory in the directory given to --output-dir.")
    parser.add_option("--no-db-output", dest="no_db_output", default=False,
                      action="store_true",
                      help="Do not output a GFF database.")

    (options, args) = parser.parse_args()

    if options.output_dir is None:
        print "Error: need --output-dir to be provided.\n"
        greeting()
        sys.exit(1)

    output_dir = options.output_dir
    output_dir = utils.pathify(options.output_dir)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.input_gff is not None:
        gff_fname = utils.pathify(options.input_gff)
        sanitize_gff(gff_fname, output_dir,
                     db_subdir=options.db_subdir)


if __name__ == "__main__":
    main()

