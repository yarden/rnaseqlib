##
## Create GFF databases using Ryan Dale's gffutils
##
import os
import sys
import time

import gffutils

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.gff.gffutils_helpers as gffutils_helpers


def greeting():
    print "gff_create_db:\n\tCreate SQLite database for GFF file.\n"
    print "Usage:"
    print "\tgff_create_db.py --input-gff mygff.gff --output-dir dirname\n"
    print "See --help for options."


def create_db(gff_fname, output_dir):
    """
    Create a GFF database for a given GFF filename.
    """
    output_basename = os.path.basename(gff_fname)
    db_fname = os.path.join(output_dir, "%s.db" %(output_basename))
    gffutils_helpers.create_db(gff_fname, db_fname)


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--input-gff", dest="input_gff", default=None, nargs=1,
                      help="Create a database for input GFF filename. Takes a " \
                      "GFF filename.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
#    parser.add_option("--gtf", dest="gtf", default=False, action="store_true",
#                      help="Output file as GTF. Default is GFF.")
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
        if not os.path.isfile(gff_fname):
            print "Error: GFF file %s does not exist." %(gff_fname)
            sys.exit(1)
        create_db(gff_fname, output_dir)
 

if __name__ == "__main__":
    main()

