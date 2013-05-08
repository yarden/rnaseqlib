##
## Make a GFF annotation of alternative splicing
## events from a set of UCSC tables.
##
import os
import sys
import time
import argparse

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.events.defineEvents as def_events

def load_ucsc_tables(tables_dir):
    """
    Load UCSC tables from a directory.
    """
    fnames = [os.path.join(tables_dir, f) \
              for f in os.listdir(tables_dir)]
    table_fnames = [fname for fname in fnames \
                    if os.path.isfile(fname)]
    return table_fnames
    

def make_annotation(args):
    """
    Make GFF annotation. Takes GFF tables directory
    and an output directory.
    """
    tables_dir = utils.pathify(args.tables_dir)
    output_dir = utils.pathify(args.output_dir)
    print "Making GFF alternative events annotation..."
    print "  - UCSC tables read from: %s" %(tables_dir)
    print "  - Output dir: %s" %(output_dir)
    t1 = time.time()
    table_fnames = load_ucsc_tables(tables_dir)
    print "Loaded %d UCSC tables." %(len(table_fnames))
    def_events.defineAllSplicing(tables_dir, output_dir,
                                 flanking=args.flanking_rule,
                                 multi_iso=args.multi_iso)
    t2 = time.time()
    print "Took %.2f minutes to make the annotation." \
          %((t2 - t1)/60.)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("tables_dir",
                        help="Directory where UCSC tables are. These are "
                        "used in making the annotation.")
    parser.add_argument("output_dir", help="Output directory.")
    parser.add_argument("--flanking-rule", default="commonshortest",
                        help="Rule to use when defining exon trios. "
                        "E.g. \'commonshortest\' to use the most common "
                        "and shortest regions are flanking exons to an "
                        "alternative trio.")
    parser.add_argument("--multi-iso", action="store_true",
                        default=False, help="If passed, generates "
                        "multi-isoform annotations. Off by default.")
    args = parser.parse_args()
    make_annotation(args)

