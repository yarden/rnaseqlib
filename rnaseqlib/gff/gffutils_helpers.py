##
## Helper functions for Ryan Dale's gffutils
## library
##
import os
import sys
import time
import glob

import tempfile

import gffutils

import rnaseqlib
import rnaseqlib.utils as utils

import rnaseqlib.gff.GFFDB as gffdb


def get_default_db_fname(gff_fname, db_dirname="gff_db"):
    """
    Look for canonical GFF database filename. If exists,
    return its path, otherwise return None.

    Looks for that has 'gff_fname's basename ending in .db
    inside a 'gff_db' subdirectory in the same
    directory where 'gff_fname' is stored.

    For example, if 'gff_fname' is /home/user/mygff.gff it will
    look for /home/user/gff_db/mygff.gff.db.
    """
    gff_fname = utils.pathify(gff_fname)
    # If the input ends in .db, assume it is the database
    if gff_fname.endswith(".db"):
        return gff_fname
    gff_basename = os.path.basename(gff_fname)
    gff_db_dir = os.path.join(os.path.dirname(gff_fname),
                              db_dirname)
    if not os.path.isdir(gff_db_dir):
        return None
    db_fname = os.path.join(gff_db_dir, "%s.db" %(gff_basename))
    if not os.path.isfile(db_fname):
        return None
    return db_fname


def get_output_db_fname(gff_fname, output_dir,
                        db_subdir="gff_db"):
    """
    Return output file for 'gff_fname' (either a GFF db or
    a regular GFF file) in 'output_dir'
    """
    gff_fname = utils.pathify(gff_fname)
    gff_basename = os.path.basename(gff_fname)
    db_fname = \
        os.path.join(output_dir, db_subdir, gff_basename)
    if not db_fname.endswith(".db"):
        db_fname += ".db"
    return db_fname


def get_temp_db_fname():
    """
    Return a named temporary db filename.
    """
    pass


def create_db(gff_fname, db_fname=None):
    """
    Create a GFF database with gffutils. Returns the name
    of the gff_fname. If no output database name db_fname is
    given, use a named temporary file.
    """
    print "Creating a GFF database..."
    print "  - Input GFF: %s" %(gff_fname)
    print "  - Output file: %s" %(db_fname)
    if os.path.isfile(db_fname):
        print "GFF database %s exists. Quitting." %(db_fname)
        sys.exit(0)
    gffutils.create_db(gff_fname, db_fname)


def sanitize_gff(db_fname, output_db_fname, output_dir,
                 event_sanitize=False):
    """
    Sanitize a GFF file. Return the revised
    GFF file and its database.
    """
    gff_out_fname = \
        os.path.join(output_dir,
                     os.path.basename(db_fname).split(".db")[0])
    if not utils.endsin_gff(gff_out_fname):
        gff_out_fname += ".gff"
    gff_out = open(gff_out_fname, "w")
    gff_db = gffdb.GFFDB(db_fname)
    gff_db.load_all_genes()
    for rec in gff_db.iter_recs():
        print rec
#    for gene_rec in gff_db.iter_by_type("gene"):
#        gene_id = gene_rec.id
#        print mRNA
#    for record in gff_db.db.all():
#        if record.start > record.stop:
#            # Switch coords of start > stop
#            record.start, record.end = record.end, record.start
    # Create database for new GFF filename
    # ...
    return gff_out_fname, output_db_fname
