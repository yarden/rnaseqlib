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


def write_rec_to_gff(gff_out, record):
    # Get our own attribute string to avoid duplicates
    gff_out.write("%s\n" %(str(record)))


def sanitize_record(record):
    if record.start > record.stop:
        record.start, record.stop = \
            record.stop, record.start
    return record


def fix_up_down_ids(parts):
    first_part, last_part = parts[0], parts[-1]
    # If it's a minus strand event, then the upstream exon
    # must have a greater coordinate than the downstream
    if (first_part.id.endswith(".dn") and \
        last_part.id.endswith(".up")) and \
       (first_part.stop <= last_part.start):
        first_part.attributes["ID"] = "%s.up" %(first_part.id[:-3])
        last_part.attributes["ID"] = "%s.dn" %(last_part.id[:-3])
    parts[0] = first_part
    parts[-1] = last_part
    return parts


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
    print "  - Writing sanitized GFF to: %s" %(gff_out_fname)
    t1 = time.time()
    for gene in gff_db.get_gene_objects():
        gene_rec = gene["gene_rec"]
        # Write gene record
        write_rec_to_gff(gff_out, sanitize_record(gene_rec))
        # Write mRNA of each gene
        for mRNA in gene["mRNAs"]:
            mRNA_obj = gene["mRNAs"][mRNA]
            write_rec_to_gff(gff_out, sanitize_record(mRNA_obj["record"]))
            parts = list([sanitize_record(part) \
                          for part in mRNA_obj["parts"]])
            fixed_parts = fix_up_down_ids(parts)
            for part in fixed_parts:
                write_rec_to_gff(gff_out, part)
    gff_out.close()
    t2 = time.time()
    print "Sanitizing took %.2f seconds" %(t2 - t1)
    return gff_out_fname, output_db_fname
