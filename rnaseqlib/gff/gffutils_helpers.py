##
## Helper functions for Ryan Dale's gffutils
## library
##
import os
import sys
import time
import glob
import string

import tempfile
#import gffutils

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils
import rnaseqlib.gff.GFFDB as gffdb

import misopy
import misopy.gff_utils as miso_gff_utils
import misopy.Gene as gene_utils


def output_gff_event_seqs(event_ids, input_fasta_fname, output_fasta_fname,
                          entry_types=None,
                          suffixes=None,
                          remove_repeats=False):
    """
    Given a set of event ids, pull out their sequences from an
    input fasta filename and output these to a separate FASTA file.

    Return the entries that were outputted.

      - entry_types: optional list of entry types that should be outputted, 
        e.g. 'exon', 'intron'. Skip all entry types not within list.
      - suffixes: optional list of suffixes that the first field of the
        FASTA name should end in. For example, if the FASTA field is:

          >event_id;part_id;entry_type

        Then event_id must end in one of the suffixes for it to be
        included.
    """
    num_events = len(event_ids)
    print "Retrieving sequences for %d events" %(num_events)
    def is_event_fasta(fasta_name):
        """
        Return true if the event is a FASTA one.
        """
        # If there's any event such that the
        # FASTA record starts with that event's name, then
        # the FASTA record should be outputted
        return len(filter(lambda e: \
                          fasta_name.startswith(e),
                          event_ids)) > 0
    kept_fasta_entries = []
    print "Fetching: ", entry_types, suffixes
    with open(output_fasta_fname, "w") as fasta_out:
        for entry in fastx_utils.get_fastx_entries(input_fasta_fname):
            fasta_name, fasta_seq = entry
            fasta_name_fields = fasta_name.split(";")
            entry_type = fasta_name_fields[2]
            if is_event_fasta(fasta_name[1:]):
                # If given entry types, check that this sequence
                # is of one of the right entry types; otherwise
                # skip it
                if (entry_types is not None) and \
                   (entry_type not in entry_types):
                    # Not of correct entry type
                    continue
                # If given suffixes, check that the first field
                # of the FASTA name ends in one of the suffixes
                if suffixes is not None:
                    if not any([fasta_name_fields[0].endswith(s) \
                                for s in suffixes]):
                        # The first FASTA name field does not end
                        # in any of the suffixes, so skip it.
                        continue
                # If asked, remove repeats from sequence
                if remove_repeats:
                    repeatless_seq = \
                        fasta_seq.translate(None, string.ascii_lowercase)
                    if len(repeatless_seq) == 0:
                        print "%s is all repeat! Not removing" %(fasta_name)
                        continue
                    fasta_seq = repeatless_seq
                fasta_out.write("%s\n" %(fasta_name))
                fasta_out.write("%s\n" %(fasta_seq))
                kept_fasta_entries.append(fasta_name)
    print "Outputted %d entries." %(len(kept_fasta_entries))
    return kept_fasta_entries
    

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


#def gffutils_write_rec_to_gff(gff_out, record):
#    # Get our own attribute string to avoid duplicates
#    gff_out.write("%s\n" %(str(record)))


def write_rec_to_gff(gff_out, record):
    # Get our own attribute string to avoid duplicates
    gff_out.write(record)


#def gffutils_sanitize_record(record):
#    if record.start > record.stop:
#        record.start, record.stop = \
#            record.stop, record.start
#    return record


def sanitize_record(record):
    if record.start > record.end:
        record.start, record.end = \
            record.end, record.start
    return record

# def gffutils_fix_up_down_ids(parts):
#     first_part, last_part = parts[0], parts[-1]
#     # If it's a minus strand event, then the upstream exon
#     # must have a greater coordinate than the downstream
#     if (first_part.id.endswith(".dn") and \
#         last_part.id.endswith(".up")) and \
#        (first_part.stop <= last_part.start):
#         first_part.attributes["ID"] = "%s.up" %(first_part.id[:-3])
#         last_part.attributes["ID"] = "%s.dn" %(last_part.id[:-3])
#     parts[0] = first_part
#     parts[-1] = last_part
#     return parts


def fix_up_down_ids(parts):
    first_part, last_part = parts[0], parts[-1]
    # If it's a minus strand event, then the upstream exon
    # must have a greater coordinate than the downstream
    if (first_part.get_id().endswith(".dn") and \
        last_part.get_id().endswith(".up")) and \
       (first_part.end <= last_part.start):
        first_part.attributes["ID"] = ["%s.up" %(first_part.get_id()[:-3])]
        last_part.attributes["ID"] = ["%s.dn" %(last_part.get_id()[:-3])]
    parts[0] = first_part
    parts[-1] = last_part
    return parts


# def gffutils_sanitize_gff(db_fname, output_db_fname, output_dir,
#                  event_sanitize=False):
#     """
#     Sanitize a GFF file. Return the revised
#     GFF file and its database.
#     """
#     gff_out_fname = \
#         os.path.join(output_dir,
#                      os.path.basename(db_fname).split(".db")[0])
#     if not utils.endsin_gff(gff_out_fname):
#         gff_out_fname += ".gff"
#     gff_out = open(gff_out_fname, "w")
#     gff_db = gffdb.GFFDB(db_fname)
#     print "  - Writing sanitized GFF to: %s" %(gff_out_fname)
#     t1 = time.time()
#     for gene in gff_db.get_gene_objects():
#         gene_rec = gene["gene_rec"]
#         if gene_rec.id != "chr6:71510891:71511070:+@chr6:71524109:71524230:+@chr6:71529635:71531580:+":
#             continue
#         print "GENE: ", gene
#         print "mRNAs: ", gene["mRNAs"]
#         print "values: ", gene["mRNAs"].values
#         # Write gene record
#         write_rec_to_gff(gff_out, sanitize_record(gene_rec))
#         # Write mRNA of each gene
#         for mRNA in gene["mRNAs"]:
#             print "writing out parts of %s" %(mRNA)
#             mRNA_obj = gene["mRNAs"][mRNA]
#             write_rec_to_gff(gff_out, sanitize_record(mRNA_obj["record"]))
#             parts = list([sanitize_record(part) \
#                           for part in mRNA_obj["parts"]])
#             print "  - parts: ", parts
#             print "  - parts parents: ", [p.attributes["Parent"] for p in parts]
# #            fixed_parts = fix_up_down_ids(parts)
#             fixed_parts = parts
#             for part in fixed_parts:
#                 write_rec_to_gff(gff_out, part)
#     gff_out.close()
#     t2 = time.time()
#     print "Sanitizing took %.2f seconds" %(t2 - t1)
#     return gff_out_fname, output_db_fname


def add_introns_to_gff(gff_filename, output_dir):
    """
    Add 'intron' entries to GFF.
    """
    output_basename = \
        utils.trim_gff_ext(os.path.basename(gff_filename))
    ext_to_use = os.path.basename(gff_filename).rsplit(".", 1)[1]
    output_filename = \
        os.path.join(output_dir,
                     "%s.with_introns.%s" %(output_basename,
                                            ext_to_use))
    print "Adding introns to GFF..."
    print "  - Input: %s" %(gff_filename)
    print "  - Output: %s" %(output_filename)
    if os.path.isfile(output_filename):
        print "Found file %s, skipping.." %(output_filename)
        return output_filename
    gff_out = miso_gff_utils.Writer(open(output_filename, "w"))
    gff_db = miso_gff_utils.GFFDatabase(from_filename=gff_filename,
                                        reverse_recs=True)
    t1 = time.time()
    genes = gene_utils.load_genes_from_gff(gff_filename)
    for gene_id in genes:
        gene_info = genes[gene_id]
        gene_tree = gene_info["hierarchy"]
        gene_obj = gene_info["gene_object"]
        gene_rec = gene_tree[gene_id]["gene"]
        # Write the GFF record
        gff_out.write(gene_rec)
        # Write out the mRNAs, their exons, and then
        # input the introns
        for mRNA in gene_obj.isoforms:
            mRNA_id = mRNA.label
            curr_mRNA = gene_tree[gene_id]["mRNAs"][mRNA_id]
            gff_out.write(curr_mRNA["record"])
            # Write out the exons
            curr_exons = gene_tree[gene_id]["mRNAs"][mRNA_id]["exons"]
            for exon in curr_exons:
                gff_out.write(curr_exons[exon]["record"])
        # Now output the introns
        for isoform in gene_obj.isoforms:
            intron_coords = []
            for first_exon, second_exon in zip(isoform.parts,
                                               isoform.parts[1::1]):
                # Intron start coordinate is the coordinate right after
                # the end of the first exon, intron end coordinate is the
                # coordinate just before the beginning of the second exon
                intron_start = first_exon.end + 1
                intron_end = second_exon.start - 1
                if intron_start >= intron_end:
                    continue
                intron_coords.append((intron_start, intron_end))
                # Create record for this intron
                intron_id = "%s:%d-%d:%s.intron" \
                    %(gene_obj.chrom,
                      intron_start,
                      intron_end,
                      gene_obj.strand)
                intron_rec = \
                    miso_gff_utils.GFF(gene_obj.chrom, gene_rec.source, "intron",
                                       intron_start, intron_end,
                                       strand=gene_obj.strand,
                                       attributes={"ID": [intron_id],
                                                   "Parent": [isoform.label]})
                gff_out.write(intron_rec)
    t2 = time.time()
    print "Addition took %.2f minutes." %((t2 - t1)/60.)


def sanitize_gff(gff_fname, output_dir, include_introns=True):
    """
    Sanitize a GFF file. Return the revised
    GFF file.
    """
    gff_out_fname = os.path.join(output_dir, os.path.basename(gff_fname))
    genes = gene_utils.load_genes_from_gff(gff_fname,
                                           include_introns=include_introns)
    t1 = time.time()
    with open(gff_out_fname, "w") as gff_out_file:
        gff_out = miso_gff_utils.Writer(gff_out_file)
        for gene in genes:
            gene_obj = genes[gene]["gene_object"]
            gene_record = genes[gene]["hierarchy"][gene]["gene"]
            gene_hierarchy = genes[gene]["hierarchy"][gene]
            # Write gene record
            write_rec_to_gff(gff_out, sanitize_record(gene_record))
            for mRNA in gene_obj.isoforms:
                mRNA_id = mRNA.label
                mRNA_obj = gene_hierarchy["mRNAs"][mRNA_id]
                mRNA_record = mRNA_obj["record"]
                # Write out the mRNA record
                write_rec_to_gff(gff_out, sanitize_record(mRNA_record))
                # Get parts of each mRNA
                parts_records = \
                    [mRNA_obj["exons"][p.label]["record"] for p in mRNA.parts]
                if gene_obj.strand == "-":
                    parts_records = fix_up_down_ids(parts_records)
                for part in parts_records:
                    write_rec_to_gff(gff_out, part)
    t2 = time.time()
    print "Sanitizing took %.2f seconds" %(t2 - t1)
    return gff_out_fname


def get_fasta_from_gff(input_gff, output_dir):
    """
    Get FASTA sequences from GFF file.
    """
    pass
