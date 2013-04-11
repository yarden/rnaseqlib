##
## Helper functions for Ryan Dale's gffutils
## library
##
import os
import sys
import time
import glob
import re

import tempfile
#import gffutils
import string
from string import maketrans

import pybedtools

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils

import misopy
import misopy.gff_utils as miso_gff_utils
import misopy.Gene as gene_utils


def coords_from_relative_regions(skipped_exon,
                                 adjacent_exon,
                                 coords):
    """
    Translate relative coordinates to actual mRNA coordinates.

    Positive coordinates are relative to 5' splice site
    Negative coordinates are relative to 3' splice site
    """
    translated_coords = []
    for coord_pair in coords:
        if coords[0] >= 0:
            # Positive coordinates are read as distance
            # from end of adjacent exon
            offset = adjacent_exon.end
        elif coords < 0:
            # Negative coordinates are read as distance
            # start of skipped exon
            offset = skipped_exon.start
        offset_coords = (offset + coord_pair[0],
                         offset + coord_pair[1])
        if offset_coords[0] > offset_coords[1]:
            # Reverse coordinates to make first smaller
            offset_coords = offset_coords[::-1]
        translated_coords.append(offset_coords)
    return translated_coords


def get_flanking_introns_coords(gene_obj):
    """
    Get coordinates of flanking intron regions around a
    skipped exon of interest.
    """
    long_mRNA, short_mRNA = \
        gene_obj.isoforms[0], gene_obj.isoforms[1]
    if long_mRNA.len < short_mRNA:
        raise Exception, "%s must have long mRNA first." \
              %(gene_obj.label)
    exons = long_mRNA.parts
    if len(exons) != 3:
        return None
    # Only extract information from exon trios, assuming
    # middle exon is skipped
    up_exon = exons[0]
    skipped_exon = exons[1]
    dn_exon = exons[2]
    
    # Get all flanking intronic sequence
    up_intron_coords = (up_exon.end + 1,
                        skipped_exon.start - 1)
    dn_intron_coords = (skipped_exon.end + 1,
                        dn_exon.start - 1)
    regions = {"up_intron": up_intron_coords,
               "dn_intron": dn_intron_coords}
    return regions

def get_event_recs_from_gene(gene_obj, gene_tree):
    """
    Given a gene object a tree of gene records,
    return the GFF records corresponding to the
    following event parts:

      upstream exon, skipped exon, downstream exon

    Return None if gene cannot be parsed
    """
    gene_parts = {}
    gene_id = gene_obj.label
    num_isoforms = len(gene_obj.isoforms)
    # Only consider two-isoform events from GFF
    if num_isoforms != 2:
        print "Skipping %s since it does not have two isoforms." \
              %(gene_id)
        return gene_parts
    # Take the long mRNA (usually the first mRNA)
    long_mRNA, short_mRNA = gene_obj.isoforms[0], gene_obj.isoforms[1]
    if long_mRNA.len < short_mRNA.len:
        print "WARNING: Needed to swap mRNAs %s and %s" \
              %(long_mRNA.label, short_mRNA.label)
        raise Exception
        # Swap mRNAs if first is shorter than second
        long_mRNA, short_mRNA = short_mRNA, long_mRNA
    long_mRNA_id = long_mRNA.label
    long_mRNA_rec = gene_tree[gene_id]["mRNAs"][long_mRNA_id]["record"]
    # Get the GFF record for the current mRNA
    long_mRNA_tree = gene_tree[gene_id]["mRNAs"][long_mRNA_id]
    # Check that the mRNA has three exons
    if len(long_mRNA.parts) != 3:
        print "Error: Need mRNA %s to have exactly three exons." \
              %(long_mRNA_id)
        print "  - Skipping.."
        return gene_parts
    # Find the upstream, skipped and downstream exons
    # Up exon
    up_exon_id = long_mRNA.parts[0].label
    up_exon_rec = \
        long_mRNA_tree["exons"][up_exon_id]
    # Skipped exon
    se_exon_id = long_mRNA.parts[1].label
    se_exon_rec = \
        long_mRNA_tree["exons"][se_exon_id]
    # Down exon
    dn_exon_id = long_mRNA.parts[2].label
    dn_exon_rec = \
        long_mRNA_tree["exons"][dn_exon_id]
    gene_parts["up_exon"] = up_exon_rec
    gene_parts["se_exon"] = se_exon_rec
    gene_parts["dn_exon"] = dn_exon_rec
    # Record long mRNA rec
    gene_parts["long_mRNA"] = long_mRNA_rec
    return gene_parts
    
    
def fetch_seq_from_gff(gff_fname, fasta_fname, output_dir,
                       with_flanking_introns=False,
                       flanking_introns_coords=None,
                       entries_to_include=["gene",
                                           "mRNA",
                                           "exon"]):
    """
    Fetch sequence from GFF file.

    Outputs:

    (1) GFF file containing an annotation of the sequences.

    (2) FASTA file with the actual sequences.

    If asked, fetch the flanking intronic sequences.

    Flanking regions are marked below:

      U: region of upstream intron
      D: region of downstream intron

             U           D

    [ U P ]-----[ S E ]-----[ D N ]

            a,b         c,d

    a,b,c,d correspond to optional flanking intron coordinates
    that determine the regions of the upstream/downstream
    introns that should be fetched:

       a, b: negative ints, position relative to 5' splice site of SE
             a < b

       c, d: positive ints, position relative to 3' splice site of SE
             c < d
    """
    # Load GFF genes
    gff_db = miso_gff_utils.GFFDatabase(from_filename=gff_fname,
                                        reverse_recs=True)
    file_basename = re.sub("\.gff3?", "",
                           os.path.basename(gff_fname))
    output_basename = "%s.event_seqs" %(file_basename)
    if flanking_introns_coords is not None:
        output_basename = "%s.flank_intronic_%s_%s_%s_%s" \
            %(output_basename,
              flanking_introns_coords[0],
              flanking_introns_coords[1],
              flanking_introns_coords[2],
              flanking_introns_coords[3])
    gff_outdir = os.path.join(output_dir, "gff_coords")
    utils.make_dir(gff_outdir)
    gff_output_fname = os.path.join(gff_outdir, "%s.gff" %(output_basename))
    fasta_output_fname = os.path.join(output_dir, "%s.fa" %(output_basename))
    print "Outputting GFF coordinates to: %s" %(gff_output_fname)
    if os.path.isfile(gff_output_fname):
        print "  - Overwriting existing file"
    print "Outputting sequences to: %s" %(fasta_output_fname)
    if os.path.isfile(fasta_output_fname):
        print "  - Overwriting existing file"
    genes = gene_utils.load_genes_from_gff(gff_fname)
    gff_out_file = open(gff_output_fname, "w")
    gff_out = miso_gff_utils.Writer(gff_out_file)
    for gene_id in genes:
        gene_info = genes[gene_id]
        gene_tree = gene_info["hierarchy"]
        gene_obj = gene_info["gene_object"]
        # GFF records to write for the current gene
        recs_to_write = []
        # For mRNA entries, extract the flanking introns of the
        # alternative exon if asked
        event_recs = get_event_recs_from_gene(gene_obj, gene_tree)
        long_mRNA_id = event_recs["long_mRNA"].get_id()
        if event_recs is None:
            continue
        # Write out up, se, and dn exons
        recs_to_write.extend([event_recs["up_exon"]["record"],
                              event_recs["se_exon"]["record"],
                              event_recs["dn_exon"]["record"]])
        if with_flanking_introns:
            introns_coords = \
                get_flanking_introns_coords(gene_obj)
            if introns_coords == None:
                raise Exception, "Cannot find flanking introns coordinates."
                sys.exit(1)
            # Fetch upstream intron sequence
            up_intron_start, up_intron_end = \
                introns_coords["up_intron"]
            up_intron_len = up_intron_end - up_intron_start + 1
            # Fetch downstream intron sequence
            dn_intron_start, dn_intron_end = \
                introns_coords["dn_intron"]
            dn_intron_len = dn_intron_end - dn_intron_start + 1
            # If given custom coordinates, use them instead of entire up/down
            # flanking intronic coordinates.
            se_exon_rec = event_recs["se_exon"]["record"]
            if flanking_introns_coords is not None:
                # (start,end) of upstream intron sequence
                a, b = \
                    int(flanking_introns_coords[0]), int(flanking_introns_coords[1])
                c, d = \
                    int(flanking_introns_coords[2]), int(flanking_introns_coords[3])
                a, b, c, d = error_check_intronic_coords(a, b, c, d,
                                                         up_intron_len, dn_intron_len)
                # Coordinates relative to 5' splice site of sequence to be fetched
                # The start of upstream intron sequence is negative from the 5' ss
                up_intron_start = se_exon_rec.start + a
                up_intron_end = se_exon_rec.start + b
                dn_intron_start = se_exon_rec.end + c
                dn_intron_end = se_exon_rec.end + d
            # Make GFF records for up/dn intronic sequences
            chrom = se_exon_rec.seqid
            source = se_exon_rec.source
            rec_type = "intron"
            strand = se_exon_rec.strand
            up_intron_str = "%s.up_intron" %(long_mRNA_id)
            up_intron_rec = \
                miso_gff_utils.GFF(chrom, source, "intron",
                              up_intron_start, up_intron_end,
                              strand=strand,
                              attributes={"ID": [up_intron_str],
                                          "Parent": [gene_obj.label]})
            dn_intron_str = "%s.dn_intron" %(long_mRNA_id)
            dn_intron_rec = \
                miso_gff_utils.GFF(chrom, source, "intron",
                                   dn_intron_start, dn_intron_end,
                                   strand=strand,
                                   attributes={"ID": [dn_intron_str],
                                               "Parent": [gene_obj.label]})
            recs_to_write.append(up_intron_rec)
            recs_to_write.append(dn_intron_rec)
        # Write out records to GFF
        for rec in recs_to_write:
            gff_out.write(rec)
    gff_out_file.close()
    # Output FASTA sequences
    output_fasta_seqs_from_gff(gff_output_fname,
                               fasta_fname,
                               fasta_output_fname)
    return fasta_output_fname


def output_fasta_seqs_from_gff(gff_fname,
                               fasta_input_fname,
                               fasta_output_fname,
                               s=True,
                               name=True,
                               use_gff_id=True):
    """
    Output FASTA sequence from GFF.
    """
    def assign_id_to_gff_field(gff_rec):
        rec_coords = "%s:%s-%s:%s" %(gff_rec[0],
                                     str(gff_rec.start),
                                     str(gff_rec.stop),
                                     gff_rec.strand)
        gff_type = gff_rec[2]
        gff_rec[2] = "%s;%s;%s" %(gff_rec.attrs["ID"],
                                  rec_coords,
                                  gff_type)
        return gff_rec
    gff_tool = pybedtools.BedTool(gff_fname)
    if use_gff_id:
        # Use the GFF ID= field to label the FASTA sequences
        gff_tool = gff_tool.each(assign_id_to_gff_field)
    seqs = gff_tool.sequence(fi=fasta_input_fname,
                             s=s,
                             fo=fasta_output_fname,
                             name=True)


def error_check_intronic_coords(a, b, c, d,
                                up_intron_len,
                                dn_intron_len,
                                trim_len=10):
    """
    Error check relative intronic coordinates used to
    determine subparts of upstream/downstream intronic
    sequences to be fetched.
    """
    assert (a < b < c < d), "a, b, c, d must be increasing in order."
    assert (a < b), "a must be less than b" %(a, b)
    assert (c < d), "c must be less than d" %(c, d)
    assert (a < 0), "a must be negative."
    assert (b < 0), "b must be negative."
    assert (c > 0), "c must be positive."
    assert (d > 0), "d must be positive."
    # Check that the sequence requested does not exceed intron
    # length
    up_diff = abs(abs(a) - up_intron_len) 
    if (abs(a) > up_intron_len):
        a -= (up_diff + trim_len)
    dn_diff = abs(d - dn_intron_len)
    if (d > dn_intron_len):
        d += (dn_diff - trim_len)
    return a, b, c, d

    
def fetch_seq(chrom_file,
              start, end,
              chrom, strand,
              revcompl=maketrans("actguACTGU|",
                                 "tgacaTGACA|")):
    """
    Fetch sequence from file handle.
    """
    seq = None
    try:
        # Convert to 0-based and fetch sequence
        chrom_file.seek(int(start) - 1)
        seq = chrom_file.read(int(end) - int(start) + 1).upper()
    except:
        print chrom, start, end, " not found"
        return None
    if strand == '-':
        # If on minus strand, take the reverse complement
        seq = seq.translate(revcompl)[::-1]
    return seq


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
    print "  - Input FASTA: %s" %(input_fasta_fname)
    print "  - Output FASTA: %s" %(output_fasta_fname)
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
