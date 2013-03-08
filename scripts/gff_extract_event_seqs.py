import os
import sys
import re
import time

import numpy as np

from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.gff_utils as gff_utils
import misopy.Gene as gene_utils

import string
from string import maketrans

import pybedtools


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
    print "UP EXON: ", up_exon
    print "DN EXON: ", dn_exon
    
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
    gff_db = gff_utils.GFFDatabase(from_filename=gff_fname,
                                   reverse_recs=True)
    file_basename = re.sub("\.gff3?", "",
                           os.path.basename(gff_fname))
    output_basename = os.path.join(output_dir,
                                   file_basename)
    gff_output_fname = "%s.gff" %(output_basename)
    fasta_output_fname = "%s.fa" %(output_basename)
    print "Outputting GFF coordinates to: %s" %(gff_output_fname)
    if os.path.isfile(gff_output_fname):
        print "  - Overwriting existing file"
    print "Outputting sequences to: %s" %(fasta_output_fname)
    if os.path.isfile(fasta_output_fname):
        print "  - Overwriting existing file"
    genes = gene_utils.load_genes_from_gff(gff_fname)
    gff_out_file = open(gff_output_fname, "w")
    gff_out = gff_utils.Writer(gff_out_file)
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
                gff_utils.GFF(chrom, source, "intron",
                              up_intron_start, up_intron_end,
                              strand=strand,
                              attributes={"ID": [up_intron_str],
                                          "Parent": [gene_obj.label]})
            dn_intron_str = "%s.dn_intron" %(long_mRNA_id)
            dn_intron_rec = \
                gff_utils.GFF(chrom, source, "intron",
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
#    output_fasta_seqs_from_gff(gff_output_fname, fasta_output_fname)


def output_fasta_seqs_from_gff(gff_fname, fasta_output_fname):
    gff_tool = pybedtools.BedTool(gff_output_fname)
    gff_tool.sortBed()
    


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
        print "\'a\' coordinate must be less than %d (upstream intron length.)" \
              %(up_intron_len)
        print "Trimming by %d nt" %(trim_len)
        # Make 'a' shorter by adding to it
        a -= (up_diff + trim_len)
    dn_diff = abs(d - dn_intron_len)
    if (d > dn_intron_len):
        print "\'d\' coordinate must be less than %d (downstream intron length." \
              %(dn_intron_len)
        print "Trimming by %d nt" %(trim_len)
        # Make 'd' shorter by subtracting from it
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
        # if flanking_introns != None:
        #     flanking_introns = map(int, flanking_introns)            
        #     # Use same regions for upstream intron and downstream intron
        #     flanking_introns = {"up": ((flanking_introns[0], flanking_introns[1]),
        #                                (flanking_introns[2], flanking_introns[3])),
        #                         "dn": ((flanking_introns[0], flanking_introns[1]),
        #                                (flanking_introns[2], flanking_introns[3]))}
        fetch_seq_from_gff(gff_filename, fasta_fname, output_dir,
                           with_flanking_introns=options.with_flanking_introns,
                           flanking_introns_coords=options.flanking_introns_coords)
        
if __name__ == '__main__':
    main()
