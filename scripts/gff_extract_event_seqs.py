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


def get_flanking_introns_coords(gff_db,
                                mRNA_record):
    """
    Get coordinates of flanking intron regions around a
    skipped exon of interest.

    Regions are marked as follows:

    A: upstream region of upstream intron
    B: downstream region of upstream intron
    C: skipped exon itself
    D: upstream region of downstream intron
    E: downstream region of downstream intron

        A   B   C   D   E

    [  ]-----[ S E ]-----[  ]
    """
    mRNA_id = mRNA_record.get_id()
    exons = gff_db.exons_by_mRNA[mRNA_id]
    if len(exons) != 3:
        return None
    # Only extract information from exon trios, assuming
    # middle exon is skipped
    skipped_exon = exons[1]
    if mRNA_record.strand == "-":
        up_exon, dn_exon = exons[2], exons[0]
    else:
        up_exon, dn_exon = exons[0], exons[2]
    # Get all flanking intronic sequence
    up_intron_coords = (up_exon.end + 1,
                        skipped_exon.start - 1)
    dn_intron_coords = (skipped_exon.end + 1,
                        dn_exon.start - 1)
    regions = {"up": up_intron_coords,
               "dn": dn_intron_coords}
#    print "GOT: ", regions
    return regions
    # up_regions = flanking_introns["up"]
    # dn_regions = flanking_introns["dn"]
    # # Return regions in upstream and downstream
    # # intronic regions
    # regions = defaultdict(list)
    # # Retrieve upstream regions
    # regions["up"] = coords_from_relative_regions(skipped_exon,
    #                                              up_exon,
    #                                              up_regions)
    # # Retrieve downstream regions
    # regions["dn"] = coords_from_relative_regions(skipped_exon,
    #                                              dn_exon,
    #                                              dn_regions)

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
    return gene_parts
    
    
def fetch_seq_from_gff(gff_fname, fasta_fname, output_dir,
                       with_flanking_introns=False,
                       entries_to_include=["gene",
                                           "mRNA",
                                           "exon"]):
    """
    Fetch sequence from GFF file.

    Outputs:

    (1) GFF file containing an annotation of the sequences.

    (2) FASTA file with the actual sequences.

    If asked, fetch the flanking intronic sequences.
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
    gff_out = gff_utils.Writer(open(gff_output_fname, "w"))
    for gene_id in genes:
        gene_info = genes[gene_id]
        gene_tree = gene_info["hierarchy"]
        gene_obj = gene_info["gene_object"]
        gene_rec = gene_tree[gene_id]["gene"]
        ##
        ## GFF records to write for the current gene
        ##
        recs_to_write = []
        # For mRNA entries, extract the flanking introns of the
        # alternative exon if asked
        event_recs = get_event_recs_from_gene(gene_obj, gene_tree)
        if event_recs is None:
            continue
        if with_flanking_introns:
            flanking_intron_coords = \
                get_flanking_introns_coords(gff_db, long_mRNA_rec)
            if flanking_intron_coords == None:
                continue
            # Fetch upstream intron sequence
            up_intron_start, up_intron_end = flanking_intron_coords["up"]
            up_intron_seq = fetch_seq(chrom_file,
                                      up_intron_start,
                                      up_intron_end,
                                      chrom,
                                      strand)
            if up_intron_seq != None:
                seq_id = ":".join(map(str, [chrom, start, end, strand]))
                seq_id += "_up_intron"
                output_line = ">%s;%s\n%s\n" %(name, seq_id, up_intron_seq)
                output_file.write(output_line)
            else:
                print "WARNING: Could not fetch upstream intron seq %s, %s, %s" \
                    %(start, end, strand)
            # Fetch downstream intron sequence
            dn_intron_start, dn_intron_end = flanking_intron_coords["dn"]
            dn_intron_seq = fetch_seq(chrom_file,
                                      dn_intron_start,
                                      dn_intron_end,
                                      chrom,
                                      strand)
            if dn_intron_seq != None:
                seq_id = ":".join(map(str, [chrom, start, end, strand]))
                seq_id += "_dn_intron"
                output_line = ">%s;%s\n%s\n" %(name, seq_id, dn_intron_seq)
                output_file.write(output_line)
            else:
                print "WARNING: Could not fetch downstream intron seq %s, %s, %s" \
                    %(start, end, strand)

        # Retrieve the flanking intron coordinates, if asked
        #     seq_id = ":".join(map(str, [chrom, start, end, strand]))
        #     output_line = ">%s;%s;ID=%s;Parent=%s\n%s\n" \
        #         %(name, seq_id,
        #           rec_id,
        #           parent_id,
        #           seq)
        #     output_file.write(output_line)

    
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
                      "(negative int), "
                      "(4) end position relative to 3 prime splice site. "
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
                           with_flanking_introns=options.with_flanking_introns)
        
if __name__ == '__main__':
    main()
