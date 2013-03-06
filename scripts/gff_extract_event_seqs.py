import os
import sys
import re
import time

from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils import utils

import misopy
import misopy.gff_utils as gff_utils

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
    
def fetch_seq_from_gff(gff_f, genome_dir, output_dir,
                       with_flanking_introns=False,
                       entries_to_include=["gene",
                                           "mRNA",
                                           "exon"]):
    """
    Fetch sequence from GFF file.

    If asked, fetch the flanking intronic sequences.
    """
    # Mapping from chromosome IDs to their filenames
    files = {}
    # Load GFF file
    #gff_in = gff_utils.Reader(open(gff_f, "r"))
    gff_in = gff_utils.GFFDatabase(from_filename=gff_f,
                                   reverse_recs=True)
    file_basename = "%s.fa" %(re.sub("\.gff3?", "",
                                     os.path.basename(gff_f)))
    output_filename = os.path.join(output_dir,
                                   file_basename)
    print "Outputting sequences to: %s" %(output_filename)
    output_file = open(output_filename, "w")
    for record in gff_in:
        # Chromosome
        chrom = record.seqid
        # Type
        source = record.type
        # Name
        name = record.type
        start, end = record.start, record.end
        strand = record.strand
        attributes = record.attributes
        rec_id = attributes["ID"][0]
        if "Parent" not in attributes:
            parent_id = ""
        else:
            parent_id = attributes["Parent"][0]
        # If it's not the kind of entry we were looking for,
        # skip it
        if source not in entries_to_include:
            continue
        # Skip irrelevant chromosomes
        if ('random' in chrom) or ('j' in chrom) or ('_' in chrom):
            continue
        # Open chromosome file if not open already
        if chrom not in files:
            chrom_filename = os.path.join(genome_dir, "%s.fa" %(chrom))
            if not os.path.isfile(chrom_filename):
                print "WARNING: %s does not exist" %(chrom_filename)
                continue
            else:
                print "Opening: %s" %(chrom_filename)
            files[chrom] = open(chrom_filename, "r")
        chrom_file = files[chrom]
        seq = fetch_seq(chrom_file,
                        start, end,
                        chrom, strand)
        if seq != None:
            seq_id = ":".join(map(str, [chrom, start, end, strand]))
            output_line = ">%s;%s;ID=%s;Parent=%s\n%s\n" \
                %(name, seq_id,
                  rec_id,
                  parent_id,
                  seq)
            output_file.write(output_line)
            # For mRNA entries, extract the flanking introns of the
            # alternative exon if asked
            if with_flanking_introns:
                flanking_intron_coords = get_flanking_introns_coords(gff_in,
                                                                     record)
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

    output_file.close()

    
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



def main():
    from optparse import OptionParser
    parser = OptionParser()
      
    parser.add_option("--fetch", dest="fetch", default=None, nargs=2,
                      help="Fetch sequence from GFF file. Takes as input: "
                      "(1) GFF filename, (2) genome directory (with no newlines).")
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

    if options.output_dir == None:
        print "Error: need --output-dir to be provided."
        sys.exit(1)

    output_dir = options.output_dir
    output_dir = utils.pathify(options.output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.fetch != None:
        gff_filename = utils.pathify(options.fetch[0])
        genome_dir = utils.pathify(options.fetch[1])
        flanking_introns_coords = options.flanking_introns_coords
        # if flanking_introns != None:
        #     flanking_introns = map(int, flanking_introns)            
        #     # Use same regions for upstream intron and downstream intron
        #     flanking_introns = {"up": ((flanking_introns[0], flanking_introns[1]),
        #                                (flanking_introns[2], flanking_introns[3])),
        #                         "dn": ((flanking_introns[0], flanking_introns[1]),
        #                                (flanking_introns[2], flanking_introns[3]))}
        fetch_seq_from_gff(gff_filename, genome_dir, output_dir,
                           with_flanking_introns=options.with_flanking_introns)

        
if __name__ == '__main__':
    main()
