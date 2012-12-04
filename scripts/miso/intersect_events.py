##
## Intersect GFF events file with a gene table
##
import os
import re
import sys
import time

import numpy as np
import pandas as p

from collections import defaultdict

import misopy
import misopy.gff_utils as gff_utils

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.tables as tables
import rnaseqlib.mapping.bedtools_utils as bedtools_utils


def gff_genes_from_events(gff_filename,
                          gene_field="gene"):
    """
    Return command that only selects genes from events
    GFF filename.
    """
    cmd = "grep -w %s %s" %(gene_field,
                            gff_filename)
    return cmd


def intersect_events_with_bed(events_filename,
                              bed_filename,
                              output_dir):
    """
    Intersect events with coordinates in BED format.
    Return the output filename of bedtools intersectBed.

    - events_filename: events filename (GFF)
    - bed_filename: coordinates to intersect with (BED)
    - output_dir: Output directory
    """
    print "Intersecting events with BED coordinates..."
    print "  - Events file: %s" %(events_filename)
    print "  - BED file: %s" %(bed_filename)
    output_dir = os.path.join(output_dir, "intersected_bed")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # Pass only genes from events GFF
    just_genes_cmd = gff_genes_from_events(events_filename)
    # Intersect the genes GFF with the BED file, taking only
    # stranded hits. 
    # use -wao: output also event features w/o overlaps
    intersect_bed_cmd = "intersectBed -a stdin -b %s -wao -s" \
      %(bed_filename)
    #    intersect_bed_cmd = "intersectBed -a stdin -b %s -wb -s" \
    #        %(bed_filename)
    basename = "%s_%s" %(os.path.basename(events_filename),
                         os.path.basename(bed_filename))
    basename = re.sub("[.]gff3?", "", basename)
    basename = basename.replace(".bed", "")
    output_filename = "%s.bed" %(os.path.join(output_dir, basename))
    print "  - Output file: %s" %(output_filename)
    if os.path.isfile(output_filename):
        print "Found %s. Skipping.." %(output_filename)
        return output_filename
    intersect_cmd = "%s | %s > %s" %(just_genes_cmd,
                                     intersect_bed_cmd,
                                     output_filename)
    print "Executing: %s" %(intersect_cmd)
    os.system(intersect_cmd)
    return output_filename


def get_events_to_genes(intersected_bed_filename,
                        gff_id_col=8,
                        gene_id_col=12,
                        na_val="NA"):
    """
    Parse the output file resulting from intersectBed
    that mapped events to genes and parse its results.

    Return a mapping from events to the genes they map to.

    - intersected_bed_filename: filename outputted by 
      intersectBed.

    - gff_id_col: the column number encoding the GFF ID= for 
      the event 
    - gene_id_col: the column number encoding the gene ID for
      the event
    - na_val: NA value, used when an event does not map to gene
    """
    events_to_genes = defaultdict(list)
    with open(intersected_bed_filename) as bed_in:
        for line in bed_in:
            fields = line.strip().split("\t")
            gff_event_id_field = fields[gff_id_col]
            gene_id = fields[gene_id_col]
            if not gff_event_id_field.startswith("ID="):
                raise Exception, "Malformed event line: %s" \
                  %(line)
            # If the event does not map to a gene, record its
            # gene id as a missing value
            if gene_id == ".":
                gene_id = na_val
            gff_event_id = gff_event_id_field.split("ID=")[1].split(";")[0]
            events_to_genes[gff_event_id].append(gene_id)
    return events_to_genes
    
        
def output_inclusive_trans_coords(gene_table, output_dir):
    """
    Output inclusive trans coords (in BED format)
    for the genes in the given gene table.

    Returns the output filename for the BED file.
    """
    print "Outputting inclusive transcript coordinates..."
    bed_basename = "%s_trans_coords.bed" %(gene_table.source)
    bed_output_fname = os.path.join(output_dir, bed_basename)
    print "  - Output file: %s" %(bed_output_fname)
    if os.path.isfile(bed_output_fname):
        print "Found %s. Skipping..." %(bed_output_fname)
        return bed_output_fname
    with open(bed_output_fname, "w") as bed_out:
        # Output the inclusive transcript coordinates for each gene
        for gene_id, gene in gene_table.genes.iteritems():
            chrom = gene.chrom
            start, end = gene.get_inclusive_trans_coords()
            name = gene_id
            score = "1"
            strand = gene.strand
            bed_line = bedtools_utils.make_bed_line(chrom, start, end,
                                                    name, score, strand)
            bed_out.write("%s\n" %(bed_line))
    return bed_output_fname
        
        

def intersect_events_with_genes(events_gff_fname,
                                gene_tables_dir,
                                output_dir,
                                genes_source="ensGene",
                                na_val="NA"):
    """
    Intersect GFF events with a genes table (also in GFF format).

    Computes the outermost transcription start/end bounds for each
    genes and then intersects the GFF events with these bounds.

    Outputs a mapping from event ID to one or more genes IDs that it
    maps to, if the event overlaps an annotated gene.

    - events_gff_fname: GFF events filename
    - gene_tables_dir: Directory with gene tables (created by --init module
      of rnaseqlib)
    - output_dir: output directory
    - genes_source: source of genes table, e.g. ensGene or refGene.
      By default, assumes input is an Ensembl table.
    """
    utils.make_dir(output_dir)
    events_basename = os.path.basename(events_gff_fname)
    events_to_genes_fname = \
      os.path.join(output_dir, "%s_to_%s.txt" \
                   %(events_basename, 
                     genes_source))
    print "Outputting events to genes..."
    print "  - Output file: %s" %(events_to_genes_fname)
    if os.path.isfile(events_to_genes_fname):
        print "Found %s. Skipping.." 
        return events_to_genes_fname
    # Load the gene table without parsing the individual genes
    gene_table = tables.GeneTable(gene_tables_dir, genes_source)
    # Create a BED file containing the most inclusive txStart/txEnd
    # for each gene in the table
    bed_coords_fname = output_inclusive_trans_coords(gene_table, 
                                                     output_dir)
    # Intersect the GFF events with this BED file of coordinates 
    # to determine what genes each event overlaps
    intersected_bed_fname = \
      intersect_events_with_bed(events_gff_fname,
                                bed_coords_fname,
                                output_dir)
    # Parse the resulting intersectBed results to get a mapping
    # from events to the genes they map to
    events_to_genes = get_events_to_genes(intersected_bed_fname)
    # Output the result to a file
    with open(events_to_genes_fname, "w") as events_to_genes_out:
        header = "event_id\tgene_id\n"
        events_to_genes_out.write(header)
        for event, genes in events_to_genes.iteritems():
            genes_str = ",".join(genes)
            output_line = "%s\t%s\n" %(event, genes_str)
            events_to_genes_out.write(output_line)


def parse_query_region(region):
    if ":" not in region:
        print "Error: malformed query region %s" %(region)
        sys.exit(1)
    parsed_region = region.split(":")
    if len(parsed_region) < 2:
        print "Error: malformed query region %s - need at least chrom, " \
            "start, and end" %(region)
        sys.exit(1)
    chrom = parsed_region[0]
    coords = parsed_region[1].split("-")
    start = int(coords[0])
    end = int(coords[1])
    strand = None
    if len(parsed_region) == 3:
        strand = parsed_region[2]
    if start > end:
        print "Error: start must be greater than end in %s" %(region)
        sys.exit(1)
    parsed_region = [chrom,
                     start,
                     end,
                     strand]
    return parsed_region
        
        
def get_events_in_region(gff_filename, region,
                         record_types=["gene"]):
    """
    Output Return all 'gene' entries in a given GFF file that
    intersect the given region.

    record_types is a list of GFF records to collect (e.g. gene, mRNA, ...)
    """
    gff_db = gff_utils.GFFDatabase(from_filename=gff_filename,
                                   reverse_recs=True)
    # Parse the query region
    parsed_region = parse_query_region(region)
    query_chrom, query_start, query_end, \
        query_strand = parsed_region
    matched_records = []
    num_recs = 0
    for record in gff_db:
        chrom = record.seqid
        # Name
        name = record.type
        start, end = int(record.start), int(record.end)
        strand = record.strand
        # Skip GFF records that don't match our record types
        if name not in record_types:
            continue
        num_recs += 1
        # Check that there is intersection
        if (query_chrom != chrom) or \
            (not utils.intersect_coords(query_start, query_end,
                                        start, end)):
            # Skip if chromosomes don't match or if there's no intersection
            continue
        # If strand is supplied in query region, check that
        # the strand matches 
        if (query_strand is not None) and \
            (strand != query_strand):
            continue
        # Must match
        record_id = record.get_id()
        print "%s" %(record_id)
        print "  - ", record
        matched_records.append(record)
    print "Looked through %d records." %(num_recs)
    return matched_records
        

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--intersect", dest="intersect", default=None, nargs=2,
                      help="Intersect events with GFF. Takes: (1) an events "
                      "filename (GFF), and (2) a gene tables directory "
                      "(produced by --init of rnaseqlib).")
    parser.add_option("--events-in-region", dest="events_in_region", default=None,
                      nargs=2,
                      help="Return all gene entries in a GFF that match a "
                      "particular region. Takes as input a GFF filename "
                      "followed by a chromosome region, e.g.: "
                      "SE.mm9.gff   chr:start:end")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    # Options that require output dir
    options_require_output_dir = [options.intersect]

    # Check that output dir is given if we're called with options
    # that need it
    output_dir = None
    for given_opt in options_require_output_dir:
        if given_opt != None:
            if options.output_dir == None:
                print "Error: need --output-dir"
                sys.exit(1)
            else:
                output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    if options.intersect != None:
        event_filename = os.path.abspath(os.path.expanduser(options.intersect[0]))
        gene_tables_dir = os.path.abspath(os.path.expanduser(options.intersect[1]))        
        intersect_events_with_genes(event_filename,
                                    gene_tables_dir,
                                    output_dir)
        
    if options.events_in_region != None:
        event_filename = os.path.abspath(os.path.expanduser(options.events_in_region[0]))
        region = options.events_in_region[1]
        get_events_in_region(event_filename, region)
    

if __name__ == "__main__":
    main()
