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

    
def get_events_to_transcripts(gene_bed_filename,
                              # column corresponding to GFF attributes
                              gff_attribs_col=8,
                              # column corresponding to transcript name
                              transcript_col=12):
    """
    Return dictionary mapping event IDs to transcript names.
    """
    events_to_transcripts = defaultdict(list)
    bed_in = open(gene_bed_filename, "r")
    for line in bed_in:
        fields = line.strip().split("\t")
        gff_attribs = fields[gff_attribs_col]
        transcript = fields[transcript_col]
        attribs = gff_attribs.split(";")
        event_id = None
        for attr in attribs:
            if attr.startswith("ID="):
                event_id = attr.split("ID=")[1]
        if event_id == None:
            raise Exception, "ID= less event %s" %(line)
        # Mapping of GFF attributes to transcripts
        events_to_transcripts[event_id].append(transcript)
    bed_in.close()
    return events_to_transcripts

        
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


# def intersect_events_with_gff(events_filename,
#                               gene_table_filename,
#                               gene_bed_filename,
#                               output_dir,
#                               settings_info,
#                               table_sources=["ensembl_genes"],
#                               gene_keys=["geneSymbol",
#                                          "name2",
#                                          "desc"],
#                               # Key to resolve redundancies in entries with (e.g. ENSMG..)
#                               name_key="name2"):
#     """
#     Intersect GFF events with gene tables.
#     Return helpful fields from table.

#     The gene table file contains information about genes
#     (like symbol, desc, etc.)

#     The gene BED file (obtained from Tables in UCSC Genome Browser)
#     is just a BED format describing the transcripts of genes, which
#     is used for fast intersection with GFF file of events,
#     performed by bedtools.

#     """
#     print "Intersecting events with genes"
#     print "  - Events GFF: %s" %(events_filename)
#     print "  - Genes BED: %s" %(gene_bed_filename)
#     if not os.path.isdir(output_dir):
#         os.makedirs(output_dir)
#     # First make a BED file corresponding to intersection of
#     # GFF events with genes BED (transcripts)
#     intersected_bed = intersect_events_with_genes_bed(events_filename,
#                                                       gene_bed_filename,
#                                                       output_dir)
#     # Get events to transcripts
#     events_to_transcripts = get_events_to_transcripts(intersected_bed)

#     # Load genes table and supplement keys
#     gene_table = gt.GeneTable(settings_info)
#     num_events_with_multiples = 0
#     for source in table_sources:
#         table_keys = gene_table.table_fields[source]
#         event_genes_info = []
#         # Get genes from the current table source
#         genes = gene_table.genes[source]
#         # For each transcript, find its associated keys
#         for event_name, transcripts in events_to_transcripts.iteritems():
#             assert (len(transcripts) != 0), \
#                 "No transcripts for %s" %(event_name)
#             seen_name_values = {}
#             # Make a dictionary mapping name keys (gene IDs) to
#             # their values
#             gene_values = defaultdict(lambda: defaultdict(list))
#             for transcript in transcripts:
#                 transcript_info = genes.ix[transcript]
#                 # For each transcript, accumulate a list of all the
#                 # desired fields
#                 name_value = transcript_info[table_keys[name_key]]
#                 for needed_key in gene_keys:
#                     value = transcript_info[table_keys[needed_key]]
#                     if (type(value) == str) and value.endswith(","):
#                         value = value[:-1]
#                     gene_values[name_value][needed_key].append(value)
#             # Select best field value for each gene
#             gene_values = select_gene_values(gene_values)
#             all_gene_values = []
#             if len(gene_values) > 1:
#                 num_events_with_multiples += 1
#             for gene, gene_info in gene_values.iteritems():
#                 all_gene_values.append([gene_info[key] for key in gene_keys])
#             all_gene_values = zip(*all_gene_values)
#             # Converting fields to strings for entry in DataFrame later on
#             df_entry = [",".join(map(str, all_gene_values[idx])) \
#                         for idx in range(len(all_gene_values))]
#             df_entry = [event_name] + df_entry
#             event_genes_info.append(df_entry)
#         df_columns = ["event_name"] + gene_keys
#         event_genes_df = p.DataFrame(event_genes_info,
#                                      columns=df_columns)
#         ##
#         ## Output mapping from events to genes
#         ##
#         events_basename = os.path.basename(events_filename)
#         output_filename = "%s.%s.to_genes.txt" %(events_basename,
#                                                  source)
#         output_filename = os.path.join(output_dir,
#                                        output_filename)
#         print "Outputting to: %s" %(output_filename)
#         event_genes_df.to_csv(output_filename,
#                               sep="\t",
#                               index=False)
#         print "Total of %d events with multiple gene mappings." \
#             %(num_events_with_multiples)
        

        
def select_gene_values(gene_values):
    """
    Given a dictionary mapping genes to a dictionary of their
    keys and a list of values, select a single value that is
    best.  Usually pick the longest value (e.g. longest description of gene.)
    """
    selected_gene_values = defaultdict(dict)
    for gene, gene_info in gene_values.iteritems():
        for key, values in gene_info.iteritems():
            # Get index of maximum string-length element
            max_ind = utils.maxi(map(len, map(str, values)))
            # Use it
            selected_gene_values[gene][key] = values[max_ind]
    return selected_gene_values


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
                      "particular region.  Takes as input a GFF filename "
                      "followed by a chromosome region, e.g.: "
                      "SE.mm9.gff   chr:start:end")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
#    parser.add_option("--settings", dest="settings", nargs=1, default=None,
#                      help="Settings filename.")
    (options, args) = parser.parse_args()

    # Options that require output dir
    options_require_output_dir = [options.intersect]
    # Options that require settings filename
    #options_require_settings = [options.intersect]

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
    # Same for settings option
    #settings_info = None
    #parsed_settings = None
    # for given_opt in options_require_settings:
    #     if given_opt != None:
    #         if options.settings is None:
    #             print "Error: need --settings"
    #             sys.exit(1)
    #         else:
    #             settings_filename = os.path.abspath(os.path.expanduser(options.settings))
    #             print "Loading settings from: %s" %(settings_filename)
    #             # Parse settings
    #             settings_info, parsed_settings = \
    #                 settings.load_settings(settings_filename)

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
