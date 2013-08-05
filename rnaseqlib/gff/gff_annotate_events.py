##
## Annotate GFF events with gene information based on a UCSC table.
##
## This script takes a GFF file with alternative event annotations
## (produced by gff_make_annotations) and a combined kgXref UCSC table
## (produced by rnaseqlib --init) and supplements the GFF event annotation
## with gene information. Specifically, it adds Ensembl gene IDs, RefSeq
## IDs, and gene symbols to the attributes field of the events GFF.
##
## This is done by converting by the UCSC table txStart/txEnd sites to
## a GFF internally and intersecting the event gene GFF coordinates with
## these txStart/txEnd coordinates (using pybedtools). The results of
## the intersection are converted.
##
import os
import sys
import time

import pandas

import rnaseqlib
import rnaseqlib.utils as utils

from collections import defaultdict

import gffutils
import gffutils.gffwriter as gffwriter
import pybedtools


def get_table_as_bedtool(table_fname):
    """
    Load UCSC table as BedTool where txStart/txEnd
    as coordinates.

    The UCSC tables (when downloaded) have a 0-based
    start coordinate, so need to add 1 to start.

    Uses the kgXref.combined ensGene table produced by
    
    """
    if "kgXref.combined" not in os.path.basename(table_fname):
        print "WARNING: Are you sure %s is a combined ensGene table?" \
              %(table_fname)
    table = pandas.read_table(table_fname, sep="\t")
    # For gene symbol: Use "value" if available, otherwise use geneSymbol
    if "value" in table.columns:
        gene_symbol_col = "value"
    else:
        gene_symbol_col = "geneSymbol"
    def get_bedtool_iter():
        for gene_num, gene_entry in table.iterrows():
            chrom = gene_entry["chrom"]
            start = int(gene_entry["txStart"]) + 1
            end = int(gene_entry["txEnd"])
            strand = gene_entry["strand"]
            # Annotation fields
            name2 = gene_entry["name2"]
            if pandas.isnull(name2):
                name2 = "NA"
            refseq_id = gene_entry["refseq"]
            if pandas.isnull(refseq_id):
                refseq_id = "NA"
            gene_symbol = gene_entry[gene_symbol_col]
            if pandas.isnull(gene_symbol):
                gene_symbol = "NA"
            attributes = \
                "ID=%s;ensg_id=%s;refseq_id=%s;gsymbol=%s;" \
                %(name2,
                  name2,
                  refseq_id,
                  gene_symbol)
            # Convert table to BedTool
            gff_entry = [chrom,
                         "genes_table",
                         "gene",
                         str(start),
                         str(end),
                         ".",
                         strand,
                         ".",
                         attributes]
            gff_interval = \
                pybedtools.create_interval_from_list(gff_entry)
            yield gff_interval
    table_bedtool = pybedtools.BedTool(get_bedtool_iter())
    return table_bedtool


def get_coords_to_gene_info(coord_intervals, table_bed):
    """
    Given an iterator of BedTool Intervals 'coord_intervals'
    and a BedTool containing table coordinates 'table_bed'
    return a mapping from events to genes.

    Map event genes to their IDs. If event_geneX
    is an entry of coords_intervals, generate:
    
       event_gene1 -> refseq  -> value
                   -> ensgene -> value
       event_gene2 -> refseq  -> ...
       ...
    """
    # Mapping from event intervals (i.e. the genes in 'coord_intervals')
    # and the gene information
    event_genes_to_info = defaultdict(lambda: defaultdict(list))
    # Intersect event genes with gene txStart/txEnd
    intersected_bed = \
        coord_intervals.intersect(table_bed, wb=True, s=True, f=1)
    for entry in intersected_bed:
        event_gene_attrs = utils.parse_attributes(entry.fields[8])
        event_gene_str = event_gene_attrs["ID"]
        gene_info_field = entry.fields[-1]
        # Strip semicolon of ID attributes
        if gene_info_field.endswith(";"):
            gene_info_field = gene_info_field[0:-1]
        # Convert attributes into dictionary
        gene_info = utils.parse_attributes(gene_info_field)
        ensgene_id = gene_info["ensg_id"]
        refseq_id = gene_info["refseq_id"]
        gene_symbol = gene_info["gsymbol"]
        # Skip null entries
        if not is_null_id(ensgene_id):
            event_genes_to_info[event_gene_str]["ensg_id"].append(ensgene_id)
        if not is_null_id(refseq_id):
            event_genes_to_info[event_gene_str]["refseq_id"].append(refseq_id)
        if not is_null_id(gene_symbol):
            event_genes_to_info[event_gene_str]["gsymbol"].append(gene_symbol)
    return event_genes_to_info


def annotate_gff_with_genes(args):
    """
    Annotate GFF with genes table.
    """
    gff_fname = utils.pathify(args.gff_filename)
    if not os.path.isfile(gff_fname):
        raise Exception, "Cannot find %s" %(gff_fname)
    table_fname = utils.pathify(args.table_filename)
    if not os.path.isfile(table_fname):
        raise Exception, "Cannot find %s" %(table_fname)
    table_bed = get_table_as_bedtool(table_fname)
    # Get BedTool for events, containing only the gene entries
    all_events_bed = pybedtools.BedTool(gff_fname)
    event_genes = \
        all_events_bed.filter(lambda entry: entry.fields[2] == "gene")
    print "Determining overlap between events and genes..."
    # Get mapping from events to gene information
    event_genes_to_info = \
      get_coords_to_gene_info(event_genes, table_bed)

    # Incorporate the gene information into the GFF and output it
    # it using gffutils
    print "Loading events into GFF database..."
    events_db = gffutils.create_db(gff_fname, ":memory:",
                                   verbose=False)
    output_fname = gff_fname 
    events_out = gffwriter.GFFWriter(output_fname,
                                     in_place=True)
    print " - Outputting annotated GFF to: %s" %(output_fname)
    def new_recs():
        for gene_recs in list(events_db.iter_by_parent_childs()):
            gene_rec = gene_recs[0]
            event_id = gene_rec.id
            # Use existing IDs if present
            if "ensgene_id" in gene_rec.attributes:
                ensgene_id = gene_rec.attributes["ensg_id"][0]
            else:
                ensgene_id = "NA"
            if "refseq_id" in gene_rec.attributes:
                refseq_id = gene_rec.attributes["refseq_id"][0]
            else:
                refseq_id = "NA"
            if "gene_symbol" in gene_rec.attributes:
                gene_symbol = gene_rec.attributes["gsymbol"][0]
            else:
                gene_symbol = "NA"
            if event_id in event_genes_to_info:
                event_info = event_genes_to_info[event_id]
                ensgene_ids = \
                    utils.unique_list(event_info["ensg_id"])
                if len(ensgene_ids) > 0 and ensgene_ids[0] != "NA":
                    ensgene_id = ",".join(ensgene_ids)
                refseq_ids = \
                    utils.unique_list(event_info["refseq_id"])
                if len(refseq_ids) > 0 and refseq_ids[0] != "NA":
                    refseq_id = ",".join(refseq_ids)
                gene_symbols = \
                    utils.unique_list(event_info["gsymbol"])
                if len(gene_symbols) > 0 and gene_symbols[0] != "NA":
                    gene_symbol = ",".join(gene_symbols)
            gene_rec.attributes["ensg_id"] = [ensgene_id]
            gene_rec.attributes["refseq_id"] = [refseq_id]
            gene_rec.attributes["gsymbol"] = [gene_symbol]
            # Yield all the gene's records
            for g in gene_recs:
                yield g
    t1 = time.time()
    print "Creating annotated GFF database..."
    annotated_db = gffutils.create_db(new_recs(), ":memory:",
                                      verbose=False)
    t2 = time.time()
    print "Creation took %.2f secs" %(t2 - t1)
    # Write to file
    print "Writing annotated GFF to file..."
    for gene_rec in annotated_db.all_features(featuretype="gene"):
        events_out.write_gene_recs(annotated_db, gene_rec.id)
    events_out.close()


def is_null_id(gene_id):
    """
    Return True if the ID is considered null.
    """
    if (gene_id == "NA") or (gene_id == "nan") or pandas.isnull(gene_id):
        return True
    return False
    

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_filename",
                        help="GFF filename to annotate with gene information.")
    parser.add_argument("table_filename",
                        help="Table contains txStart/txEnd sites and the "
                        "gene fields.")
    parser.add_argument("--in-place", default=False, action="store_true",
                        help="If passed, outputs annotation in place (i.e. "
                        "overwriting the passed in file.) Also sanitizes "
                        "the GFF.")
    args = parser.parse_args()
    annotate_gff_with_genes(args)


if __name__ == "__main__":
    main()
