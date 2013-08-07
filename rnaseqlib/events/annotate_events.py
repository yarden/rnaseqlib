##
## Annotate GFF events with gffutils
##
import os
import sys
import time

import pandas

import rnaseqlib
import rnaseqlib.utils as utils

import gffutils
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
    def get_bedtool_iter():
        seen_genes = {}
        for gene_num, gene_entry in table.iterrows():
            chrom = gene_entry["chrom"]
            start = int(gene_entry["start"]) + 1
            end = int(gene_entry["end"])
            strand = gene_entry["strand"]
            # Annotation fields
            name2 = gene_entry["name2"]
            knownGene_name = gene_entry["knownGene_name"]
            gene_symbol = gene_entry["value"]
            # Skip genes we've seen already
            if name2 in seen_genes:
                continue
            # Record gene
            seen_genes[name2] = True
            attributes = \
                "ID=%s;ensgene_id=;knowngene_id=;gsymbol=;" \
                %(name2,
                  name2,
                  knownGene_name,
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
    print "TABLE BED: ", table_bed
    

def main():
    pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_filename",
                        help="GFF filename to annotate with gene information.")
    parser.add_argument("table_filename",
                        help="Table contains txStart/txEnd sites and the "
                        "gene fields.")
    parser.add_argument("--flanking-rule", default="commonshortest",
                        help="Rule to use when defining exon trios. "
                        "E.g. \'commonshortest\' to use the most common "
                        "and shortest regions are flanking exons to an "
                        "alternative trio.")
    parser.add_argument("--in-place", default=False, action="store_true",
                        help="If passed, outputs annotation in place (i.e. "
                        "overwriting the passed in file.) Also sanitizes "
                        "the GFF.")
    args = parser.parse_args()
    make_annotation(args)
