##
## Extract lengths of all GFF files
##
import os
import sys
import time

import numpy as np
import pandas

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.Gene as gene_utils


def extract_lens_from_gff(gff_fname, output_dir):
    entries = []
    output_basename = "%s.lens" %(os.path.basename(gff_fname))
    output_fname = os.path.join(output_dir, output_basename)
    print "Extracting lengths from GFF file..."
    print "  - Input GFF: %s" %(gff_fname)
    print "  - Output file: %s" %(output_fname)
    if os.path.isfile(output_fname):
        print "Overwriting %s" %(output_fname)
    gff_genes = gene_utils.load_genes_from_gff(gff_fname)
    for gene_id in gff_genes:
        gene = gff_genes[gene_id]["gene_object"]
        # Get the length of each isoform
        iso_lens = []
        iso_labels = [] 
        iso_exon_lens = []
        genomic_coords = []
        for iso in gene.isoforms:
            iso_lens.append(str(iso.len))
            iso_labels.append(iso.label)
            exon_lens = [str(exon.len) for exon in iso.parts]
            iso_exon_lens.append(exon_lens)
            genomic_coords.append([iso.genomic_start,
                                   iso.genomic_end])
        genomic_coords = np.array(genomic_coords)
        genomic_lens = \
            map(str, list(genomic_coords[:, 1] - genomic_coords[:, 0] + 1))
        entry = \
            {"event_name":
             gene.label,
             "mRNA_lens":
             ",".join(iso_lens),
             "mRNA_labels":
             ",".join(iso_labels),
             "exon_lens":
             ";".join([",".join(exons) for exons in iso_exon_lens]),
             "genomic_lens":
             ",".join(genomic_lens)}
        entries.append(entry)
    entries_df = pandas.DataFrame(entries)
    entries_df.to_csv(output_fname,
                      cols=["event_name", "mRNA_labels",
                            "mRNA_lens", "exon_lens",
                            "genomic_lens"],
                      sep="\t",
                      index=False)


#def make_gff_db(gff_fname, output_dir):
#    """
#    Index the GFF if its database doesn't already exist.
#    """
#    if not os.path.isdir(output_dir):
#        os.makedirs(output_dir)
#    if not os.path.isfile(gff_fname):
#        print "GFF %s not found."
#        sys.exit(1)

def greeting():
    print "gff_extract_lens:\n\tExtract lengths from GFF file"
    print "See --help for options."


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--input-gff", dest="input_gff", default=None, nargs=1,
                      help="Extract lengths from GFF file.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.output_dir is None:
        print "Error: need --output-dir to be provided.\n"
        greeting()
        sys.exit(1)

    output_dir = options.output_dir
    output_dir = utils.pathify(options.output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.input_gff is not None:
        gff_fname = utils.pathify(options.input_gff)
        extract_lens_from_gff(gff_fname, output_dir)


if __name__ == "__main__":
    main()

