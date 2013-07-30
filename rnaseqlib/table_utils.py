##
## Utilities for manipulating init tables
##
import os
import sys
import time

from collections import defaultdict

import pandas

import pybedtools

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.tables as tables


def trans_to_gene_from_table(tables_dir):
    """
    Return mapping from transcript to gene given
    a pipeline init tables directory.
    """
    table_fname = os.path.join(tables_dir, "ensGene.kgXref.combined.txt")
    if not os.path.isfile(table_fname):
        raise Exception, "Cannot find combined table %s" %(table_fname)
    trans_to_gene = {}
    # Map transcripts to genes
    table_df = pandas.read_table(table_fname, sep="\t")
    for row, entry in table_df.iterrows():
        trans_to_gene[entry["name"]] = entry["name2"]
    return trans_to_gene
    

def output_intron_table(tables_dir,
                        intron_gff_fname,
                        output_dir):
    """
    Output a table of introns. Just adds length
    and gene information to each entry.
    """
    print "Outputting intron table from %s" %(intron_gff_fname)
    output_basename = os.path.basename(intron_gff_fname).rsplit(".", 1)[0]
    utils.make_dir(output_dir)
    output_fname = os.path.join(output_dir, "%s.gff" %(output_basename))
    print "  - Output file: %s" %(output_fname)
    if not os.path.isfile(intron_gff_fname):
        raise Exception, "Cannot find %s" %(intron_gff_fname)
    trans_to_gene = trans_to_gene_from_table(tables_dir)
    table_fname = os.path.join(tables_dir, "ensGene.kgXref.combined.txt")
    table_df = pandas.read_table(table_fname, sep="\t")
    trans_to_gene = trans_to_gene_from_table(tables_dir)
    intron_entries = pybedtools.BedTool(intron_gff_fname)
    output_file = open(output_fname, "w")
    for entry in intron_entries:
        transcripts = entry.attrs["Parent"].split(",")
        genes_str = \
            ",".join([trans_to_gene[trans] for trans in transcripts])
        entry.attrs["gene_id"] = genes_str
        entry.attrs["region_len"] = str(len(entry))
        output_file.write(str(entry))
    output_file.close()
    return output_fname
    
    
def output_utr_table(tables_dir,
                     utr_gff_fname,
                     output_dir,
                     choice_rule="longest"):
    """
    Output a UTR table (one UTR per gene) given a
    UTR GFF file. Possible rules for choosing the
    UTR ('choice_rule'):

      - longest, uses longest UTR
      - shortest, uses shortest UTR

    Outputs a GFF file.
    """
    print "Outputting UTR table from %s" %(utr_gff_fname)
    output_basename = os.path.basename(utr_gff_fname).rsplit(".", 1)[0]
    utils.make_dir(output_dir)
    output_fname = os.path.join(output_dir, "%s.%s.gff" %(output_basename,
                                                          choice_rule))
    print "  - Output file: %s" %(output_fname)
    if not os.path.isfile(utr_gff_fname):
        raise Exception, "Cannot find %s" %(utr_gff_fname)
    # Load table
    table_fname = os.path.join(tables_dir, "ensGene.kgXref.combined.txt")
    table_df = pandas.read_table(table_fname, sep="\t")
    trans_to_gene = {}
    # Map transcripts to genes
    for row, entry in table_df.iterrows():
        trans_to_gene[entry["name"]] = entry["name2"]
    # Mapping from gene ID to a dictionary mapping each
    # UTR to its length
    genes_to_utr_lens = defaultdict(lambda: defaultdict(int))
    print "Computing lengths of UTRs.."
    gff_utrs = pybedtools.BedTool(utr_gff_fname)
    # Compute lengths of UTRs for each gene
    for entry in gff_utrs:
        # Get transcript that UTR belongs to
        trans_id = entry.attrs["Parent"]
        # Get UTR id
        utr_id = entry.attrs["ID"]
        # Get the gene it corresponds to
        gene_id = trans_to_gene[trans_id]
        # Compute length of UTRs
        # Length of UTR
        utr_len = len(entry)
        genes_to_utr_lens[gene_id][utr_id] = utr_len
    # Select UTR for each gene
    gene_to_chosen_utr = {}
    for gene in genes_to_utr_lens:
        all_utrs = genes_to_utr_lens[gene].items()
        utr_lens = [curr_utr[1] for curr_utr in all_utrs]
        if choice_rule == "longest":
            utr_indx = utils.max_item(utr_lens)[0]
            chosen_utr = all_utrs[utr_indx]
            gene_to_chosen_utr[gene] = chosen_utr
        elif choice_rule == "shortest":
            utr_indx = utils.min_item(utr_lens)[0]
            chosen_utr = all_utrs[utr_indx]
            gene_to_chosen_utr[gene] = chosen_utr
        else:
            raise Exception, "Unsupported choice rule %s" %(choice_rule)
    # Now select the relevant entries for outputting. Also
    # add relevant information about genes/length
    gff_utrs = pybedtools.BedTool(utr_gff_fname)
    gff_out = open(output_fname, "w")
    for entry in gff_utrs:
        # Current UTR id
        curr_utr_id = entry.attrs["ID"]
        # Current UTR's transcript
        curr_utr_trans = entry.attrs["Parent"]
        # Get the current UTR's gene
        curr_utr_gene = trans_to_gene[curr_utr_trans]
        # If this UTR is the chosen UTR, output it
        if gene_to_chosen_utr[curr_utr_gene][0] == curr_utr_id:
            # Look up the gene ID it belongs to
            curr_gene_id = trans_to_gene[curr_utr_trans]
            entry.attrs["gene_id"] = curr_gene_id
            entry.attrs["region_len"] = \
                str(gene_to_chosen_utr[curr_utr_gene][1])
            gff_out.write("%s" %(str(entry)))
    gff_out.close()
    return output_fname


def output_table_seqs(table_gff_fname, fi_fname, output_dir):
    """
    Output table sequences to a file.
    """
    print "Outputting sequences from GFF table..."
    print "  - Input GFF table: %s" %(table_gff_fname)
    print "  - Genome FASTA index: %s" %(fi_fname)
    print "  - Output dir: %s" %(output_dir)
    utils.make_dir(output_dir)
    table_basename = os.path.basename(table_gff_fname).rsplit(".", 1)[0]
    output_fname = os.path.join(output_dir, "%s.fa" %(table_basename))
    print "  - Output file: %s" %(output_fname)
    if os.path.isfile(output_fname):
        print "Found %s. Skipping..." %(output_fname)
        return output_fname
    entries = pybedtools.BedTool(table_gff_fname)
    def fields2name(f):
        """
        replace GFF featuretype field with the attributes field.
        """
        #f[2] = f[-1]
        custom_field = "%s:%s-%s:%s" %(f.chrom, f.start, f.stop, f.strand)
        f[2] = "%s;%s" %(custom_field, f[-1])
        return f
    # Output sequences as FASTA
    try:
        entries.each(fields2name).sequence(fi=fi_fname, fo=output_fname,
                                           s=True, name=True)
    except pybedtools.helpers.BEDToolsError as s:
        pass
    return output_fname
    

def main():
    # Load tables from the init dir.
    tables_dir = os.path.expanduser("~/jaen/test/mm9/ucsc/")
    utr_fname = os.path.expanduser("~/jaen/pipeline_init/mm9/ucsc/utrs/ensGene.3p_utrs.gff")
    output_fname = "./testtable.gff"
    #output_utr_table(tables_dir, utr_fname, output_fname)
    #output_fname = "./smallutrs.gff"
    #utrs = pybedtools.BedTool(output_fname)
    fi_fname = os.path.expanduser("~/jaen/test/mm9/genome/mm9.fa")
    output_table_seqs(output_fname, fi_fname, ".")
    


if __name__ == "__main__":
    main()
    
