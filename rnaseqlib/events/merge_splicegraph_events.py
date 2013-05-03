##
## Merge SpliceGraph events with new annotation
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.gff_utils as gff_utils
import misopy.Gene as gene_utils

import pandas
import pybedtools

from collections import defaultdict

def merge_se(splicegraph_gff_fname, gff_fname, output_fname,
             coords_diff_cutoff=10):
    """
    Merge skipped exons. Takes as splicegraph GFF filename
    and a new GFF filename.
    """
#    splicegraph_db = gene_utils.load_genes_from_gff(splicegraph_gff_fname)
#    gff_db = load_genes_from_gff()
#    gff_out = open(output_fname, "w")
#    merged_gff_db = gff_utils.Writer(gff_out)

    # Load SpliceGraph skipped exons
    splicegraph_in = pybedtools.BedTool(splicegraph_gff_fname)
    splicegraph_exons = \
        splicegraph_in.filter(lambda x: x.attrs["ID"].endswith(".se"))
    # New annotation's skipped exons
    new_in = pybedtools.BedTool(gff_fname)
    new_exons = new_in.filter(lambda x: x.attrs["ID"].endswith(".se"))
    # Intersect splicegraph exons with new exons
    intersected_gff = splicegraph_exons.intersect(new_exons,
                                                  wao=True,
                                                  s=True)
    print intersected_gff.head()
    # Compile the overlaps for each exon
    exons_to_overlaps = defaultdict(list)
    for exon in intersected_gff:
        curr_overlap = int(exon.fields[-1])
        exons_to_overlaps[exon.attrs["ID"]].append(curr_overlap)
    # If the maximum overlap between the SpliceGraph exon and
    # all exons in the new GFF annotation is LESS than 'coords_diff_cutoff'
    # then keep the SpliceGraph exon in the annotation
    # Name of SpliceGraph event trios to include in merged annotation
    splicegraph_trios_to_add = []
    print "MERGING %s and %s" %(splicegraph_gff_fname, gff_fname)
    for exon_id in exons_to_overlaps:
        if max(exons_to_overlaps[exon_id]) < coords_diff_cutoff:
            trio_id = exon_id.rsplit(".", 1)[0]
            splicegraph_trios_to_add.append(trio_id)
    num_sg_trios = len(splicegraph_trios_to_add)
    print "Added %d trios from SpliceGraph" %(num_sg_trios)
    #gff_out.close()
    

def merge_mxe():
    pass


def merge_a3ss():
    pass


def merge_a5ss():
    pass


def merge_afe():
    pass


def merge_ale():
    pass


def merge_events(genome,
                 event_type,
                 splicegraph_events_dir,
                 new_events_dir,
                 output_dir):
    """
    Merge events.
    """
    sg_gff_fname = os.path.join(splicegraph_events_dir,
                                genome,
                                "%s.%s.gff3" %(event_type, genome))
    if not os.path.isfile(sg_gff_fname):
        raise Exception, "Cannot find %s" %(sg_gff_fname)
    new_gff_fname = os.path.join(new_events_dir,
                                 genome,
                                 "commonshortest",
                                 "%s.gff3" %(event_type))
    if not os.path.isfile(new_gff_fname):
        raise Exception, "Cannot find %s" %(new_gff_fname)
    output_dir = os.path.join(output_dir, genome)
    utils.make_dir(output_dir)
    output_gff_fname = \
        os.path.join(output_dir, "%s.%s.gff3" %(event_type, genome))
    print "Merging %s.." %(event_type)
    if event_type.startswith("SE"):
        merge_se(sg_gff_fname,
                 new_gff_fname,
                 output_gff_fname)
                 

def main():
    splicegraph_events_dir = os.path.expanduser("~/jaen/gff-events")
    new_events_dir = os.path.expanduser("~/jaen/new-gff-events/")
    output_dir = os.path.expanduser("~/jaen/new-gff-events/merged-events/")
    # TODO: add AFE/ALE
    event_types = ["SE"]#, "MXE", "A3SS", "A5SS", "RI"]
    print "Merging events..."
    for genome in ["mm9"]:#["mm9", "hg18", "hg19"]:
        print "Genome: %s" %(genome)
        for event_type in event_types:
            print "Event type: %s" %(event_type)
            merge_events(genome,
                         event_type,
                         splicegraph_events_dir,
                         new_events_dir,
                         output_dir)


if __name__ == "__main__":
    main()
    
    
    
    
