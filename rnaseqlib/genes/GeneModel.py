##
## Gene Model
##
import os
import sys
import time

import numpy
from numpy import *

import misopy
import misopy.gff_utils as gff_utils

from collections import namedtuple

class Gene:
    """
    Representation of a gene model.

    By convention, all coordinates will be a 1-based start (GFF conventions.)
    """
    __slots__ = ['transcripts', 'chrom', 'strand', 'label', 'gene_symbol',
                 'parts', 'const_exons', 'has_cds']
    def __init__(self, transcripts, chrom, strand,
                 label=None,
                 gene_symbol=None):
        self.transcripts = transcripts
        self.chrom = chrom
        self.strand = strand
        self.label = label
        self.gene_symbol = gene_symbol
        self.const_exons = []
        self.has_cds = False
        # Determine if the gene has at least one
        # CDS containing transcript
        for t in transcripts:
            if t.cds_coords is not None:
                self.has_cds = True
                break
        if not self.has_cds: print "has no cds: %s" %(self.label)
        # Get all parts from all transcripts
        self.parts = []
        for trans in self.transcripts:
            self.parts.extend(trans.parts)
        # Get all CDS parts from all transcripts
        self.cds_parts = []
        for trans in self.transcripts:
            self.cds_parts.extend(trans.cds_parts)
        

    def compute_const_exons(self, base_diff=6, cds_only=False):
        """
        Get constitutive exons.
        """
        self.const_exons = []
        num_trans = len(self.transcripts)
        # If we have only one transcript then all
        # exons are constitutive
        if num_trans == 1:
            self.const_exons = self.transcripts[0].parts
            return self.transcripts[0].parts
        first_trans_exons = self.transcripts[0].parts
        # If there's no CDS-containing transcript, then
        # there are no constitutive CDS only exons
        if cds_only and (not self.has_cds):
            return self.const_exons
        for exon in first_trans_exons:
            # Compare the first exon of the transcript
            # to all other transcripts' exons
            const_exon = True
            # Whether we found the CDS or not
            found_cds = False
            for curr_trans in self.transcripts[1:]:
                if not cds_only:
                    # Consider all exons
                    trans_parts = curr_trans.parts
                else:
                    # Take only exons that fall in the CDS
                    trans_parts = curr_trans.get_cds_parts()
                    if len(trans_parts) == 0:
                        continue
                    else:
                        found_cds = True
                # Compute the start coord difference and end coord difference
                # between our exon and each exon in the current transcript
                start_end_diffs = array([(abs(exon.start - curr_exon.start),
                                          abs(exon.end - curr_exon.end)) \
                                          for curr_exon in trans_parts])
                # The exon is NOT considered constitutive if there are no exons
                # in the transcripts whose start/end diff with the current exon
                # is less than or equal to 'base_diff'
                status = all(start_end_diffs <= base_diff, axis=1)
                if all(status == False):
                    const_exon = False
                    # Determined exon is not constitutive, so skip to next exon
                    break
            if cds_only:
                # If we're asked for CDS only exons, check that the CDS was
                # found and that the exon is constitutive
                if found_cds and const_exon:
                    self.const_exons.append(exon)
            elif const_exon:
                # If not CDS only, just check that we have constitutive
                # exons
                self.const_exons.append(exon)
        return self.const_exons
        

    def __repr__(self):
        return "GeneModel(%s, %s, %s)" %(self.label,
                                         self.chrom,
                                         self.strand)
                


class Transcript:
    """
    Transcript of a gene.
    """
    __slots__ = ['parts', 'chrom', 'strand', 'label',
                 'cds_start', 'cds_end', 'cds_coords',
                 'cds_parts', 'parent']
    def __init__(self, parts, chrom, strand,
                 label=None,
                 cds_start=None,
                 cds_end=None,
                 parent=None):
        self.gene = None
        self.chrom = chrom
        self.strand = strand
        self.parts = parts
        self.label = label
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.cds_coords = (self.cds_start,
                           self.cds_end)
        self.parent = parent
        self.cds_parts = self.get_cds_parts()

        
    def __repr__(self):
        parts_str = ",".join(p.__str__() for p in self.parts)
        return "Transcript(%s, %s, %s, label=%s, parent=%s)" %(parts_str,
                                                               self.chrom,
                                                               self.strand,
                                                               str(self.label),
                                                               self.parent)

    def get_cds_parts(self, min_cds_len=10):
        """
        Compute the parts that are in the CDS.  Trim UTR containing
        exons to start/end at the CDS start/end, respectively.

        If the CDS length is less than 'min_cds_len' nucleotides,
        skip it altogether.
        """
        self.cds_parts = []
        cds_len = self.cds_end - self.cds_start + 1
        if (self.cds_start is None) or \
           (cds_len < min_cds_len):
            return self.cds_parts
        # Compute the parts that are in the CDS
        for part in self.parts:
            # Skip parts that end before the CDS or
            # start after the CDS
            if (part.end <= self.cds_start) or \
               (part.start >= self.cds_end):
                continue
            cds_part_label = "cds.%s" %(part.label)
            cds_part = None
            # Check for special case that CDS is entirely contained
            # within the exon
            if (part.start <= self.cds_start) and \
               (part.end >= self.cds_end):
                # If so, make the part be the CDS itself
                cds_part = Part(self.cds_start, self.cds_end,
                                chrom=part.chrom,
                                strand=part.chrom,
                                label=cds_part_label,
                                parent=part.parent)
                self.cds_parts = [cds_part]
                break
            elif (part.start <= self.cds_start) and \
                 (part.end > self.cds_start):
                # If the part is overlapping the CDS start make it
                # start at the CDS
                cds_part = Part(self.cds_start, part.end,
                                chrom=part.chrom,
                                strand=part.chrom,
                                label=cds_part_label,
                                parent=part.parent)
            elif (part.start <= self.cds_end) and \
                 (part.end > self.cds_end):
                # If the part is overlapping the CDS end make it
                # end at the CDS
#                print "CDS PART IS OVERLAPPING THE END, so.."
                cds_part = Part(part.start, self.cds_end,
                                chrom=part.chrom,
                                strand=part.chrom,
                                label=cds_part_label,
                                parent=part.parent)
            elif (part.start >= self.cds_start) and \
                 (part.end <= self.cds_end):
                # If the part is totally contained within CDS,
                # add it as is
#                print "part totally within CDS, =>"
                cds_part = Part(part.start, part.end,
                                chrom=part.chrom,
                                strand=part.chrom,
                                label=cds_part_label,
                                parent=part.parent)
            if cds_part is not None:
                self.cds_parts.append(cds_part)
        self.cds_parts = tuple(self.cds_parts)
        return self.cds_parts
        
    
class Part:
    """
    Part of a transcript.
    """
    __slots__ = ['start', 'end', 'chrom', 'strand', 'label', 'parent']
    def __init__(self, start, end, chrom=None, strand=None,
                 label=None, parent=None):
        self.start = start
        self.end = end
        self.chrom = chrom
        self.strand = strand
        self.label = label
        self.parent = parent
        # If label not given, set default label
        if self.label is None:
            self.label = "%s:%s-%s:%s" %(chrom,
                                         start,
                                         end,
                                         strand)

    def __repr__(self):
        return "Part(%s, %s)" %(str(self.start),
                                str(self.end))
    

    def __str__(self):
        return self.__repr__()


def output_parts_as_gff(gff_out, parts, chrom, strand,
                        source=".",
                        rec_type="exon",
                        gene_id="NA",
                        na_val="NA"):
    """
    Output a set of parts to GFF.
    """
    for part in parts:
        attributes = {
            'ID': ["%s.%s" %(rec_type, part.label)],
            'Parent': [part.parent],
            'gene_id': [gene_id],
            }
        rec_start, rec_end = part.start, part.end
        gff_rec = gff_utils.GFF(chrom, source, rec_type,
                                rec_start, rec_end,
                                attributes=attributes,
                                strand=strand)
        gff_out.write(gff_rec)


##
## Coordinate utilities
##
def part_contained_in(part, coords):
    """
    Return True if the first part is contained
    with in the given coords.
    """
    if (part.start >= coords[0]) and \
       (part.end <= coords[1]):
       return True
    return False


        

    
        
    
