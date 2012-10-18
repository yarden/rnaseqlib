##
## Gene Model
##
import os
import sys
import time

import numpy
from numpy import *

from collections import namedtuple

class Gene:
    """
    Representation of a gene model.

    By convention, all coordinates will be a 1-based start (GFF conventions.)
    """
    __slots__ = ['transcripts', 'chrom', 'strand', 'label', 'gene_symbol',
                 'const_exons', 'has_cds']
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
        

    def compute_const_exons(self, base_diff=6, cds_only=False):
        """
        Get constitutive exons.
        """
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
            for curr_trans in self.transcripts[1:]:
                if not cds_only:
                    # Consider all exons
                    trans_parts = curr_trans.parts
                else:
                    # Take only exons that fall in the CDS
                    trans_parts = curr_trans.get_cds_parts()
                    if trans_parts == ():
                        continue
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
            if const_exon:
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
        self.label = None
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.cds_coords = (self.cds_start,
                           self.cds_end)
        self.cds_parts = None
        self.parent = parent

    def __repr__(self):
        parts_str = ",".join(p.__str__() for p in self.parts)
        return "Transcript(%s, %s, %s, parent=%s)" %(parts_str,
                                                     self.chrom,
                                                     self.strand,
                                                     self.parent)

    def get_cds_parts(self):
        # Compute the parts that are in the CDS
        self.cds_parts = tuple(filter(lambda p: \
                                      part_contained_in(p, self.cds_coords),
                                      self.parts))
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


        

    
        
    
