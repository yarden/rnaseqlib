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
    """
    __slots__ = ['transcripts', 'chrom', 'strand', 'label', 'gene_symbol',
                 'const_exons']
    def __init__(self, transcripts, chrom, strand,
                 label=None,
                 gene_symbol=None):
        self.transcripts = transcripts
        self.chrom = chrom
        self.strand = strand
        self.label = label
        self.gene_symbol = gene_symbol
        self.const_exons = []
        

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
        for exon in first_trans_exons:
            # Compare the first exon of the transcript
            # to all other transcripts' exons
            const_exon = True
            for curr_trans in self.transcripts[1:]:
                # Compute the start coord difference and end coord difference
                # between our exon and each exon in the current transcript
                start_end_diffs = array([(abs(exon.start - curr_exon.start),
                                          abs(exon.end - curr_exon.end)) \
                                          for curr_exon in curr_trans.parts])
                # The exon is NOT considered constitutive if there are no exons
                # in the transcripts whose start/end diff with the current exon
                # is less than or equal to 'base_diff'
                status = all(start_end_diffs <= base_diff, axis=1)
                if all(status == False):
                    const_exon = False
                    # Determined exon is not constitutive, so skip to next exon                    
                    break
            # Exon is constitutive
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
                 'cdsStart', 'cdsEnd', 'parent']
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
        self.cds_start = None
        self.cds_end = None
        self.parent = parent

    def __repr__(self):
        parts_str = ",".join(p.__str__() for p in self.parts)
        return "Transcript(%s, %s, %s, parent=%s)" %(parts_str,
                                                     self.chrom,
                                                     self.strand,
                                                     self.parent)
    
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
            self.label = "%s_%s_%s_%s" %(start, end,
                                         chrom, strand)
        

    def __repr__(self):
        return "Part(%s, %s)" %(str(self.start),
                                str(self.end))
    

    def __str__(self):
        return self.__repr__()

    
        
    
