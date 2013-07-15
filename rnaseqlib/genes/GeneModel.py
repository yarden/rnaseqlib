##
## Gene Model
##
import os
import sys
import time


import numpy
from numpy import *
import numpy as np

import rnaseqlib
import rnaseqlib.utils as utils

import misopy
import misopy.gff_utils as gff_utils

import operator
from collections import namedtuple

class Gene:
    """
    Representation of a gene model.

    By convention, all coordinates will be a 1-based start (GFF convention.)
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
        if not self.has_cds:
            print "has no cds: %s" %(self.label)
        # Get all parts from all transcripts
        self.parts = []
        for trans in self.transcripts:
            self.parts.extend(trans.parts)
        # Get all CDS parts from all transcripts
        self.cds_parts = []
        for trans in self.transcripts:
            self.cds_parts.extend(trans.cds_parts)

            
    def get_inclusive_trans_coords(self):
        """
        Return the most inclusive transcription start/end coordinates,
        i.e. the lowest start coordinate across all transcripts and
        the highest end coordinate across all transcripts.
        """
        trans_starts = [transcript.start for transcript in self.transcripts]
        trans_ends = [transcript.end for transcript in self.transcripts]
        inclusive_coords = [min(trans_starts), max(trans_ends)]
        return inclusive_coords
        

    def get_cds_transcripts(self):
        """
        Return only transcripts with CDS regions.
        """
        cds_transcripts = []
        for transcript in self.transcripts:
            if transcript.has_cds:
                cds_transcripts.append(transcript)
        return cds_transcripts

            
    def get_parts(self, cds_only=False):
        """
        Get all parts from all transcripts.
        """
        seen_parts = {}
        parts = []
        for trans in self.transcripts:
            if cds_only:
                trans_parts = trans.get_cds_parts()
            else:
                trans_parts = trans.parts
            for exon in trans_parts:
                # Skip parts that we've collected
                # already
                if (exon.start, exon.end) in seen_parts:
                    continue
                parts.append(exon)
                seen_parts[(exon.start, exon.end)] = True
        return parts
        
        
    def compute_const_exons(self,
                            base_diff=6,
                            cds_only=False,
                            frac_const=.7,
                            min_exon_space=40):
        """
        Get constitutive exons.

        - base_diff: the number of nucleotides by which
          an exon can differ from exons in a transcript to
          be considered as occurring in that transcript

        - cds_only: use CDS only exons or not

        - frac_const: the fraction of transcripts in which an exon
          must occur to be considered constitutive.  If 1, then the
          exon must occur in all transcripts.  Only used
          when no truly constitutive exons are available.
        """
        self.const_exons = []
        transcripts = []
        if cds_only:
            # If asked for CDS-only but there's no CDS,
            # then quit
            if not self.has_cds:
                return self.const_exons
            transcripts = self.get_cds_transcripts()
        else:
            transcripts = self.transcripts
        num_trans = len(transcripts)
        # If we have only one transcript then all
        # exons are constitutive
        frac_str = "NA"
        if num_trans == 1:
            if cds_only:
                self.const_exons = transcripts[0].get_cds_parts()
            else:
                self.const_exons = transcripts[0].parts
            frac_str = ",".join(len(self.const_exons) * ["1"])
            return self.const_exons, frac_str
        # If there's one or more fully constitutive exons, use these only
        self.const_exons, frac_str = \
            self.get_const_exons_from_transcripts(transcripts,
                                                  frac_const,
                                                  cds_only,
                                                  base_diff,
                                                  min_exon_space=min_exon_space)
        return self.const_exons, frac_str

    
    def get_const_exons_from_transcripts(self,
                                         transcripts,
                                         frac_const,
                                         cds_only,
                                         base_diff,
                                         min_exon_space=40):
        """
        Get (approximately) constitutive exons in set of transcripts.
        Return the constitutive exons and the fraction of
        transcripts they appear in.

        Parameters:
        -----------

        transcripts: set of transcripts
        frac_const: fraction of transcripts an exon must be in to be
        constitutive
        base_diff: base difference allowed when considering an exon
        to be 'included' in transcript
        """
        num_trans = len(transcripts)
        # Iterate through each exon and determine what fraction
        # of the transcripts it appears in.
        exons = self.get_parts(cds_only=cds_only)
        # Fraction of transcripts that each exon occurs in
        exons_fracs = []
        frac_str = "NA"
        if len(exons) == 0:
            return exons, frac_str
        for exon in exons:
            # Calculate if the exon is in or out of each
            # transcript
            exon_status = array([trans.has_part(exon,
                                                base_diff=base_diff,
                                                cds_only=cds_only) \
                                 for trans in transcripts])
            # Compute the fraction of transcripts in which the exon
            # occurs.
            in_transcripts_frac = \
                len(where(exon_status == True)[0]) / float(num_trans)
            exons_fracs.append(in_transcripts_frac)
        ## Get the exons that are most approximately constitutive, i.e.
        ## occur in high fraction of transcripts
        exons_dict = \
            dict([(exon, frac) for exon, frac in zip(exons, exons_fracs)])
        # Sort exons by how constitutive they are
        sorted_exons = \
            sorted(exons_dict.iteritems(), key=operator.itemgetter(1),
                   reverse=True)
        # Get maximally constitutive exons first
        max_frac = max(exons_fracs)
        const_exons = \
            [e[0] for e in sorted_exons if e[1] == max_frac]
        # Fraction of transcripts chosen exons appear in 
        fracs = len(const_exons) * [max_frac]
        # Compute sum of lengths of exons
        sum_lens = np.sum([exon.len for exon in const_exons])
        if sum_lens < min_exon_space:
            print "NOTE: sum of constitutive exon lengths for %s is only %d" \
                  %(self.label, sum_lens)
            # Try to get next exon in list if there are any
            if len(sorted_exons) >= 2:
                next_exon, next_exon_frac = sorted_exons[1]
                print "Fetching next exon which appears in %.2f fraction " \
                      "of transcripts" %(next_exon_frac)
                const_exons.append(next_exon)
                fracs.append(next_exon_frac)
                sum_lens += next_exon.len
        # Return constitutive exons and a string describing the fraction of
        # transcripts they appear in
        frac_str = ",".join(["%.2f" %(f) for f in fracs])
        return const_exons, frac_str


#     def old_compute_const_exons(self,
#                             base_diff=6,
#                             cds_only=False,
#                             frac_const=1):
#         """
#         Get constitutive exons.

#         - base_diff: the number of nucleotides by which
#           an exon can differ from exons in a transcript to
#           be considered as occurring in that transcript

#         - cds_only: use CDS only exons or not

#         - frac_const: the fraction of transcripts in which an exon
#           must occur to be considered constitutive.  If 1, then the
#           exon must occur in all transcripts.
#         """
#         self.const_exons = []
#         transcripts = []
#         if cds_only:
#             # If asked for CDS-only but there's no CDS,
#             # then quit
#             if not self.has_cds:
#                 return self.const_exons
#             transcripts = self.get_cds_transcripts()
#         else:
#             transcripts = self.transcripts
#         num_trans = len(transcripts)
#         # If we have only one transcript then all
#         # exons are constitutive
#         if num_trans == 1:
#             if cds_only:
#                 self.const_exons = transcripts[0].get_cds_parts()
#             else:
#                 self.const_exons = transcripts[0].parts
#             return self.const_exons
# #        print "CDS TRANSCRIPTS and their CDS PARTS: "
# #        if cds_only:
# #            for t in transcripts:
# #                print "=> %s" %(t.label), t.get_cds_parts(), " n = %d" %(len(t.get_cds_parts()))
# #            raise Exception
#         # If we're asked to deal only with CDS exons, take
#         # only CDS exons if first transcript
#         if not cds_only:
#             first_trans_exons = transcripts[0].parts
#         else:
#             first_trans_exons = transcripts[0].get_cds_parts()
#         for exon in first_trans_exons:
#             # Compare the first exon of the transcript
#             # to all other transcripts' exons
#             const_exon = True
#             # Whether we found the CDS or not
#             found_cds = False
#             for curr_trans in transcripts[1:]:
#                 if not cds_only:
#                     # Consider all exons
#                     trans_parts = curr_trans.parts
#                 else:
#                     # Take only exons that fall in the CDS
#                     trans_parts = curr_trans.get_cds_parts()
#                     if len(trans_parts) == 0:
#                         continue
#                     else:
#                         found_cds = True
#                 # Compute the start coord difference and end coord difference
#                 # between our exon and each exon in the current transcript
#                 start_end_diffs = array([(abs(exon.start - curr_exon.start),
#                                           abs(exon.end - curr_exon.end)) \
#                                           for curr_exon in trans_parts])
#                 print "start_end_diffs: "
#                 print start_end_diffs
#                 # The exon is NOT considered constitutive if there are no exons
#                 # in the transcripts whose start/end diff with the current exon
#                 # is less than or equal to 'base_diff'
#                 status = all(start_end_diffs <= base_diff, axis=1)
#                 exons_in_transcripts.append(status)
# #                if all(status == False):
# #                    const_exon = False
# #                    # Determined exon is not constitutive, so skip to next exon
# #                    break
#             # Compute the fraction of transcripts in which the exon
#             # occurs.  Add 1 to account for the fact that the
#             # exon occurs in the current transcript.
#             in_transcripts_frac = \
#                 (len(where(status == True)) + 1) / float(num_trans)
#             print "in transcripts frac: ", in_transcripts_frac

#             if cds_only:
#                 # If we're asked for CDS only exons, check that the CDS was
#                 # found and that the exon is constitutive
#                 if found_cds and const_exon:
#                     self.const_exons.append(exon)
#             elif const_exon:
#                 # If not CDS only, just check that we have constitutive
#                 # exons
#                 self.const_exons.append(exon)
#         return self.const_exons
        

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
        # Start/end of transcript defined by start/end
        # of first and last exons, respectively
        self.start = parts[0].start
        self.end = parts[-1].end
        self.parts = parts
        self.label = label
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.cds_coords = (self.cds_start,
                           self.cds_end)
        self.parent = parent
        self.cds_parts = self.get_cds_parts()
        self.has_cds = False
        if len(self.cds_parts) > 0:
            self.has_cds = True
            
        
    def __repr__(self):
        parts_str = ",".join(p.__str__() for p in self.parts)
        return "Transcript(%s, %s, %s, label=%s, parent=%s, has_cds=%s)" \
            %(parts_str,
              self.chrom,
              self.strand,
              str(self.label),
              self.parent,
              str(self.has_cds))

    
    def has_part(self, part,
                 cds_only=False,
                 base_diff=0):
        """
        Return True if the part exists in the transcript.

        - base_diff: the number of nucleotides by which
          an exon can differ from exons in a transcript to
          be considered as occurring in that transcript

        - cds_only: whether to use CDS only parts of the
          transcript
        """
        trans_parts = []
        if cds_only:
            trans_parts = self.get_cds_parts()
        else:
            trans_parts = self.parts
        if len(trans_parts) == 0:
            # If no parts found, assume that the exon
            # is not present in the transcript
            return False
        # Compute the start coord difference and end coord difference
        # between our exon and each exon in the current transcript
        start_end_diffs = array([(abs(part.start - curr_exon.start),
                                  abs(part.end - curr_exon.end)) \
                                  for curr_exon in trans_parts])
        # The exon is NOT considered constitutive if there are no exons
        # in the transcripts whose start/end diff with the current exon
        # is less than or equal to 'base_diff'
        status = all(start_end_diffs <= base_diff, axis=1)
        if all(status == False):
            return False
        return True

    
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
            cds_part = None
            # Check for special case that CDS is entirely contained
            # within the exon
            if (part.start <= self.cds_start) and \
               (part.end >= self.cds_end):
                # If so, make the part be the CDS itself
                cds_part = Part(self.cds_start, self.cds_end,
                                chrom=part.chrom,
                                strand=part.strand,
                                parent=part.parent)
                cds_part.label = "cds.%s:%s-%s:%s" %(cds_part.chrom,
                                                     str(cds_part.start),
                                                     str(cds_part.end),
                                                     part.strand)
                self.cds_parts = [cds_part]
                break
            elif (part.start <= self.cds_start) and \
                 (part.end > self.cds_start):
                # If the part is overlapping the CDS start make it
                # start at the CDS
                cds_part = Part(self.cds_start, part.end,
                                chrom=part.chrom,
                                strand=part.strand,
                                parent=part.parent)
            elif (part.start <= self.cds_end) and \
                 (part.end > self.cds_end):
                # If the part is overlapping the CDS end make it
                # end at the CDS
                cds_part = Part(part.start, self.cds_end,
                                chrom=part.chrom,
                                strand=part.strand,
                                parent=part.parent)
            elif (part.start >= self.cds_start) and \
                 (part.end <= self.cds_end):
                # If the part is totally contained within CDS,
                # add it as is
                cds_part = Part(part.start, part.end,
                                chrom=part.chrom,
                                strand=part.strand,
                                parent=part.parent)
            # Set label for CDS part to match CDS coordinates
            cds_part.label = "cds.%s:%s-%s:%s" %(cds_part.chrom,
                                                 str(cds_part.start),
                                                 str(cds_part.end),
                                                 part.strand)
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
        self.len = end - start + 1


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


        

    
        
    
