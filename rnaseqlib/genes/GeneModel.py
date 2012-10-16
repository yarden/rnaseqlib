##
## Gene Model
##
class Gene:
    """
    Representation of a gene model.
    """
    def __init__(self, transcripts, chrom, strand,
                 label=None):
        self.transcripts = transcripts
        self.chrom = chrom
        self.strand = strand
        self.label = label
        

    def get_const_exons(self, base_diff=6):
        """
        Get constitutive exons.
        """
        num_trans = len(self.transcripts)
        # If we have only one transcript then all
        # exons are constitutive
        if num_trans == 1:
            return self.transcripts[0].parts
        first_trans_exons = self.transcripts[0].parts
        for exon in first_trans_exons:
            # Compare the first exon of the transcript
            # to all other transcripts' exons
            const_exon = True
            for curr_trans in self.transcripts[1:]:
                start_end_diffs = [(abs(exon.start - curr_exon.start),
                                    abs(exon.end - curr_exon.end)) \
                                    for curr_exon in curr_trans.parts]
                print "start end diffs: "
                print start_end_diffs
                raise Exception
                # If all of the current exon start/end diffs
                # are greater than 'base_diff', then the exon is
                # not considered constitutive
                # ....
                # if all(start_end_diffs, axis=1) > base_diff
                # ...
                const_exon = False


    def __repr__(self):
        return "GeneModel(%s, %s, %s)" %(self.label,
                                         self.chrom,
                                         self.strand)
                


class Transcript:
    """
    Transcript of a gene.
    """
    def __init__(self, parts, chrom, strand,
                 label=None):
        self.gene = None
        self.chrom = chrom
        self.strand = strand
        self.parts = parts
        self.label = None

    def __repr__(self):
        parts_str = ",".join(p.__str__() for p in self.parts)
        return "Transcript(%s)" %(parts_str)
    

class Interval:
    """
    Representation of an interval. (Like an exon.)
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
	assert(self.start <= self.end)
	self.len = self.end - self.start + 1
	assert(self.len >= 1)
	
    def __repr__(self):
        return "Interval([%d, %d])" %(self.start, self.end)
    
    def __eq__(self, interval):
	if interval == None: return False
        return self.start == interval.start and self.end == interval.end

    def __ne__(self, interval):
        return not self.__eq(interval)

    def __lt__(self, interval):
        if self.start == interval.start:
            return self.end < interval.end
        return self.start < interval.start

    def contains(self, start, end):
	if self.start <= start and self.end >= end:
	    return True
	return False

    def intersects(self, other):
        if (self.start < other.end
            and self.end > other.start):
            return True
        return False
        
        
class Part(Interval):
    """
    Part of a transcript.
    """
    def __init__(self, start, end, chrom, strand,
                 label=None):
        Interval.__init__(self, start, end)
        self.chrom = chrom
        self.strand = strand
        self.label = label
        

    def __repr__(self):
        return "Part(%s, %s)" %(str(self.start),
                                str(self.end))
    

    def __str__(self):
        return self.__repr__()

    
        
    
