##
## Gene Model
##
class Gene:
    """
    Representation of a gene model.
    """
    def __init__(self, chrom, strand):
        self.exons = []
        self.transcripts = []
        self.chrom = chrom
        self.strand = strand
        

    def get_const_exons(self, base_diff=6):
        """
        Get constitutive exons.
        """
        pass



class Transcript:
    """
    Transcript of a gene.
    """
    def __init__(self, parts, gene):
        self.gene = gene
        self.parts = parts
    

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
        
        
class Part:
    """
    Part of a transcript.
    """
    def __init__(self, start, end, chrom, strand):
        Interval.__init__(self, start, end)
        self.chrom = chrom
        self.strand = strand

    
        
    
