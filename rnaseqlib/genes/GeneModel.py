

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
    Transcript class.
    """
    def __init__(self, parts, gene):
        self.gene = gene
        self.parts = parts
    

class Part:
    """
    Transcript part.
    """
    def __init__(self, start, end, chrom, strand):
        self.start = start
        self.end = end
        self.chrom = chrom
        self.strand = strand

    
        
    
