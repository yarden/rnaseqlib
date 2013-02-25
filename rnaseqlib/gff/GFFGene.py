##
## Representation of a GFF gene for Ryan Dale's gffutils
##
import gffutils

from collections import OrderedDict


class GFFGene:
    """
    Representation of GFF gene.
    """
    def __init__(self, gene_id, db):
        self.gene_id = gene_id
        self.db = db
        # Map from mRNAs to their children
        #self.mRNAs_to_exons = OrderedDict()
        

    def get_mRNAs(self):
        """
        Get all the mRNA records of gene.
        """
        for mRNA in self.db.children(self.gene_id, level=1):
            mRNA.attributes["gene_id"] = self.gene_id
            yield mRNA


    def get_parts(self):
        """
        Get all parts (e.g. exon/intron entries) of gene.
        """
        for parts in self.db.children(self.gene_id, level=2):
            part.attributes["gene_id"] = self.gene_id
            yield part


    def get_mRNA_parts(self, mRNA_id):
        """
        Get all the parts (exon/intron entries) of a particular
        mRNA by its ID.
        """
        for mRNA_part in self.db.children(mRNA_id):
            gene_parent = self.db.parents(mRNA_id).next()
            mRNA_part.attributes["gene_id"] = gene_parent.id
            yield mRNA_part


    def make_gene(self):
        """
        Make a representation of a gene.
        Too expensive; using iterator approach instead.
        """
        pass
        #self.recs = [r for r in self.db.children(self.gene_id, level=1)]
        #for mRNA in self.db.children(self.gene_id, level=1):
        #    # Exons of the mRNA
        #    self.mRNAs_to_exons[mRNA.id] = \
        #        [e for e in self.db.children(mRNA)]
