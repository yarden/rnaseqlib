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


    def get_gene_rec(self):
        gene_rec = self.db[self.gene_id]
        gene_rec.attributes["gene_id"] = self.gene_id
        return gene_rec
        

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


    def make_gene_object(self):
        """
        Make a representation of a gene.
        Too expensive; using iterator approach instead.
        """
        mRNAs_to_parts = OrderedDict()
        for mRNA in self.get_mRNAs():
            # Exons of the mRNA
            mRNAs_to_parts[mRNA.id] = {"parts": self.get_mRNA_parts(mRNA.id),
                                       "record": mRNA}
        gene_obj = {"gene_rec": self.get_gene_rec(),
                    "mRNAs": mRNAs_to_parts}
        return gene_obj
