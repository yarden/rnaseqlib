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
        self.mRNAs_to_exons = OrderedDict()
        self.make_gene()


    def make_gene(self):
        """
        Make a representation of a gene.
        """
        for mRNA in self.db.children(self.gene_id, level=1):
            # Exons of the mRNA
            self.mRNAs_to_exons[mRNA.id] = \
                [e for e in self.db.children(mRNA)]
