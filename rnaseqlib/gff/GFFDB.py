##
## GFFDB: Database based on Ryan Dale's gffutils library
##
import os
import sys
import time

import gffutils

import rnaseqlib.gff.GFFGene as gffgene


class GFFDB:
    """
    GFF database.
    """
    def __init__(self, db):
        self.db = gffutils.FeatureDB(db)
        self.num_genes = 0
        self.genes = []


    def __repr__(self):
        return "GFFDB(num_genes=%d)" %(self.num_genes)


    def __str__(self):
        return self.__repr__()


    def iter_by_type(self, feature_type):
        """
        Return all genes.
        """
        for entry in self.db.all():
            if entry.featuretype == feature_type:
                yield entry
        

    def get_gene(self, gene_id):
        """
        Get a gene representation.
        """
        gene = gffgene.GFFGene(gene_id, self.db)
        return gene


    def load_all_genes(self):
        """
        Load all genes into memory as GFFGene objects.
        """
        print "Loading genes from GFF database."
        t1 = time.time()
        self.genes = []
        num_genes = 0
        for gene_rec in self.iter_by_type("gene"):
            self.genes.append(self.get_gene(gene_rec.id))
            num_genes += 1
        t2 = time.time()
        print "Loaded %d genes in %.2f seconds." %(num_genes,
                                                   (t2 - t1))
        self.num_genes = num_genes
