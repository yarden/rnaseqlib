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
        self.load_all_genes()


    def __repr__(self):
        return "GFFDB(num_genes=%d)" %(self.num_genes)


    def __str__(self):
        return self.__repr__()


    def iter_recs(self):
        """
        Iterate through all records in gene-centric way.

        Returns records that include the 'gene_id' attribute.
        """
        for gene in self.genes:
            # Return gene record
            gene_rec = self.db[gene.gene_id]
            gene_rec.attributes["gene_id"] = gene.gene_id
            yield gene_rec
            # Return each mRNA
            for mRNA in gene.get_mRNAs():
                yield mRNA
                # Return each mRNA's parts
                for part in gene.get_mRNA_parts(mRNA.id):
                    yield part


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


    def get_gene_objects(self):
        """
        Load expensive representation of genes as OrderedDicts.
        """
        self.gene_objects = []
        for gene in self.genes:
            yield gene.make_gene_object()
