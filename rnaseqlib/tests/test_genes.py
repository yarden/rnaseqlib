##
## Unit testing for gene models and calculation of features of genes
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.tests
import rnaseqlib.tests.test_utils as test_utils
import rnaseqlib.genes.GeneModel as gene_model
import rnaseqlib.tables as tables


class TestGenes:
    """
    Test features of gene models.
    """
    def __init__(self):
        pass


    def test_const_exons(self):
        """
        Test calculation of constitutive exons
        """
        print "Testing constitutive exons"
        # Load gene table
        gt = tables.GeneTable(os.path.join(test_utils.TESTDIR, "hg19"), "ensGene")
        ensGene_fname = \
            test_utils.load_test_data(os.path.join("hg19",
                                                   "ensGene.hg19.ENSG00000153944.txt"))
        # Load ensembl genes table into GeneModel objects
        gt.get_ensGene_by_genes(ensGene_fname)
        print "Loaded genes: "
        for g in gt.genes:
            print g, " => ", gt.genes[g]
        assert "ENSG00000153944" in gt.genes, "Could not load ENSG00000153944"
        gene = gt.genes["ENSG00000153944"]
        print "Computing constitutive exons..."
        const_exons = gene.compute_const_exons(base_diff=6, frac_const=0.3)
        print "  - const_exons: ", const_exons
        assert len(const_exons) > 0, "Could not find constitutive exons!"
        print "Computing CDS constitutive exons..."
        cds_const_exons = gene.compute_const_exons(cds_only=True)
        print "  - cds_const_exons: ", cds_const_exons
        

if __name__ == "__main__":
    test_g = TestGenes()
    test_g.test_const_exons()
