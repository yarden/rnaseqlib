##
## Unit testing for rnaseqlib
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.events as events
import rnaseqlib.events.defineEvents as def_events
import rnaseqlib.events.SpliceGraph as sgraph

import gffutils

CURR_DIR = os.path.dirname(os.path.abspath(__file__))

def data_filename(fn):
    return os.path.join(CURR_DIR, "test_data", fn)


def test_old_ri():
    table_fnames = \
        [data_filename("ri-test/ensGene.txt")]
#         data_filename("ri-test/knownGene.txt"),
#         data_filename("ri-test/refGene.txt")]
    DtoA_F, AtoD_F, DtoA_R, AtoD_R = \
        def_events.prepareSplicegraph(*table_fnames)
    ri_fname = data_filename("ri-test/RI.gff3")
    print "Making RI..."
    def_events.RI(DtoA_F, AtoD_F, DtoA_R, AtoD_R, ri_fname,
                  multi_iso=False)

def test_ri():
    table_fname = \
        os.path.expanduser("./test_data/ri-test/ensGene.txt")
    tables = [table_fname]
    sg = sgraph.SpliceGraph(tables)
    gff_out = gffutils.gffwriter.GFFWriter("./ri_test.gff3")
    sgraph.define_RI(sg, gff_out)
    gff_out.close()


def test_afe():
    pass

def test_afe():
    pass
            

def main():
    test_ri()
    #test_old_ri()


if __name__ == "__main__":
    main()
