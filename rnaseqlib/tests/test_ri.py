##
## Unit testing for rnaseqlib
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.events.defineEvents as def_events


CURR_DIR = os.path.dirname(os.path.abspath(__file__))

def data_filename(fn):
    return os.path.join(CURR_DIR, "test-data", fn)


def test_ri():
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
            

def main():
    test_ri()


if __name__ == "__main__":
    main()
