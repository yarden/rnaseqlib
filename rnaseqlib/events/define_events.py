##
## Main script for defining AS events from a series of databases
##

import misopy
import misopy.gff_utils as gff_utils
import misopy.Gene as gene_utils

##
## UCSC AltEvents: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownAlt.txt.gz
## Ensembl events:  Biomart API
## Swiss Institute of Bioinformatics events: 
##


def define_events(event_type):
    """
    Define alternative splicing events.

    Output a GFF3 file with event representations and
    related gene annotations/domain in attributes field.
    """
    pass

def define_se():
    """
    Define skipped exons
    """
    pass

def define_multi_se():
    """
    Define multiply skipped exons.
    """
    pass

def define_a3ss():
    """
    A3SS
    """
    pass

def define_a5ss():
    """
    A5SS
    """
    pass

def define_ri():
    """
    Retained introns.
    """
    pass

def define_afe():
    """
    AFE
    """
    pass

def define_ale():
    """
    ALE
    """
    pass

def define_mxe():
    """
    MXE
    """
    pass

def define_tandem_utr():
    """
    Tandem UTRs.
    """
    pass
    

def main():
    pass

if __name__ == "__main__":
    main()
