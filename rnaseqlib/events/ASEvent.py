##
## Class defining an alternative splicing event
##

import misopy
import misopy.gff_utils as gff_utils
import misopy.Gene as gene_utils

class Event:
    """
    Simple class defining an alternative splicing event
    """
    def __init__(self):
        self.event_type = None
        # GFF entry corresponding to event
        self.gff_entry = None
        # Gene object (from MISO) that represents the event
        self.gene_obj = None
        # Ensembl gene corresponding to event
        self.ensembl_id = None
        # Supporting evidence / transcripts / ESTs
        # for the event
        self.evidence = None
