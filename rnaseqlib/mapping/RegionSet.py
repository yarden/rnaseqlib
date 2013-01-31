##
## RegionSet: class for holding regions and their
## lengths
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

class RegionSet:
    """
    Representation of a set of regions, optionally indexed by
    transcripts.  Stores the length of regions and their
    number of occurrences.
    """
    def __init__(self):
        self.regions = None
        # Map from region to their lengths
        self.region_lens = defaultdict(int)
        # Map from region to their counts
        self.region_counts = defaultdict(int)
        # Map from transcripts to regions
        self.transcripts_to_regions = defaultdict(list)
        # Map from region to region types
        self.region_to_type = {}


    def add_region(self, region):
        pass


    def add_region_by_transcript(self, region, transcript,
                                 region_type=None):
        """
        Add region by transcript.
        """
        # Record mapping from transcript to region
        self.transcripts_to_regions[transcript].append(region)
        # Record mapping from region to type
        if region_type is not None:
            self.region_to_type[region] = region_type
        
        pass

    
