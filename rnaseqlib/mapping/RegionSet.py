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
        """
        Add region to the class.
        """
        # Parse the region length
        if region not in self.region_lens:
            chrom, start, end, strand = utils.parse_coords(region)
            region_len = end - start + 1
            # Record region length
            self.region_lens[region] = region_len
        # Record number of counts per region
        self.region_counts[region] += 1


    def add_region_by_transcript(self, region, transcript,
                                 region_type=None):
        """
        Add region by transcript.
        """
        # Record mapping from region to type
        if region_type is not None:
            self.region_to_type[region] = region_type
        # Record mapping from transcript to region
        self.transcripts_to_regions[transcript].append(region)
        self.add_region(region)


    def get_transcript_regions(transcript, region_type=None):
        """
        Return all the regions in the transcript by their types.
        Returns their counts and lengths.
        """
        transcript_regions = self.transcripts_to_regions[transcript]
        regions = []
        for region in transcript_regions:
            if (region_type is not None):
                curr_region_type = self.region_to_type[region_type]
                # If we're given a region type and this region is
                # not of the right type, then skip ahead
                if curr_region_type != region_type:
                    continue
            # Accumulate the count and length of each region
            region_count = self.region_counts[region]
            region_len = self.region_lens[region]
            regions.append((region_count, region_len))
        print "REGIONS: ", regions
