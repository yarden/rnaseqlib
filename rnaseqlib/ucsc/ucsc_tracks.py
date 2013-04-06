##
## Utilities for making UCSC tracks
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

import pybedtools


def make_ucsc_trackfile(bigWig_urls, output_fname,
                        track_type="bigWig"):
    """
    Make a bigWig track file.

    - bigWig_urls: list of tuples corresponding to the bigWig label
      and its URL
    - output_fname: Output filename.
    - type: type of track (e.g. bigWig, bedGraph, etc.)

   Example track line:
   
   track type=bigWig name="Example One" description="A bigWig file" bigDataUrl=http://genome.ucsc.edu/goldenPath/help/examples/bigWigExample.bw
    """
    rgb_colors = [ \
                   # black
                   "0,0,0",
                   # red
                   "255,0,0",
                   # green
                   "0,255,0",
                   # blue
                   "0,0,255"]
    with open(output_fname) as output_file:
        for bigWig_label, bigWig_url in bigWig_urls:
            if not bigWig_url.lower().startswith("http://"):
                print "Error: URLs must start with http://"
                sys.exit(1)
            track_line = \
                "track type=%s name=\"%s\" description=\"%s\" bigDataUrl=%s\n" \
                %(track_type,
                  bigWig_label,
                  bigWig_label,
                  bigWig_url)
            output_file.write(track_line)
