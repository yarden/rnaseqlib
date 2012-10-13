##
## Downloading of genomic sequences related to pipeline
## 

# Basic Downloading
#
# efetch.fcgi?db=<database>&id=<uid_list>&rettype=<retrieval_type>
# &retmode=<retrieval_mode>
#
# Input: List of UIDs (&id); Entrez database (&db); Retrieval type (&rettype); Retrieval mode (&retmode)
# Output: Formatted data records as specified
# Example: Download nuccore GIs 34577062 and 24475906 in FASTA format
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=34577062,24475906&rettype=fasta&retmode=text

import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.cluster_utils.cluster as cluster

from rnaseqlib.init.genome_urls import *
import rnaseqlib.init.download_utils as download_utils


def download_genome_seq(genome,
                        output_dir):
    """
    Download genome sequence files from UCSC.
    """
    print "Downloading genome sequence files for %s" %(genome)
    print "  - Output dir: %s" %(output_dir)
    output_dir = os.path.join(output_dir, "genome")
    if os.path.isdir(output_dir):
        print "Directory %s exists; skipping downloading of genome." \
            %(output_dir)
        return None
    utils.make_dir(output_dir)
    # Change to output directory
    os.chdir(output_dir)
    ##
    ## Download the genome sequence files
    ##
    genome_url = "%s/%s/chromosomes/" %(UCSC_GOLDENPATH,
                                        genome)
    # Fetch all chromosome sequence files
    print "Downloading genome files from: %s" %(genome_url)
    download_utils.wget(download_cmd)
    print "Downloading took %.2f minutes" %((t2 - t1)/60.)
    ##
    ## Uncompress the files
    ##
    print "Uncompressing files..."
    uncompress_cmd = "gunzip %s/*.gz" %(output_dir)
    t1 = time.time()
    os.system(uncompress_cmd)
    t2 = time.time()
    print "Uncompressing took %.2f minutes" %((t2 - t1)/60.)



