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
import glob

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.cluster_utils.cluster as cluster
import rnaseqlib.fasta_utils as fasta_utils

from rnaseqlib.init.genome_urls import *
import rnaseqlib.init.download_utils as download_utils


##
## Mapping from organisms to misc. sequence NCBI
## accession IDs.
##
NCBI_MISC_SEQS = {"human": {"chrRibo": "U13369.1",
                            "chrMito": None},
                  "mouse": {"chrRibo": "BK000964.1",
                            "chrMito": None}}

def download_ncbi_fasta(access_id, output_dir):
    """
    Download NCBI FASTA file by accession number and
    label them as access.fasta in the given output directory.
    """
    ncbi_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text" \
        %(access_id)
    url_filename = download_utils.download_url(ncbi_url,
                                               output_dir,
                                               basename="%s.fa" %(access_id),
                                               binary=False)
    return url_filename
    

def download_genome_seq(genome,
                        output_dir):
    """
    Download genome sequence files from UCSC.
    """
    print "Downloading genome sequence files for %s" %(genome)
    print "  - Output dir: %s" %(output_dir)
    output_dir = utils.pathify(os.path.join(output_dir, "genome"))
    utils.make_dir(output_dir)
    dir_files = os.listdir(output_dir)
    # Change to output directory
    os.chdir(output_dir)
    ##
    ## Download the genome sequence files
    ##
    genome_url = "%s/%s/chromosomes/" %(UCSC_GOLDENPATH_FTP,
                                        genome)
    # Fetch all chromosome sequence files
    if len(dir_files) >= 1:
        print "Directory %s exists and contains files; " \
              "skipping download of genome..." \
              %(output_dir)
    else:
        download_utils.wget(os.path.join(genome_url, "*"))
        # Remove random chromosome contigs
        for fname in glob.glob(os.path.join(output_dir, "*.fa.gz")):
            if "_" in os.path.basename(fname):
                print "Deleting: %s" %(fname)
                os.remove(fname)
        ##
        ## Uncompress the files
        ##
        print "Uncompressing files..."
        uncompress_cmd = "gunzip %s/*.gz" %(output_dir)
        print "  - Uncompress cmd: %s" %(uncompress_cmd)
        t1 = time.time()
        ret_val = os.system(uncompress_cmd)
        if ret_val != 0:
            print "Error: Cannot uncompress files in %s" %(output_dir)
            sys.exit(1)
        t2 = time.time()
        print "Uncompressing took %.2f minutes" %((t2 - t1)/60.)
    # Create a single genome FASTA file by concatenating the
    # chromosomes together
    genome_output_fname = \
        os.path.join(output_dir, "%s.fa" %(genome))
    if not os.path.isfile(genome_output_fname):
        print "Concatenating genome chromosomes into one file..."
        print "  - Output file: %s" %(genome_output_fname)
        t1 = time.time()
        concat_chrom_cmd = "cat %s/*.fa > %s" %(output_dir,
                                                genome_output_fname)
        print "  - Concat cmd: %s" %(concat_chrom_cmd)
        ret_val = os.system(concat_chrom_cmd)
        if ret_val != 0:
            print "Error: Could not concatenate genome chromosomes."
            sys.exit(1)
        # Create an index for resulting genome file
        print "Indexing genome file..."
        samtools_index_cmd = "samtools faidx %s" %(genome_output_fname)
        print "  - Index cmd: %s" %(samtools_index_cmd)
        ret_val = os.system(samtools_index_cmd)
        if ret_val != 0:
            print "Error: Could not index genome file."
            sys.exit(1)
        t2 = time.time()
        print "Concatenation and indexing took %.2f minutes" \
            %((t2 - t1)/60.)


def download_misc_seqs(genome, output_dir):
    """
    Download assorted sequences related to genome.
    """
    # Mapping from sequence label (e.g. rRNA)
    # to accession numbers
    organism = None
    if genome.startswith("hg"):
        organism = "human"
    elif genome.startswith("mm"):
        organism = "mouse"
    else:
        print "Error: Unsupported genome."
        sys.exit(1)
    # Fetch the accession numbers for the organism's
    # misc sequences and download them
    misc_seqs = NCBI_MISC_SEQS[organism]
    ncbi_outdir = os.path.join(output_dir, "ncbi")
    misc_outdir = os.path.join(output_dir, "misc")
    utils.make_dir(ncbi_outdir)
    utils.make_dir(misc_outdir)
    for seq_label, access_id in misc_seqs.iteritems():
        if access_id is None:
            continue
        output_filename = os.path.join(misc_outdir, "%s.fa" %(seq_label))
        if os.path.isfile(output_filename):
            print "%s exists. Skipping download.." %(seq_label)
            continue
        print "Downloading: %s (NCBI: %s)" %(seq_label,
                                             access_id)
        url_filename = download_ncbi_fasta(access_id, ncbi_outdir)
        fasta_in = fasta_utils.read_fasta(url_filename)
        fasta_out = open(output_filename, "w")
        print "  - Writing to: %s" %(output_filename)
        # Fetch first FASTA record
        rec = fasta_in.next()
        curr_label, fasta_seq = rec
        # Output it with the required label
        new_rec = (">%s" %(seq_label), fasta_seq)
        fasta_utils.write_fasta(fasta_out, [new_rec])
