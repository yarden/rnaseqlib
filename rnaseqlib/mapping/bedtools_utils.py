##
## Basic wrappers for bedtools
##
def intersect_bam_with_gff_cmd(bam_filename,
                               gff_filename,
                               output_filename,
                               stranded=False):
    """
    Intersect a BAM file with a GFF filename.

    Requires 'tagBam' to be present.
    """
    if stranded:
        raise Exception, "Not implemented."
