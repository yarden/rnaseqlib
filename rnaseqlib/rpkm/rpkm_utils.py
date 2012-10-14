##
## Utilities for computing RPKM
##

def output_rpkm(sample, gff_filename, output_dir):
    """
    Output RPKM table per sample.
    """
    pass
    
    
def rpkm_from_gff_aligned_bam(bam_filename, 
                              num_total_reads,
                              output_dir):
    """
    Given a BAM file aligned by bedtools (with 'gff' field),
    compute RPKM for each region, incorporating relevant
    optional fields from gff.
    """
    bam_file = pysam.Samfile(bam_filename, "rb")

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    output_filename = os.path.join(output_dir,
                                   "%s.rpkm.gff" \
                                   %(os.path.basename(bam_filename)))

    source = os.path.basename(output_filename)

    print "Computing RPKM from BAM aligned to GFF..."
    print "  - BAM: %s" %(bam_filename)
    print "  - Output dir: %s" %(output_dir)

    loaded_gff = False
    ref_gff_recs = None

    last_chrom = None

    # Map of gff region to read counts
    region_to_count = defaultdict(int)
    region_to_len = defaultdict(int)

    num_reads = 0

    for bam_read in bam_file:
        curr_chrom = bam_file.getrname(bam_read.tid)
        
        # Load the GFF records for this chromosome if not already loaded
        # if last_chrom != curr_chrom:
        #     print "Now on chrom %s" %(curr_chrom)
        #     ref_gff_recs = load_gff_recs_by_chrom(gff_file,
        #                                           curr_chrom)
        try:
            # Read aligns to region of interest
            gff_aligned_regions = bam_read.opt("YB")
            parsed_regions = gff_aligned_regions.split("gff:")[1:]

            # Compile region counts and lengths
            for region in parsed_regions:
                region_to_count[region] += 1

                # Get region length
                coord_field = region.split(",")[0].split(":")[1]
                region_start, region_end = coord_field.split("-")
                region_start, region_end = int(region_start), \
                                           int(region_end)
                region_len = region_end - region_start
                region_to_len[region] = region_len

                if num_reads % 500000 == 0:
                    print "through %d reads.." %(num_reads)
                num_reads += 1
        except KeyError:
            gff_aligned_region = None
        last_chrom = curr_chrom

    output_rpkm_as_gff(source, region_to_count, region_to_len,
                       num_total_reads, output_filename)
