##
## CLIP utilities
##
import os
import sys
import time
import subprocess

import pybedtools

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.fastx_utils as fastx_utils
import rnaseqlib.mapping.bedtools_utils as bedtools_utils


def trim_clip_adaptors(fastq_filename,
                       adaptors_filename,
                       output_dir,
                       logger,
                       min_read_len=5):
    """
    Trim CLIP adaptors using 'cutadapt'.
    """
    logger.info("Trimming CLIP adaptors from: %s" %(fastq_filename))
    cutadapt_path = utils.which("cutadapt")
    if cutadapt_path is None:
        logger.critical("Could not find \'cutadapt\' on the path. " \
                        "Please install \'cutadapt\' or make the installed " \
                        "version available on path.")
    output_basename = \
        utils.trim_fastq_ext(os.path.basename(fastq_filename))
    output_filename = os.path.join(output_dir,
                                   "%s_trimmed.fastq.gz" \
                                   %(output_basename))
    if os.path.isfile(output_filename):
        logger.info("SKIPPING: %s already exists!" \
                    %(output_filename))
        return output_filename
    logger.info("  - Outputting trimmed sequences to: %s" \
                %(output_filename))
    # Load adaptors to pass to 'cutadapt'
    if not os.path.isfile(adaptors_filename):
        logger.critical("Could not find adaptors file %s" \
                        %(adaptors_filename))
        sys.exit(1)
    adaptors_in = open(adaptors_filename, "r")
    # Substitute newlines with spaces
    adaptors = adaptors_in.read().strip().replace("\n", " ")
    adaptors_in.close()
    cutadapt_cmd = "%s %s %s -o %s -m %d -q 3 > %s.log" %(cutadapt_path,
                                                          adaptors,
                                                          fastq_filename,
                                                          output_filename,
                                                          min_read_len,
                                                          output_filename)
    logger.info("Executing: %s" %(cutadapt_cmd))
    t1 = time.time()
    os.system(cutadapt_cmd)
    t2 = time.time()
    logger.info("Trimming took %.2f mins." %((t2 - t1)/60.))
    return output_filename


def collapse_clip_reads(sample, output_dir, logger):
    """
    Collapse CLIP reads. Uses fastx_collapser.
    """
    logger.info("Collapsing CLIP reads for %s" %(sample.label))
    t1 = time.time()
    collapsed_seq_filename = \
        fastx_utils.fastx_collapse_fastq(sample.rawdata.reads_filename,
                                         output_dir,
                                         logger)
    if collapsed_seq_filename is None:
        logger.critical("Collapsing of CLIP reads failed.")
        sys.exit(1)
    t2 = time.time()
    logger.info("Collapsing took %.2f minutes." %((t2 - t1)/60.))
    return collapsed_seq_filename


def check_clip_utils(logger,
                     required_utils=["cutadapt",
                                     "fastx_collapser"]):
    """
    Check that necessary utilities are available.
    """
    logger.info("Checking that utilities required for CLIP are available..")
    for program in required_utils:
        program_path = utils.which(program)
        if program_path is None:
            logger.critical("Could not access: %s" %(program))
            logger.critical("Make %s avaialble and try again." %(program))
            sys.exit(1)
    logger.info("Found CLIP utilities.")


def filter_clusters(logger, clusters_bed_fname, output_dir,
                    num_reads=8,
                    depth=0,
                    min_size=12,
                    max_size=100):
    """
    Filter clusters by number of reads in them and/or depth.
    """
    if not clusters_bed_fname.endswith(".bed"):
        logger.critical("Error: clusters filename %s must end in .bed" \
                        %(clusters_bed_fname))
        sys.exit(1)
    bed_basename = \
        os.path.basename(clusters_bed_fname).rsplit(".bed", 1)[0]
    filtered_clusters_fname = \
        os.path.join(output_dir, "%s.num_reads_%d.depth_%d.bed" \
                     %(bed_basename,
                       num_reads,
                       depth))
    num_passing_filter = 0
    with open(filtered_clusters_fname, "w") as clusters_out:
        with open(clusters_bed_fname, "r") as clusters_in:
            for line in clusters_in:
                fields = line.strip().split("\t")
                # BED length does not include the stop
                cluster_len = int(fields[2]) - int(fields[1]) 
                cluster_reads = int(fields[4])
                # Check that it meets the number of reads filter
                if cluster_reads < num_reads:
                    continue
                # Check that it meets the depth filter
                cluster_depth = cluster_reads / float(cluster_len)
                if cluster_depth < depth:
                    continue
                # Check that it meets the cluster size filter
                if cluster_len < min_size:
                    continue
                if cluster_len > max_size:
                    continue
                clusters_out.write(line)
                num_passing_filter += 1
    if num_passing_filter == 0:
        logger.critical("0 clusters pass filter.")
        sys.exit(1)
    return filtered_clusters_fname


def output_clip_clusters(logger,
                         bed_filename,
                         output_filename,
                         genes_gff_fname,
                         cluster_dist=0,
                         skip_junctions=True):
    """
    Output clusters in BED-like format:

      chrom, start, end, cluster_id, score, strand

    where:

      score is the number of reads in the cluster.
      strand is determined by the gene it aligns to.

    In contrast to merge, cluster does not flatten the cluster of
    intervals into a new meta-interval; instead, it assigns an unique
    cluster ID to each record in each cluster.

    Return the output filename or None if it went wrong.
    """
    if cluster_dist != 0:
        raise Exception, "Cluster distance must be 0 for now."
    logger.info("Outputting CLIP clusters...")
    logger.info("  - BED input: %s" %(bed_filename))
    logger.info("  - Output file: %s" %(output_filename))
    logger.info("  - Cluster dist: %d" %(cluster_dist))
    # Find clusters in BED representation of reads:
    #  (1) Run clusterBed
    #  (2) Merge the cluster
    # Cluster the BED 
    clusterBed_cmd = "clusterBed -i %s" %(bed_filename)
    # Parse the resulting clusters, skipping junction reads if asked
    cluster_proc = subprocess.Popen(clusterBed_cmd, shell=True,
                                    stdin=sys.stdin,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    # Open stream to mergeBed command and write clusters to it
    mergeBed_cmd = "mergeBed -i - -nms -scores collapse"
    output_file = open(output_filename, "w")
    merge_proc = subprocess.Popen(mergeBed_cmd, shell=True,
                                  stdin=subprocess.PIPE,
                                  stdout=output_file)
    # Process each cluster line and feed it to mergeBed
    for line in iter(cluster_proc.stdout.readline, ""):
        fields = line.strip().split("\t")
        cigar = fields[6]
        if skip_junctions and ("N" in cigar):
            # Skip junctions if asked
            continue
        # Make the score of BED be the cluster ID (cluster number)
        fields[3], fields[4] = \
            fields[-1], fields[3]
        processed_line = "%s\n" %("\t".join(fields))
        merge_proc.stdin.write(processed_line)
    merge_proc.communicate()
    output_file.close()
    ##
    ## Post process the merged clusters to find gene strand
    ##
    # Rename the current file
    orig_merged_fname = "%s.orig" %(output_filename)
    os.rename(output_filename, orig_merged_fname)
    # Get rid of semi-colon separated cluster names
    with open(orig_merged_fname, "r") as orig_merged:
        with open(output_filename, "w") as output_file:
            for line in orig_merged:
                merged_fields = line.strip().split("\t")
                # Remove semicolon separated read names
                merged_fields[3] = merged_fields[3].split(";")[0]
                merged_line = "%s\n" %("\t".join(merged_fields))
                output_file.write(merged_line)
    # Map clusters to genes
    genes_merged_fname = "%s.genes" %(output_filename)
    os.rename(output_filename, genes_merged_fname)
    output_clip_clusters_with_genes(genes_merged_fname,
                                    genes_gff_fname,
                                    output_filename)
    # Remove temporary files
    os.remove(orig_merged_fname)
    os.remove(genes_merged_fname)
    return output_filename


def output_clip_clusters_with_genes(input_fname,
                                    genes_gff_fname,
                                    output_fname):
    """
    Output clusters with gene annotations.
    """
    def get_gff_gene_id(gff_entry):
        # Keep only ID of GFF gene entry
        gff_entry.attrs = "ID=%s" %(gff_entry.attrs["ID"])
        return gff_entry
    input_bed = pybedtools.BedTool(input_fname)
    # Take gene features only from GFF
    genes_gff = \
        pybedtools.BedTool(genes_gff_fname).filter(lambda x: \
                                                   x.fields[2] == "gene")
    # Keep only gene IDs
    genes_gff_filtered = genes_gff.each(get_gff_gene_id)
    # Intersect clusters with genes GFF annotation
    intersected_bed = input_bed.intersect(genes_gff_filtered, wao=True)
    # Collapse the genes that each cluster aligns to so that
    # you have one line per cluster
    group_cols = [1, 2, 3, 4]
    combined_cols = [5, 9, 10, 12, 14, 15]
    grouped_bed = \
        intersected_bed.groupby(g=group_cols,
                                c=combined_cols,
                                ops=["collapse"] * len(combined_cols))
    output_file = open(output_fname, "w")
    for line in grouped_bed:
        chrom, start, end = \
            line.fields[0], line.fields[1], line.fields[2]
        cluster_id = line.fields[3]
        read_ids = line.fields[4]
        # Compute number of reads in cluster
        num_reads = len(read_ids.split(","))
        target_chrom = line.fields[5]
        gene_ids = "."
        target_strands = line.fields[7]
        if target_strands == ".":
            # No gene could be assigned for cluster, so
            # assume + strand
            target_strands = "+"
        else:
            gff_gene_fields = line.fields[8].split(",")
            genes = []
            for gff_attrs in gff_gene_fields:
                curr_gene = dict([tuple(gff_gene.split("=")) \
                                  for gff_gene in gff_attrs.split(";")])
                genes.append(curr_gene["ID"])
            gene_ids = ",".join(genes)
        # Assign strands: if there are multiple ones,
        # just take the first
        assigned_strand = target_strands.split(",")[0]
        # Output clusters
        cluster_fields = [chrom,
                          start,
                          end,
                          cluster_id,
                          num_reads,
                          assigned_strand,
                          gene_ids]
        cluster_line = "%s\n" %("\t".join(map(str, cluster_fields)))
        output_file.write(cluster_line)
    output_file.close()


# $ grep gene ~/jaen/test/mm9/ucsc/ensGene.gff3 | head -n 200000 | intersectBed -a ~/jaen/Musashi-seq/clip-pipeline-output/analysis/clusters/clip_KH2MSI1_NoDox/accepted_hits.ribosub.sorted.clusters.bed -b stdin -wao | head -n 500 | groupBy -g 1,2,3,4 -c 9,10,12,14,15 -o collapse,collapse,collapse,collapse,collapse | grep 6266729
# chr1	6202501	6202530	317	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29
# chr1	6204921	6204950	318	6204693,6196278	6205373,6266729	-,+	Name=ENSMUSG00000073741;ID=ENSMUSG00000073741;Alias=ENSMUSG00000073741,Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29,29
# chr1	6205260	6205289	319	6204693,6196278	6205373,6266729	-,+	Name=ENSMUSG00000073741;ID=ENSMUSG00000073741;Alias=ENSMUSG00000073741,Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29,29
# chr1	6220114	6220143	320	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29
# chr1	6230757	6230786	321	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29
# chr1	6232787	6232816	322	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29
# chr1	6233566	6233607	323	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	41
# chr1	6233875	6233904	324	6196278	6266729	+	Name=ENSMUSG00000025907;ID=ENSMUSG00000025907;Alias=ENSMUSG00000025907	29


def intersect_clusters_with_gff(logger,
                                clusters_fname,
                                gff_filenames,
                                output_dir):
    """
    Intersect CLIP clusters with a set of GFF filenames (e.g. GFFs containing
    events.
    """
    event_clusters_fnames = []
    logger.info("Intersecting clusters %s with GFFs.." \
                %(clusters_fname))
    for gff_fname in gff_filenames:
        if not os.path.isfile(gff_fname):
            logger.critical("Cannot find events GFF file %s" \
                            %(gff_fname))
            continue
        gff_label = os.path.basename(utils.trim_gff_ext(gff_fname))
        event_clusters_fname = \
            os.path.join(output_dir,
                         "%s.clusters.bed" %(gff_label))
        logger.info("  - Processing GFF %s" %(gff_label))
        # Intersect the clusters for this sample with each
        # GFF file
        intersectBed_cmd = \
            "%s -a %s -b %s -loj -f 1 > %s" \
            %(bedtools_utils.intersectBed_path,
              clusters_fname,
              gff_fname,
              event_clusters_fname)
        if os.path.isfile(event_clusters_fname):
            logger.info("Found %s, skipping.." %(event_clusters_fname))
            event_clusters_fnames.append(event_clusters_fname)
        logger.info("  - Executing: %s" %(intersectBed_cmd))
        ret_val = os.system(intersectBed_cmd)
        if ret_val != 0:
            logger.critical("Cannot intersect clusters with %s" \
                            %(gff_fname))
        else:
            event_clusters_fnames.append(event_clusters_fname)
    return event_clusters_fnames

