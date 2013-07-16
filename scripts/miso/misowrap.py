##
## misowrap: a wrapper to run MISO on a set of samples
## and processing its output.
##
import os
import sys
import csv
import time
import glob
import itertools
from collections import defaultdict

import pandas

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.miso.misowrap_settings as misowrap_settings
import rnaseqlib.miso.PsiTable as pt
import rnaseqlib.miso.MISOWrap as mw
import rnaseqlib.miso.miso_utils as miso_utils
import rnaseqlib.cluster_utils.cluster as cluster
import rnaseqlib.pandas_utils as pandas_utils

import argh
from argcomplete.completers import EnvironCompleter
from argh import arg


@arg("settings", help="misowrap settings filename.")
@arg("logs-outdir", help="Directory where to place logs.")
@arg("--delay", help="Delay between execution of cluster jobs")
@arg("--dry-run", help="Dry run: do not submit or execute jobs.")
def summarize(settings,
              logs_outdir,
              delay=5,
              dry_run=False):
    """
    Summarize samples in MISO directory.
    """
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="summarize")
    bam_files = misowrap_obj.bam_files
    sample_labels = misowrap_obj.sample_labels
    print "Summarizing MISO output..."
    print "  - Output dir: %s" %(output_dir)
    for sample_label in sample_labels:
        sample_basename = sample_label[0]
        sample_dir_path = \
            utils.pathify(os.path.join(misowrap_obj.miso_outdir,
                                       sample_basename))
        print "Processing: %s" %(sample_basename)
        if not os.path.isdir(sample_dir_path):
            print "Skipping non-directory: %s" %(sample_dir_path)
        # List all event directories in the sample
        event_dirs = os.listdir(sample_dir_path)
        for event_dirname in event_dirs:
            event_dir_path = utils.pathify(os.path.join(sample_dir_path,
                                                        event_dirname))
            if not os.path.isdir(event_dir_path):
                print "Skipping non-dir: %s" %(event_dir_path)
            print "Processing event type: %s" %(event_dirname)
            summary_cmd = \
                "%s --summarize-samples %s %s --summary-label %s" \
                %(misowrap_obj.summarize_miso_cmd,
                  event_dir_path,
                  event_dir_path,
                  sample_basename)
            job_name = "summarize_%s_%s" %(sample_basename,
                                           os.path.basename(event_dirname))
            print "Executing: %s" %(summary_cmd)
            if misowrap_obj.use_cluster:
                misowrap_obj.my_cluster.launch_job(summary_cmd,
                                                   job_name,
                                                   ppn=1)
            else:
                os.system(summary_cmd)
            

###
### Add @ headers here
###
@arg("settings", help="misowrap settings filename.")
@arg("logs-outdir", help="Directory where to place logs.")
@arg("--delay", help="Delay between execution of cluster jobs")
@arg("--dry-run", help="Dry run: do not submit or execute jobs.")
def compare(settings,
            logs_outdir,
            delay=5,
            dry_run=False):
    """
    Run a MISO samples comparison between all pairs of samples.
    """
    settings_filename = utils.pathify(settings)
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="compare")
    bam_files = misowrap_obj.bam_files
    sample_labels = misowrap_obj.sample_labels
    read_len = misowrap_obj.read_len
    overhang_len = misowrap_obj.overhang_len
    miso_bin_dir = misowrap_obj.miso_bin_dir
    miso_output_dir = misowrap_obj.miso_outdir
    comparison_groups = misowrap_obj.comparison_groups
    comparisons_dir = misowrap_obj.comparisons_dir
    utils.make_dir(comparisons_dir)
    misowrap_obj.logger.info("Running MISO comparisons...")
    ##
    ## Compute comparisons between all pairs
    ## in a sample group
    ##
    for comp_group in comparison_groups:
        sample_pairs = utils.get_pairwise_comparisons(comp_group)
        print "  - Total of %d comparisons" %(len(sample_pairs))
        for sample1, sample2 in sample_pairs:
            # For each pair of samples, compare their output
            # along each event type
            misowrap_obj.logger.info("Comparing %s %s" %(sample1,
                                                         sample2))
            # Directories for each sample
            sample1_dir = os.path.join(miso_output_dir,
                                       sample1)
            sample2_dir = os.path.join(miso_output_dir,
                                       sample2)
            for event_type in misowrap_obj.event_types:
                sample1_event_dir = os.path.join(sample1_dir,
                                                 event_type)
                sample2_event_dir = os.path.join(sample2_dir,
                                                 event_type)
                job_name = "compare_%s_%s_%s" %(sample1,
                                                sample2,
                                                event_type)
                event_comparisons_dir = \
                    os.path.join(comparisons_dir,
                                 event_type)
                compare_cmd = "%s --compare-samples %s %s %s " \
                    "--comparison-labels %s %s" \
                    %(misowrap_obj.compare_miso_cmd,
                      sample1_event_dir,
                      sample2_event_dir,
                      event_comparisons_dir,
                      sample1,
                      sample2)
                misowrap_obj.logger.info("Executing: %s" %(compare_cmd))
                if misowrap_obj.use_cluster:
                    if not dry_run:
                        misowrap_obj.my_cluster.launch_job(compare_cmd,
                                                           job_name,
                                                           ppn=1)
                        time.sleep(delay)
                else:
                    if not dry_run:
                        os.system(compare_cmd)


@arg("settings", help="misowrap settings filename.")
@arg("logs-outdir", help="Directory where to place logs.")
@arg("--use-cluster", help="Use cluster to submit jobs.")
@arg("--base-delay", help="Base delay to use (in seconds).")
@arg("--batch-delay",
     help="Delay to use between batches of jobs (in seconds).")
@arg("--delay-every-n-jobs",
     help="Number of jobs after which a delay is imposed.")
@arg("--dry-run", help="Dry run: do not submit or execute jobs.")
def run(settings, logs_outdir,
        use_cluster=True,
        base_delay=10,
        # Batch delay (20 mins by default)
        batch_delay=60*20,
        delay_every_n_jobs=30,
        dry_run=False,
        event_types=None):
    """
    Run MISO on a set of samples.
    """
    settings_filename = utils.pathify(settings)
    if event_types is not None:
        print "Only running MISO on event types: ", event_types
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="run")
    output_dir = misowrap_obj.miso_outdir
    bam_files = misowrap_obj.bam_files
    read_len = misowrap_obj.read_len
    overhang_len = misowrap_obj.overhang_len
    events_dir = misowrap_obj.miso_events_dir
    single_end = True
    if misowrap_obj.insert_lens_dir is not None:
        insert_lens_dir = misowrap_obj.insert_lens_dir
        misowrap_obj.logger.info("Running in paired-end mode...")
        misowrap_obj.logger.info(" - Insert length directory: %s" \
                                 %(insert_lens_dir))
        single_end = False
    else:
        misowrap_obj.logger.info("Running in single-end mode...")        
    run_events_analysis = misowrap_obj.run_events_cmd
    event_types_dirs = \
        miso_utils.get_event_types_dirs(misowrap_obj.settings_info)
    miso_settings_filename = misowrap_obj.miso_settings_filename
    n = 0
    for bam_input in bam_files:
        bam_filename, sample_label = bam_input
        bam_filename = utils.pathify(bam_filename)
        misowrap_obj.logger.info("Processing: %s" %(bam_filename))
        for event_type_dir in event_types_dirs:
            event_type = os.path.basename(event_type_dir)
            if event_types is not None:
                if event_type not in event_types:
                    print "Skipping event type: %s" %(event_type)
                    continue
            print "  - Using event dir: %s" %(event_type_dir)
            miso_cmd = "%s" %(run_events_analysis)
            bam_basename = os.path.basename(bam_filename)
            # Output directory for sample
            sample_output_dir = os.path.join(output_dir, 
                                             sample_label,
                                             event_type)
            # Pass sample to MISO along with event
            miso_cmd += " --run %s %s" %(event_type_dir,
                                         bam_filename)
            if not single_end:
                insert_len_filename = \
                    os.path.join(insert_lens_dir,
                                 "%s.insert_len" %(bam_basename))
                misowrap_obj.logger.info("Reading paired-end parameters " \
                                         "from file...")
                misowrap_obj.logger.info("  - PE file: %s" \
                                         %(insert_len_filename))
                pe_params = miso_utils.read_pe_params(insert_len_filename)
                # Paired-end parameters
                miso_cmd += " --paired-end %.2f %.2f" %(pe_params["mean"],
                                                        pe_params["sdev"])
            # Read length
            miso_cmd += " --read-len %d" %(read_len)
            # Overhang length
            miso_cmd += " --overhang-len %d" %(overhang_len)
            # Prefilter?
            if misowrap_obj.prefilter_miso:
                miso_cmd += " --prefilter"
            # Output directory
            miso_cmd += " --output-dir %s" %(sample_output_dir)
            # Use cluster
            if misowrap_obj.use_cluster:
                miso_cmd += " --use-cluster"
                miso_cmd += " --chunk-jobs %d" %(misowrap_obj.chunk_jobs)
            # Settings
            miso_cmd += " --settings %s" %(miso_settings_filename)
            misowrap_obj.logger.info("Executing: %s" %(miso_cmd))
            job_name = "%s_%s" %(sample_label, event_type)
            if use_cluster:
                if not dry_run:
                    misowrap_obj.my_cluster.launch_job(miso_cmd,
                                                       job_name,
                                                       ppn=1)
                if n == delay_every_n_jobs:
                    # Larger delay everytime we've submitted n jobs
                    misowrap_obj.logger.info("Submitted %d jobs, now waiting %.2f mins." \
                                             %(n, batch_delay / 60.))
                    time.sleep(batch_delay)
                    n = 0
                time.sleep(base_delay)
            else:
                if not dry_run:
                    os.system(miso_cmd)
            n += 1


@arg("settings", help="misowrap settings filename.")
@arg("logs_outdir", help="directory where to place logs.")
@arg("--dry-run", help="Dry run. Do not execute any jobs or commands.")
def filter(settings,
           logs_outdir,
           dry_run=False):
    """
    Output a set of filtered MISO comparisons.
    """
    settings_filename = utils.pathify(settings)
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="filter")
    misowrap_obj.logger.info("Filtering MISO events...")
    psi_table = pt.PsiTable(misowrap_obj)
    psi_table.output_filtered_comparisons()


@arg("settings", help="misowrap settings filename.")
@arg("logs_outdir", help="directory where to place logs.")
@arg("--dry-run", help="Dry run. Do not execute any jobs or commands.")
def compute_insert_lens(settings,
                        output_dir,
                        dry_run=False):
    """
    Compute insert lengths for all samples.
    """
    settings_filename = utils.pathify(settings)
    logs_outdir = utils.pathify(logs_outdir)
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="insert_lens")
    const_exons_gff = misowrap_obj.const_exons_gff
    if not os.path.isfile(const_exons_gff):
        print "Error: %s const exons GFF does not exist." \
            %(const_exons_gff)
        sys.exit(1)
    pe_utils_path = misowrap_obj.pe_utils_cmd 
    insert_len_output_dir = os.path.join(output_dir, "insert_lens")
    num_bams = len(misowrap_obj.bam_files)
    
    print "Computing insert lengths for %d files" %(num_bams)
    for bam_filename, sample_name in misowrap_obj.bam_files:
        print "Processing: %s" %(bam_filename)
        insert_len_cmd = "%s --compute-insert-len %s %s --output-dir %s" \
            %(pe_utils_path,
              bam_filename,
              const_exons_gff,
              insert_len_output_dir)
        print "Executing: %s" %(insert_len_cmd)
        job_name = "%s_insert_len" %(sample_name)
        if misowrap_obj.use_cluster:
            misowrap_obj.my_cluster.launch_job(insert_len_cmd,
                                               job_name,
                                               ppn=1)
        else:
            os.system(insert_len_cmd)


@arg("settings", help="misowrap settings filename.")
@arg("logs-outdir", help="Directory where to place logs.")
@arg("--delay", help="Delay between execution of cluster jobs")
@arg("--dry-run", help="Dry run: do not submit or execute jobs.")
def combine_comparisons(settings,
                        logs_outdir,
                        common_cols=["isoforms",
                                     "chrom",
                                     "strand",
                                     "mRNA_starts",
                                     "mRNA_ends",
                                     "gene_id",
                                     "gene_symbol"],
                        delay=5,
                        dry_run=False,
                        NA_VAL="NA"):
    """
    Output combined MISO comparisons. For each event type,
    combine the MISO comparisons for the relevant groups
    based on the 'comparison_groups' in the misowrap
    settings file.
    """
    settings_filename = utils.pathify(settings)
    logs_outdir = utils.pathify(logs_outdir)
    utils.make_dir(logs_outdir)
    misowrap_obj = mw.MISOWrap(settings_filename,
                               logs_outdir,
                               logger_label="combine_comparisons")
    comparisons_dir = misowrap_obj.comparisons_dir    
    if not os.path.isdir(comparisons_dir):
        misowrap_obj.logger.critical("Comparisons directory %s not found. " \
                                     %(comparisons_dir))
        sys.exit(1)
    # Comparison types to combine: unfiltered comparisons and filtered comparisons
    # (if available)
    unfiltered_comp_dir = os.path.join(comparisons_dir,
                                       "combined_comparisons")
    filtered_comp_dir = os.path.join(comparisons_dir,
                                     "filtered_events")
    dirs_to_process = [unfiltered_comp_dir, filtered_comp_dir]
    comparison_groups = misowrap_obj.comparison_groups
    for curr_comp_dir in dirs_to_process:
        if not os.path.isdir(curr_comp_dir):
            print "Comparisons directory %s not found, skipping" %(curr_comp_dir)
            continue
        # For each event type, output the sample comparisons
        for event_type in misowrap_obj.event_types:
            # Collection of MISO comparison dataframes (to be merged later)
            # for the current event type
            comparison_dfs = []
            comparison_labels = []
            event_dir = os.path.join(curr_comp_dir, event_type)
            if not os.path.isdir(event_dir):
                misowrap_obj.logger.info("Cannot find event type %s dir, " \
                                         "skipping..." %(event_type))
                continue
            # Look only at sample comparisons within each sample group
            for comp_group in comparison_groups:
                sample_pairs = utils.get_pairwise_comparisons(comp_group)
                misowrap_obj.logger.info("  - Total of %d comparisons" \
                                         %(len(sample_pairs)))
                for sample1, sample2 in sample_pairs:
                    # Load miso_bf file for the current comparison
                    # and join it to the combined df
                    comparison_name = "%s_vs_%s" %(sample1, sample2)
                    bf_data = miso_utils.load_miso_bf_file(event_dir,
                                                           comparison_name,
                                                           substitute_labels=True)
                    if bf_data is None:
                        misowrap_obj.logger.warning("Could not find comparison %s" \
                                                    %(comparison_name))
                        continue
                    comparison_dfs.append(bf_data)
                    comparison_labels.append(comparison_name)
            # Merge the comparison dfs together
            print "Merging comparisons for %s" %(event_type)
            combined_df = pandas_utils.combine_dfs(comparison_dfs)
            output_filename = os.path.join(output_dir,
                                           "%s.miso_bf" %(event_type))
            misowrap_obj.logger.info("Outputting %s results to: %s" \
                                     %(event_type, output_filename))
            if not dry_run:
                combined_df.to_csv(output_filename,
                                   float_format="%.4f",
                                   sep="\t",
                                   na_rep=NA_VAL,
                                   index=True)


def greeting(parser=None):
    print "misowrap: wrapper for running MISO and parsing its results.\n"
    if parser is not None:
        parser.print_help()
            
        
def main():
    greeting()
    
    argh.dispatch_commands([
        run,
        summarize,
        compare,
        combine_comparisons,
    ])
    # from optparse import OptionParser
    # parser = OptionParser()
      
    # parser.add_option("--run", dest="run", nargs=1, default=None,
    #                   help="Run MISO on a set of events. "
    #                   "Takes a settings filename.")
    # parser.add_option("--summarize", dest="summarize", nargs=1, default=None,
    #                   help="Run MISO summarize on a set of samples. "
    #                   "Takes a settings filename.")
    # parser.add_option("--compare", dest="compare", nargs=1, default=None,
    #                   help="Run MISO sample comparisons on all pairwise "
    #                   "comparisons. Takes a settings filename.")
    # parser.add_option("--filter", dest="filter", nargs=1,
    #                   default=None,
    #                   help="Filter a set of MISO comparisons. "
    #                   "Takes a settings filename.")
    # parser.add_option("--compute-insert-lens", dest="compute_insert_lens",
    #                   nargs=1, default=None,
    #                   help="Compute insert lengths for a set of BAM files. " 
    #                   "Takes a settings filename.")
    # parser.add_option("--combine-comparisons", dest="combine_comparisons",
    #                   nargs=2, default=None,
    #                   help="Combine MISO comparisons into a single file, "
    #                   "based on the \'comparison_groups\' parameter in the "
    #                   "settings file. Takes a MISO comparisons directory and "
    #                   "a settings filename.")
    # parser.add_option("--output-dir", dest="output_dir", default=None,
    #                   help="Output directory.")
    # parser.add_option("--event-types", dest="event_types", default=None, type="str",
    #                   help="Optional: comma-separated list of event types to run on "
    #                   "only, e.g. A3SS,SE")
    # (options, args) = parser.parse_args()


    # if options.run != None:
    #     settings_filename = utils.pathify(options.run)
    #     event_types = None
    #     if options.event_types is not None:
    #         event_types = options.event_types.split(",")
    #     run_miso_on_samples(settings_filename, output_dir,
    #                         event_types=event_types)

    # if options.summarize != None:
    #     settings_filename = utils.pathify(options.summarize)
    #     summarize_miso_samples(settings_filename, output_dir)
        
    # if options.compare != None:
    #     settings_filename = utils.pathify(options.compare)
    #     compare_miso_samples(settings_filename, output_dir)

    # if options.filter != None:
    #     settings_filename = utils.pathify(options.filter)
    #     filter_events(settings_filename, output_dir)

    # if options.compute_insert_lens != None:
    #     settings_filename = utils.pathify(options.compute_insert_lens)
    #     compute_insert_lens(settings_filename, output_dir)

    # if options.combine_comparisons != None:
    #     comparisons_dir = options.combine_comparisons[0]
    #     settings_filename = utils.pathify(options.combine_comparisons[1])
    #     combine_comparisons(settings_filename,
    #                         comparisons_dir,
    #                         output_dir)
        

if __name__ == '__main__':
    main()

