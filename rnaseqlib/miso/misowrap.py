#!/usr/bin/env python
import os
import sys
import time
import glob
import pandas
import itertools
from collections import defaultdict

import yklib
import yklib.settings as settings
import yklib.cluster as cluster

def read_pe_params(insert_len_filename):
    """
    Get paired-end parameters from .insert_len file.
    """
    insert_len_filename = os.path.abspath(os.path.expanduser(insert_len_filename))
    if not os.path.isfile(insert_len_filename):
        print "Error: %s not a file." %(insert_len_filename)
        sys.exit(1)

    insert_file = open(insert_len_filename, "r")
    fields = insert_file.readline()[1:].strip().split(",")
    pe_params = {}
    for field in fields:
        k, v = field.split("=")
        pe_params[k] = float(v)
    insert_file.close()
    return pe_params


def load_miso_bf_file(comparisons_dir, comparison_name):
    """
    Load MISO information for a comparison name.
    """
    sample_comparison_dir = os.path.join(comparisons_dir, comparison_name)
    bf_filename = get_bf_filename(sample_comparison_dir)
    if bf_filename is None or (not os.path.isfile(bf_filename)):
        return None
    miso_bf_data = pandas.read_table(bf_filename, sep="\t")
    return miso_bf_data
    

def get_event_types_dirs(settings_info):
    """
    Return event types.
    """
    miso_events_dir = os.path.abspath(os.path.expanduser(settings_info["settings"]["miso_events_dir"]))
    event_types_dirs = [os.path.join(miso_events_dir, dirname) \
                        for dirname in os.listdir(miso_events_dir)]
    return event_types_dirs


def summarize_miso_samples(settings_filename,
                           miso_output_dir):
    """
    Summarize samples in MISO directory.
    """
    settings_info, parsed_settings = \
        settings.load_settings(settings_filename)
    miso_samples_dirs = os.listdir(miso_output_dir)
    miso_dir = settings_info["settings"]["miso_dir"]
    sample_labels = settings_info["data"]["sample_labels"]
    print "Summarizing MISO output..."
    summarize_cmd = os.path.join(miso_dir, "run_miso.py")

    for sample_label in sample_labels:
        sample_basename = sample_label[0]
        sample_dir_path = os.path.abspath(os.path.join(miso_output_dir,
                                                       sample_basename))
        print "Processing: %s" %(sample_basename)
        if not os.path.isdir(sample_dir_path):
            print "Skipping non-directory: %s" %(sample_dir_path)
        # List all event directories in the sample
        event_dirs = os.listdir(sample_dir_path)
        for event_dirname in event_dirs:
            event_dir_path = os.path.abspath(os.path.join(sample_dir_path, event_dirname))
            if not os.path.isdir(event_dir_path):
                print "Skipping non-dir: %s" %(event_dir_path)
            print "Processing event type: %s" %(event_dirname)
            summary_cmd = "%s --summarize-samples %s %s --summary-label %s" \
                %(summarize_cmd,
                  event_dir_path,
                  event_dir_path,
                  sample_basename)
            job_name = "summarize_%s_%s" %(sample_basename,
                                           os.path.basename(event_dirname))
            print "Executing: %s" %(summary_cmd)
            print "event_dir_path: %s" %(event_dir_path)
            cluster.run_on_cluster(summary_cmd, job_name,
                                   event_dir_path)


def get_summary_filename(sample_dir):
    """
    Get summary filename from directory.
    """
    summary_dir = os.path.join(sample_dir, "summary")
    if not os.path.isdir(summary_dir):
        raise Exception, "%s not a summary dir." %(summary_dir)
    summary_files = glob.glob(os.path.join(summary_dir,
                                           "*.miso_summary"))
    summary_files = [os.path.join(summary_dir, fname) \
                     for fname in summary_files]
    if len(summary_files) > 1:
        raise Exception, "Warning: more than 1 summary file for %s" \
            %(summary_dir)
    return summary_files[0]

    
def get_bf_filename(pairwise_comparison_dir):
    """
    Return a Bayes factor filename from a
    pairwise comparisons directory.
    """
    pairwise_comparison_dir = os.path.abspath(os.path.expanduser(pairwise_comparison_dir))
    bf_dir = os.path.join(pairwise_comparison_dir, "bayes-factors")
    if not os.path.isdir(bf_dir):
        print "WARNING: Could not get BF dir %s" %(bf_dir)
        return None
    bf_filename = glob.glob(os.path.join(bf_dir,
                                         "*.miso_bf"))
    if len(bf_filename) > 1:
        print "Error: Multiple BF filenames in %s" %(bf_dir)
        return None
    bf_filename = bf_filename[0]
    return bf_filename
    

def get_comparisons_dirs(comparisons_dir):
    """
    Get all comparisons directories.
    """
    comparisons_dirs = []
    comparisons_dir = os.path.abspath(os.path.expanduser(comparisons_dir))
    candidate_dirs = glob.glob(os.path.join(comparisons_dir, "*_vs_*"))
    for dirname in candidate_dirs:
        if not os.path.isdir(dirname):
            continue
        comparisons_dirs.append(dirname)
    return comparisons_dirs


def get_pairwise_from_sets(first_samples, second_samples):
    seen_pairs = []
    for sample_pair in itertools.product(first_samples,
                                         second_samples):
        sample1, sample2 = sample_pair
        if (sample_pair in seen_pairs) or (sample1 == sample2) or \
            ((sample2, sample1) in seen_pairs):
            continue
        seen_pairs.append(sample_pair)
    return seen_pairs
    
        
def get_pairwise_comparisons(samples):
    """
    Return pairwise comparisons between samples.
    """
    return get_pairwise_from_sets(samples, samples)
    

def compare_miso_samples(settings_filename,
                         miso_output_dir,
                         event_types=["A3SS",
                                      "A5SS",
                                      "AFE",
                                      "ALE",
                                      "MXE",
                                      "RI",
                                      "SE",
                                      "SE_noAceView",
                                      "SE_shortest_noAceView",
                                      "MXE_shortest_noAceView",
                                      "A3SS_shortest_noAceView",
                                      "A5SS_shortest_noAceView",
                                      "TandemUTR",
                                      "TandemUTR_3pseq"]):
    """
    Run a MISO samples comparison between all pairs of samples.
    """
    settings_info, parsed_settings = \
        settings.load_settings(settings_filename)
    miso_samples_dirs = os.listdir(miso_output_dir)
    miso_dir = settings_info["settings"]["miso_dir"]
    sample_labels = settings_info["data"]["sample_labels"]
    print "Running MISO comparisons..."
    run_miso_cmd = os.path.join(miso_dir, "run_miso.py")
    comparisons_dir = os.path.join(miso_output_dir,
                                   "comparisons")
    if not os.path.isdir(comparisons_dir):
        os.makedirs(comparisons_dir)

    sample_pairs = get_pairwise_comparisons(sample_labels)
    print "Running total of %d comparisons" %(len(sample_pairs))

    for sample1, sample2 in sample_pairs:
        # For each pair of samples, compare their output
        # along each event type
        print "Comparing %s %s" %(sample1, sample2)
        sample1_name = sample1[0]
        sample2_name = sample2[0]
        # Directories for each sample
        sample1_dir = os.path.join(miso_output_dir, sample1_name)
        sample2_dir = os.path.join(miso_output_dir, sample2_name)
        for event_type in event_types:
            print "Processing %s..." %(event_type)
            sample1_event_dir = os.path.join(sample1_dir, event_type)
            sample2_event_dir = os.path.join(sample2_dir, event_type)
            job_name = "compare_%s_%s_%s" %(sample1_name,
                                            sample2_name,
                                            event_type)
            event_comparisons_dir = os.path.join(comparisons_dir, event_type)
            compare_cmd = "%s --compare-samples %s %s %s " \
                "--comparison-labels %s %s" %(run_miso_cmd,
                                              sample1_event_dir,
                                              sample2_event_dir,
                                              event_comparisons_dir,
                                              sample1_name,
                                              sample2_name)
            print "Executing: %s" %(compare_cmd)
            cluster.run_on_cluster(compare_cmd, job_name,
                                   event_comparisons_dir)
            
        # for event_dirname in event_dirs:
        #     event_dir_path = os.path.abspath(os.path.join(sample_dir_path, event_dirname))
        #     if not os.path.isdir(event_dir_path):
        #         print "Skipping non-dir: %s" %(event_dir_path)
        #     print "Processing event type: %s" %(event_dirname)
        #     compare_cmd = "%s --summarize-samples %s %s" %(summarize_cmd,
        #                                                    event_dir_path,
        #                                                    event_dir_path)
        #     job_name = "summarize_%s_%s" %(sample_basename,
        #                                    os.path.basename(event_dirname))
        #     print "Executing: %s" %(summary_cmd)
        #     print "event_dir_path: %s" %(event_dir_path)
        #     cluster.run_on_cluster(summary_cmd, job_name,
        #                            event_dir_path)
            

def run_miso_on_samples(settings_filename, output_dir,
                        chunk_jobs=600,
                        use_cluster=True):
    """
    Run MISO on a set of samples.
    """
    ##
    ## TODO: make it so it uses the miso_output_dir variable from
    ## settings file instead of requiring parameter
    ##
    settings_info, parsed_settings = \
                   settings.load_settings(settings_filename)
    bam_files = settings_info["data"]["bam_files"]
    read_len = int(settings_info["settings"]["read_len"])
#    overhang_len = int(settings_info["settings"]["overhang_len"])    
    miso_dir = settings_info["settings"]["miso_dir"]
    events_dir = settings_info["settings"]["miso_events_dir"]
    single_end = False
    if "insert_lens_dir" not in settings_info["data"]:
        single_end = True
    else:
        insert_lens_dir = settings_info["data"]["insert_lens_dir"]
    run_events_analysis = os.path.join(miso_dir,
                                       "run_events_analysis.py")
    event_types_dirs = get_event_types_dirs(settings_info)
    miso_settings_filename = settings_info["settings"]["miso_settings_filename"]
    
    for bam_filename in bam_files:
        print "Processing: %s" %(bam_filename)
        for event_type_dir in event_types_dirs:
            event_type = os.path.basename(event_type_dir)
            print "  - Using event dir: %s" %(event_type_dir)
            miso_cmd = "%s" %(run_events_analysis)
            bam_basename = os.path.basename(bam_filename)
            sample_label = bam_basename.split(".")[0]
            # Output directory for sample
            sample_output_dir = os.path.join(output_dir, sample_label,
                                             event_type)
           # Pass sample to MISO along with event
            miso_cmd += " --compute-genes-psi %s %s" %(event_type_dir,
                                                       bam_filename)
            if not single_end:
                insert_len_filename = os.path.join(insert_lens_dir,
                                                   "%s.insert_len" %(bam_basename))
                print "Reading paired-end parameters from file..."
                pe_params = read_pe_params(insert_len_filename)
                # Paired-end parameters
                miso_cmd += " --paired-end %.2f %.2f" %(pe_params["mean"],
                                                        pe_params["sdev"])
            # Read length
            miso_cmd += " --read-len %d" %(read_len)
            # Output directory
            miso_cmd += " --output-dir %s" %(sample_output_dir)
            # Use cluster
            miso_cmd += " --use-cluster"
            miso_cmd += " --chunk-jobs %d" %(chunk_jobs)
            # Settings
            miso_cmd += " --settings %s" %(miso_settings_filename)
            print "Executing: %s" %(miso_cmd)
            job_name = "%s_%s" %(sample_label, event_type)
            if use_cluster:
                cluster.run_on_cluster(miso_cmd, job_name,
                                       output_dir)
            else:
                os.system(miso_cmd)

                
def compute_insert_lengths(settings_filename, output_dir):
    settings_info, parsed_settings = \
                   settings.load_settings(settings_filename)
    miso_dir = settings_info["settings"]["miso_dir"]
    bam_files = settings_info["data"]["bam_files"]
    const_exons_gff = settings_info["data"]["const_exons_gff"]
    const_exons_gff = os.path.abspath(os.path.expanduser(const_exons_gff))

    if not os.path.isfile(const_exons_gff):
        print "Error: %s const exons GFF does not exist." \
            %(const_exons_gff)
        sys.exit(1)

    miso_dir = os.path.abspath(os.path.expanduser(miso_dir))
    pe_utils_path = os.path.join(miso_dir, "pe_utils.py")
    insert_len_output_dir = os.path.join(output_dir, "insert_lens")
    num_bams = len(bam_files)
    
    print "Computing insert lengths for %d files" %(num_bams)
    for bam_filename in bam_files:
        print "Processing: %s" %(bam_filename)
        insert_len_cmd = "%s --compute-insert-len %s %s --output-dir %s" \
            %(pe_utils_path,
              bam_filename,
              const_exons_gff,
              insert_len_output_dir)
        print "Executing: %s" %(insert_len_cmd)
        sample_name = os.path.basename(bam_filename)
        job_name = sample_name.split(".bam")[0]
        cluster.run_on_cluster(insert_len_cmd, job_name,
                               insert_len_output_dir,
                               cmd_name="bsub")

        
def main():
    from optparse import OptionParser
    parser = OptionParser()
      
    parser.add_option("--run", dest="run", nargs=1, default=None,
                      help="Run MISO on a set of events. Takes a settings filename.")
    parser.add_option("--summarize", dest="summarize", nargs=1, default=None,
                      help="Run MISO summarize on a set of samples. Takes a settings filename.")
    parser.add_option("--compare", dest="compare", nargs=1, default=None,
                      help="Run MISO sample comparisons on all pairwise comparisons. "
                      "Takes a settings filename.")
    parser.add_option("--compute-insert-lens", dest="compute_insert_lens", nargs=1,
                      default=None,
                      help="Compute insert lengths for a set of BAM files. " \
                      "takes a settings filename.")
    parser.add_option("--output-dir", dest="output_dir", default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.output_dir == None:
        print "Error: need --output-dir."
        sys.exit(1)
        
    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.run != None:
        settings_filename = os.path.abspath(os.path.expanduser(options.run))
        run_miso_on_samples(settings_filename, output_dir)

    if options.summarize != None:
        settings_filename = os.path.abspath(os.path.expanduser(options.summarize))
        summarize_miso_samples(settings_filename, output_dir)
        
    if options.compare != None:
        settings_filename = os.path.abspath(os.path.expanduser(options.compare))
        compare_miso_samples(settings_filename, output_dir)

    if options.compute_insert_lens != None:
        settings_filename = os.path.abspath(os.path.expanduser(options.compute_insert_lens))
        compute_insert_lengths(settings_filename, output_dir)
        

if __name__ == '__main__':
    main()

