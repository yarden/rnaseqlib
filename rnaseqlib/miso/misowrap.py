##
## misowrap: a wrapper to running MISO on a set of samples
## and processing its output
##
import os
import sys
import time
import glob
import pandas
import itertools
from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.miso.miso_utils as miso_utils

class MISOWrap:
    """
    Object containing information about a set of samples to be
    processed by MISO and their MISO output.
    """
    def __init__(self, settings_filename):
        self.settings_filename = settings_filename
        # Load settings
        # ...


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

    sample_pairs = utils.get_pairwise_comparisons(sample_labels)
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
    event_types_dirs = miso_utils.get_event_types_dirs(settings_info)
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
                pe_params = miso_utils.read_pe_params(insert_len_filename)
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

