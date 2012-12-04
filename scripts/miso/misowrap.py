##
## misowrap: a wrapper to run MISO on a set of samples
## and processing its output.
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
import rnaseqlib.miso.misowrap_settings as misowrap_settings
import rnaseqlib.miso.miso_utils as miso_utils
import rnaseqlib.cluster_utils.cluster as cluster


##
## Alternative event types to consider
##
EVENT_TYPES = ["A3SS",
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
               "TandemUTR_3pseq"]

class MISOWrap:
    """
    Object containing information about a set of samples to be
    processed by MISO and their MISO output.
    """
    def __init__(self, settings_filename, output_dir):
        self.settings_filename = settings_filename
        self.settings_info = None
        # Main output directory
        self.output_dir = utils.pathify(output_dir)
        utils.make_dir(self.output_dir)
        # MISO output directory (where raw output is)
        self.miso_outdir = None
        # Comparisons output directory
        self.comparisons_outdir = None
        # BAM files to process
        self.bam_files = None
        # Sample labels
        self.sample_labels = None
        # Insert length directory (for paired-end samples)
        self.insert_lens_dir = None
        # Logs output directory
        self.logs_outdir = None
        # Logger object
        self.logger = None
        # Cluster submission object
        self.my_cluster = None
        # Event types to process
        self.event_types = None
        # Whether to submit jobs to cluster
        self.use_cluster = False
        # Load settings
        self.load_settings()

        
    def load_settings(self):
        """
        Load settings for misowrap.
        """
        settings_info, parsed_settings = \
                  misowrap_settings.load_misowrap_settings(self.settings_filename)
        self.settings_info = settings_info
        # Load basic settings about data
        self.read_len = self.settings_info["settings"]["readlen"]
        self.overhang_len = self.settings_info["settings"]["overhanglen"]
        self.miso_bin_dir = \
          utils.pathify(self.settings_info["settings"]["miso_bin_dir"])
        self.miso_settings_filename = \
          utils.pathify(self.settings_info["settings"]["miso_settings_filename"])
        self.miso_events_dir = \
          utils.pathify(self.settings_info["settings"]["miso_events_dir"])
        self.miso_outdir = \
          utils.pathify(self.settings_info["settings"]["miso_output_dir"])
        # Load data-related parameters
        self.bam_files = self.settings_info["data"]["bam_files"]
        if "insert_lens_dir" in self.settings_info["data"]:
            self.insert_lens_dir = \
              utils.pathify(self.settings_info["data"]["insert_lens_dir"])
        # Set output directories
        self.comparisons_dir = os.path.join(self.output_dir, 
                                            "comparisons")
        self.logs_outdir = os.path.join(self.output_dir,
                                        "misowrap_logs")
        # Create necessary directories
        utils.make_dir(self.miso_outdir)
        utils.make_dir(self.comparisons_dir)
        utils.make_dir(self.logs_outdir)
        if "cluster_type" in self.settings_info["settings"]:
            self.use_cluster = True
            self.cluster_type = self.settings_info["settings"]["cluster_type"]
            self.chunk_jobs = self.settings_info["settings"]["chunk_jobs"]
        if self.use_cluster:
            # Load cluster object if given a cluster type
            self.load_cluster()
        # Create a logger object 
        self.logger = utils.get_logger("misowrap",
                                       self.logs_outdir)
        # Load event types
        self.load_event_types()


    def load_event_types(self):
        """
        Read a list of event types that were defined in the 
        events indexed directory.
        """
        if not os.path.isdir(self.miso_events_dir):
            print "Error: %s is not a directory. Need a directory " \
                  "with indexed events." %(self.miso_events_dir)
            sys.exit(1)
        self.event_types = []
        for fname in os.listdir(self.miso_events_dir):
            if os.path.isdir(fname):
                event_name = os.path.basename(dirname)
                self.event_types.append(event_name)

        
    def load_cluster(self):
        """
        Load cluster submission object for the particular
        pipeline settings we were given.
        """
        self.my_cluster = \
          cluster.Cluster(self.cluster_type,
                          self.output_dir,
                          self.logger)
        

    def summarize_miso_samples(self):
        """
        Summarize samples in MISO directory.
        """
        miso_samples_dirs = os.listdir(miso_output_dir)
        miso_bin_dir = self.settings_info["settings"]["miso_bin_dir"]
        sample_labels = self.settings_info["data"]["sample_labels"]
        print "Summarizing MISO output..."
        summarize_cmd = os.path.join(miso_bin_dir, "run_miso.py")

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
                #                cluster.run_on_cluster(summary_cmd, job_name,
                #                                       event_dir_path)

    def __repr__(self):
        repr_str = "MISOWrap(settings=%s, output_dir=%s)" \
          %(self.settings_filename, 
            self.output_dir)
        return repr_str


    def __str__(self):
        return self.__repr__()

            
def compare_miso_samples(settings_filename,
                         miso_output_dir,
                         event_types=EVENT_TYPES):
    """
    Run a MISO samples comparison between all pairs of samples.
    """
    settings_info, parsed_settings = \
        settings.load_settings(settings_filename)
    miso_samples_dirs = os.listdir(miso_output_dir)
    miso_dir = settings_info["settings"]["miso_bin_dir"]
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
            #            cluster.run_on_cluster(compare_cmd, job_name,
            #                                   event_comparisons_dir)
            
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
                        use_cluster=True):
    """
    Run MISO on a set of samples.
    """
    misowrap_obj = MISOWrap(settings_filename, output_dir)
    bam_files = misowrap_obj.bam_files
    read_len = misowrap_obj.read_len
    overhang_len = misowrap_obj.overhang_len
    miso_bin_dir = misowrap_obj.miso_bin_dir
    events_dir = misowrap_obj.miso_events_dir
    single_end = False
    if misowrap_obj.insert_lens_dir is not None:
        print "Running in single-end mode..."
        single_end = True
    else:
        insert_lens_dir = misowrap_obj.insert_lens_dir
        print "Running in paired end mode..."
        print "  - Insert length directory: %s" %(insert_lens_dir)
    run_events_analysis = os.path.join(miso_bin_dir,
                                       "run_events_analysis.py")
    event_types_dirs = miso_utils.get_event_types_dirs(misowrap_obj.settings_info)
    miso_settings_filename = misowrap_obj.miso_settings_filename
    for bam_input in bam_files:
        sample_label, bam_filename = bam_input
        bam_filename = utils.pathify(bam_filename)
        print "Processing: %s" %(bam_filename)
        for event_type_dir in event_types_dirs:
            event_type = os.path.basename(event_type_dir)
            print "  - Using event dir: %s" %(event_type_dir)
            miso_cmd = "%s" %(run_events_analysis)
            bam_basename = os.path.basename(bam_filename)
            # Output directory for sample
            sample_output_dir = os.path.join(output_dir, 
                                             sample_label,
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
            if misowrap_obj.use_cluster:
                miso_cmd += " --use-cluster"
                miso_cmd += " --chunk-jobs %d" %(misowrap_obj.chunk_jobs)
            # Settings
            miso_cmd += " --settings %s" %(miso_settings_filename)
            print "Executing: %s" %(miso_cmd)
            job_name = "%s_%s" %(sample_label, event_type)
            if use_cluster:
                pass
                #                cluster.run_on_cluster(miso_cmd, job_name,
                #                                       output_dir)
            else:
                pass
                #os.system(miso_cmd)

                
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
        #        cluster.run_on_cluster(insert_len_cmd, job_name,
        #                               insert_len_output_dir,
        #                               cmd_name="bsub")

        
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

