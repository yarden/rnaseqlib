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
import rnaseqlib.miso.miso_utils as miso_utils
import rnaseqlib.cluster_utils.cluster as cluster


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
        self.comparison_groups = None
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
        # run_miso cmd
        self.run_miso_cmd = None
        # run_events_analysis cmd
        self.run_events_cmd = None
        # Constitutive exons GFF file: used to compute
        # the insert length distribution
        self.const_exons_gff = None
        # Load settings
        self.load_settings()
        ##
        ## Load annotation of events, like a map
        ## events to genes.
        ##
        self.events_to_genes = None
        self.load_events_to_genes()


    def load_events_to_genes(self,
                             source="ensGene",
                             delimiter="\t"):
        """
        Load mapping from events to genes.

        Expects a directory with files named
        according to events, e.g.:
        
          SE.mm9.gff3_to_ensGene.txt        
        """
        if "events_to_genes_dir" not in self.settings_info["settings"]:
            return
        events_to_genes_dir = \
            self.settings_info["settings"]["events_to_genes_dir"]
        events_to_genes_dir = utils.pathify(events_to_genes_dir)
        print "Loading events to genes mapping from: %s" \
            %(events_to_genes_dir)
        # If we're given mapping from events to genes, load
        # these and index them by event type.
        if not os.path.isdir(events_to_genes_dir):
            print "Error: %s not a directory."
            sys.exit(1)
        basename_card = "*_to_%s.txt" %(source)
        events_to_genes_files = \
            glob.glob(os.path.join(events_to_genes_dir,
                                   basename_card))
        if len(events_to_genes_files) == 0:
            print "Error: %s directory contains no %s files." \
                %(events_to_genes_dir,
                  basename_card)
            sys.exit(1)
        self.events_to_genes = defaultdict(lambda: defaultdict(list))
        for fname in events_to_genes_files:
            # Extract event type based on filename
            event_type = os.path.basename(fname).split(".")[0]
            with open(fname, "r") as events_file:
                events_entries = csv.DictReader(events_file,
                                                delimiter=delimiter)
                for entry in events_entries:
                    event_id = entry["event_id"]
                    # Parse genes into a list
                    genes = entry["gene_id"].split(",")
                    # Index events by their type and then by
                    # their ID
                    self.events_to_genes[event_type][event_id].extend(genes)

                
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
        # Sample labels
        self.sample_labels = self.settings_info["data"]["sample_labels"]
        # Set output directories
        self.comparisons_dir = os.path.join(self.output_dir, 
                                            "comparisons")
        self.comparison_groups = self.settings_info["data"]["comparison_groups"]
        self.logs_outdir = os.path.join(self.output_dir,
                                        "misowrap_logs")
        # Create necessary directories
        utils.make_dir(self.miso_outdir)
        utils.make_dir(self.comparisons_dir)
        utils.make_dir(self.logs_outdir)
        if "cluster_type" in self.settings_info["settings"]:
            self.use_cluster = True
            self.cluster_type = \
                self.settings_info["settings"]["cluster_type"]
            self.chunk_jobs = \
                self.settings_info["settings"]["chunk_jobs"]
        if self.use_cluster:
            print "Loading cluster information."
            # Load cluster object if given a cluster type
            self.load_cluster()
        # Create a logger object 
        self.logger = utils.get_logger("misowrap",
                                       self.logs_outdir)
        # Load event types
        self.load_event_types()
        # Set path to MISO scripts
        self.run_miso_cmd = os.path.join(self.miso_bin_dir,
                                         "run_miso.py")
        self.run_events_cmd = os.path.join(self.miso_bin_dir,
                                           "run_events_analysis.py")
        self.pe_utils_cmd = os.path.join(self.miso_bin_dir,
                                         "pe_utils.py")
        # Files related to gene tables
        self.tables_dir = \
            os.path.join(self.settings_info["pipeline-files"]["init_dir"],
                         "ucsc")
        if not os.path.isdir(self.tables_dir):
            print "Error: %s directory does not exist." \
                %(self.tables_dir)
            sys.exit(1)
        self.const_exons_gff = os.path.join(self.tables_dir,
                                            "exons",
                                            "const_exons",
                                            "ensGene.const_exons.gff")
        if not os.path.isfile(self.const_exons_gff):
            print "Error: Const. exons GFF %s does not exist." \
                %(self.const_exons_gff)
            sys.exit(1)


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
            fname = os.path.join(self.miso_events_dir, fname)
            if os.path.isdir(fname):
                event_name = os.path.basename(fname)
                self.event_types.append(event_name)
        if len(self.event_types) == 0:
            print "WARNING: Unable to load event types from %s" \
                %(self.miso_events_dir)
        # Load MISO event filters
        self.load_event_filters()
        

    def load_event_filters(self,
                           filter_types=["atleast_inc",
                                         "atleast_exc",
                                         "atleast_sum",
                                         "atleast_const"]):
        """
        Load event type count filters.
        """
        self.event_filters = defaultdict(lambda: defaultdict(int))
        # Load filter settings if they are present
        for event_type in self.event_types:
            if event_type in self.settings_info:
                event_settings = self.settings_info[event_type]
                for filter_type in filter_types:
                    if filter_type not in event_settings:
                        continue
                    # Record the filter count
                    count_filter = int(event_settings[filter_type])
                    self.event_filters[event_type][filter_type] = count_filter

        
    def load_cluster(self):
        """
        Load cluster submission object for the particular
        pipeline settings we were given.
        """
        self.my_cluster = \
          cluster.Cluster(self.cluster_type,
                          self.output_dir,
                          self.logger)


    def __repr__(self):
        repr_str = "MISOWrap(settings=%s, output_dir=%s)" \
          %(self.settings_filename, 
            self.output_dir)
        return repr_str


    def __str__(self):
        return self.__repr__()


def summarize_miso_samples(settings_filename,
                           output_dir):
    """
    Summarize samples in MISO directory.
    """
    misowrap_obj = MISOWrap(settings_filename,
                            output_dir)
    bam_files = misowrap_obj.bam_files
    sample_labels = misowrap_obj.sample_labels
    print "Summarizing MISO output..."
    print "  - Output dir: %s" %(output_dir)
    run_miso_cmd = misowrap_obj.run_miso_cmd
    for sample_label in sample_labels:
        print "sample label: ", sample_label
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
                %(run_miso_cmd,
                  event_dir_path,
                  event_dir_path,
                  sample_basename)
            job_name = "summarize_%s_%s" %(sample_basename,
                                           os.path.basename(event_dirname))
            print "Executing: %s" %(summary_cmd)
            if misowrap_obj.use_cluster:
                misowrap_obj.my_cluster.launch_job(summary_cmd, job_name,
                                                   ppn=1)
            else:
                os.system(summary_cmd)
            

def compare_miso_samples(settings_filename,
                         output_dir):
    """
    Run a MISO samples comparison between all pairs of samples.
    """
    misowrap_obj = MISOWrap(settings_filename,
                            output_dir)
    bam_files = misowrap_obj.bam_files
    sample_labels = misowrap_obj.sample_labels
    read_len = misowrap_obj.read_len
    overhang_len = misowrap_obj.overhang_len
    miso_bin_dir = misowrap_obj.miso_bin_dir
    miso_output_dir = misowrap_obj.miso_outdir
    comparison_groups = misowrap_obj.comparison_groups
    comparisons_dir = misowrap_obj.comparisons_dir
    utils.make_dir(comparisons_dir)
    print "Running MISO comparisons..."
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
            print "Comparing %s %s" %(sample1,
                                      sample2)
            # Directories for each sample
            sample1_dir = os.path.join(miso_output_dir,
                                       sample1)
            sample2_dir = os.path.join(miso_output_dir,
                                       sample2)
            for event_type in misowrap_obj.event_types:
                print "Processing %s..." %(event_type)
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
                    %(misowrap_obj.run_miso_cmd,
                      sample1_event_dir,
                      sample2_event_dir,
                      event_comparisons_dir,
                      sample1,
                      sample2)
                print "Executing: %s" %(compare_cmd)
                if misowrap_obj.use_cluster:
                    misowrap_obj.my_cluster.launch_job(compare_cmd,
                                                       job_name,
                                                       ppn=1)
                else:
                    os.system(compare_cmd)


def run_miso_on_samples(settings_filename, output_dir,
                        use_cluster=True):
    """
    Run MISO on a set of samples.
    """
    misowrap_obj = MISOWrap(settings_filename, output_dir)
    bam_files = misowrap_obj.bam_files
    read_len = misowrap_obj.read_len
    overhang_len = misowrap_obj.overhang_len
    events_dir = misowrap_obj.miso_events_dir
    single_end = False
    if misowrap_obj.insert_lens_dir is not None:
        print "Running in single-end mode..."
        single_end = True
    else:
        insert_lens_dir = misowrap_obj.insert_lens_dir
        print "Running in paired end mode..."
        print "  - Insert length directory: %s" %(insert_lens_dir)
    run_events_analysis = misowrap_obj.run_events_cmd
    event_types_dirs = \
        miso_utils.get_event_types_dirs(misowrap_obj.settings_info)
    miso_settings_filename = misowrap_obj.miso_settings_filename
    for bam_input in bam_files:
        bam_filename, sample_label = bam_input
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
                insert_len_filename = \
                    os.path.join(insert_lens_dir,
                                 "%s.insert_len" %(bam_basename))
                print "Reading paired-end parameters from file..."
                pe_params = miso_utils.read_pe_params(insert_len_filename)
                # Paired-end parameters
                miso_cmd += " --paired-end %.2f %.2f" %(pe_params["mean"],
                                                        pe_params["sdev"])
            # Read length
            miso_cmd += " --read-len %d" %(read_len)
            # Overhang length
            miso_cmd += " --overhang-len %d" %(overhang_len)
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
                misowrap_obj.my_cluster.launch_job(miso_cmd,
                                                   job_name,
                                                   ppn=1)
            else:
                os.system(miso_cmd)


def filter_events(settings_filename,
                  output_dir):
    """
    Output a set of filtered MISO events.
    """
    misowrap_obj = MISOWrap(settings_filename,
                            output_dir)
    print "Filtering MISO events..."
    psi_table = pt.PsiTable(misowrap_obj)
    psi_table.output_filtered_comparisons()

                
def compute_insert_lens(settings_filename,
                        output_dir):
    """
    Compute insert lengths for all samples.
    """
    misowrap_obj = MISOWrap(settings_filename,
                            output_dir)
    const_exons_gff = misowrap_obj.const_exons_gff
    if not os.path.isfile(const_exons_gff):
        print "Error: %s const exons GFF does not exist." \
            %(const_exons_gff)
        sys.exit(1)

    pe_utils_path = misowrap_obj.pe_utils_cmd 
    insert_len_output_dir = os.path.join(output_dir, "insert_lens")
    num_bams = len(misowrap_obj.bam_files)
    
    print "Computing insert lengths for %d files" %(num_bams)
    for bam_filename in misowrap_obj.bam_files:
        print "Processing: %s" %(bam_filename)
        insert_len_cmd = "%s --compute-insert-len %s %s --output-dir %s" \
            %(pe_utils_path,
              bam_filename,
              const_exons_gff,
              insert_len_output_dir)
        print "Executing: %s" %(insert_len_cmd)
        sample_name = os.path.basename(bam_filename)
        job_name = sample_name.split(".bam")[0]
        if misowrap_obj.use_cluster:
            misowrap_obj.my_cluster.launch_job(insert_len_cmd, job_name,
                                               ppn=1)
        else:
            os.system(insert_len_cmd)


def greeting(parser=None):
    print "misowrap: wrapper for running MISO and parsing its results.\n"
    if parser is not None:
        parser.print_help()
            
        
def main():
    from optparse import OptionParser
    parser = OptionParser()
      
    parser.add_option("--run", dest="run", nargs=1, default=None,
                      help="Run MISO on a set of events. "
                      "Takes a settings filename.")
    parser.add_option("--summarize", dest="summarize", nargs=1, default=None,
                      help="Run MISO summarize on a set of samples. "
                      "Takes a settings filename.")
    parser.add_option("--compare", dest="compare", nargs=1, default=None,
                      help="Run MISO sample comparisons on all pairwise "
                      "comparisons. Takes a settings filename.")
    parser.add_option("--filter", dest="filter", nargs=1,
                      default=None,
                      help="Filter a set of MISO events. "
                      "Takes a settings filename.")
    parser.add_option("--compute-insert-lens", dest="compute_insert_lens",
                      nargs=1, default=None,
                      help="Compute insert lengths for a set of BAM files. " 
                      "takes a settings filename.")
    parser.add_option("--output-dir", dest="output_dir", default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    greeting()

    if options.output_dir == None:
        print "Error: need --output-dir.\n"
        parser.print_help()
        sys.exit(1)
        
    output_dir = utils.pathify(options.output_dir)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.run != None:
        settings_filename = utils.pathify(options.run)
        run_miso_on_samples(settings_filename, output_dir)

    if options.summarize != None:
        settings_filename = utils.pathify(options.summarize)
        summarize_miso_samples(settings_filename, output_dir)
        
    if options.compare != None:
        settings_filename = utils.pathify(options.compare)
        compare_miso_samples(settings_filename, output_dir)

    if options.filter != None:
        settings_filename = utils.pathify(options.filter)
        filter_events(settings_filename, output_dir)

    if options.compute_insert_lens != None:
        settings_filename = utils.pathify(options.compute_insert_lens)
        compute_insert_lens(settings_filename, output_dir)
        

if __name__ == '__main__':
    main()

