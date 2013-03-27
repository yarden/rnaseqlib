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
    def __init__(self, settings_filename, output_dir,
                 logger_label=None):
        self.settings_filename = settings_filename
        self.settings_info = None
        self.logger_label = None
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
        self.comparison_groups = \
            self.settings_info["data"]["comparison_groups"]
        self.logs_outdir = os.path.join(self.output_dir,
                                        "misowrap_logs")
        # Create necessary directories
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
        if self.logger_label is None:
            self.logger_label = "misowrap"
        else:
            self.logger_label = "misowrap_%s" %(logger_label)
        self.logger = utils.get_logger(self.logger_label,
                                       self.logs_outdir)
        # Whether to prefilter MISO events
        # Set general default settings
        if "prefilter_miso" not in settings_info["settings"]:
            # By default, set it so that MISO events are not
            # prefiltered
            settings_info["settings"]["prefilter_miso"] = False
        self.prefilter_miso = \
            self.settings_info["settings"]["prefilter_miso"]
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

