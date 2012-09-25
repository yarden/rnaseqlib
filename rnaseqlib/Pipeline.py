import os
import sys
import time
import settings

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.mapping.mapper_wrappers as mapper_wrappers
import rnaseqlib.ribo.ribo_utils as ribo_utils
import cluster_utils.cluster as cluster


class Sample:
    """
    Sample to run on.
    """
    def __init__(self, label, seq_filename,
                 settings_info=None):
        self.label = label
        # Original sequence file associated with sample
        self.seq_filename = seq_filename
        # Post-processed sequence file, e.g. trimmed reads
        # for CLIP or Ribo-Seq
        # By default, just identical to raw sequence filename
        self.reads_filename = seq_filename
        # Bowtie mapping for sample
        self.bowtie_filename = None
        # Tophat mapping for sample
        self.tophat_filename = None
        # BAM filename
        self.bam_filename = None
        self.settings_info = settings_info
        self.group = None
        self.sample_type = None
        if self.settings_info is not None:
            self.sample_type = self.settings_info["pipeline"]["data_type"]

    def __str__(self):
        return "Sample(%s, %s, %s)" %(self.label,
                                      self.sample_type,
                                      self.seq_filename)


class Pipeline:
    """
    Pipeline for RNA processing.
    """
    def __init__(self,
                 settings_filename,
                 log_output_dir):
        """
        Initialize pipeline.
        """
        # Output directory for logging pipeline activity
        self.log_output_dir = log_output_dir
        # Output directory for actual pipeline output
        self.output_dir = None
        # Load settings file
        self.settings_filename = settings_filename
        # Load settings
        self.sequence_filenames = None
        self.parsed_settings = None
        self.settings_info = None
        self.data_type = None
        # Paired-end or not
        self.is_paired_end = None
        self.sample_to_group = None
        self.group_to_samples = None
        self.samples = []
        # Cluster objects to use
        self.my_cluster = None
        # Check settings are correct
        self.load_pipeline_settings()
        self.load_cluster()
        self.check_settings()
        # Load samples
        self.load_pipeline_samples()
        # Pipeline output subdirectories
        self.pipeline_outdirs = {}
        self.init_outdirs()
        

    def init_outdirs(self):
        """
        Create the output directories for the pipeline.

        Structure is:

        output_dir
          - rawdata: trimmed reads, etc.
          - mapping: mapped data files
          - qc: quality control output
          - analysis: analysis output
        """
        toplevel_dirs = ["rawdata",
                         "mapping",
                         "qc",
                         "analysis"]
        print "Initializing the pipeline output directories."
        utils.make_dir(self.output_dir)
        for dirname in toplevel_dirs:
            dirpath = os.path.join(self.output_dir, dirname)
            print " - Creating: %s" %(dirpath)
            utils.make_dir(dirpath)
            self.pipeline_outdirs[dirname] = dirpath
        

        
    def check_settings(self):
        if (self.settings_info == None) \
            or self.parsed_settings == None:
            print "Error: No settings loaded!"
            sys.exit(1)
        ##
        ## TODO: Error check that the necessary parameters
        ## are given here
        ##
        ## ...

            
    def load_pipeline_settings(self):
        """
        Load the settings filename
        """
        if not os.path.isfile(self.settings_filename):
            print "Error: %s is not a settings filename." \
                %(self.settings_filename)
            sys.exit(1)
        self.settings = settings.load_settings(self.settings_filename)
        self.settings_info, self.parsed_settings = self.settings
        # Determine if we're in paired-end mode
        self.is_paired_end = False
        if self.settings_info["mapping"]["paired"]:
            self.is_paired_end = True
        # Load the sequence files
        self.load_sequence_files()
        # Load the directory where pipeline output should go
        self.output_dir = utils.pathify(self.settings_info["data"]["outdir"])
        # Compile flags
        print "Loaded pipeline settings (source: %s)." \
            %(self.settings_filename)


    def load_cluster(self):
        """
        Load cluster submission object for the particular
        pipeline settings we were given.
        """
        self.my_cluster = cluster.Cluster(self.settings_info)
        

    def load_groups(self, settings):
        """
        If paired-end, set sample groups.
        """
        if (settings == None) or \
            ("sample_groups" not in settings["data"]["sample_groups"]):
            return
        sample_groups = settings["data"]["sample_groups"]
        sample_to_group = {}
        group_to_samples = defaultdict(list)
        # Map sample to its group
        for sample, group in sample_groups:
            sample_to_group[sample] = group
            # Map each group to its samples
            group_to_samples[group].append(sample)
        self.sample_to_group = sample_to_group
        self.group_to_samples = group_to_samples


    def load_sequence_files(self):
        """
        Load sequence files from settings file.
        """
        if self.settings_info is None:
            print "Error: cannot load sequence files if settings " \
                "are not loaded."
            sys.exit(1)
        seq_files = self.settings_info["data"]["sequence_files"]
        # Get the absolute path names, with the prefix input directory,
        # for each sequence file
        sequence_filenames = []
        input_dir = os.path.abspath(self.settings_info["data"]["indir"])
        for seq_entry in seq_files:
            if len(seq_entry) != 2:
                print "Error: Must provide a sequence filename and a " \
                    "sample label for each entry."
                sys.exit(1)
            fname, seq_label = seq_entry
            seq_fname = os.path.join(input_dir, fname)
            if not os.path.isfile(seq_fname):
                print "Error: Cannot find sequence file %s" %(seq_fname)
                sys.exit(1)
            sequence_filenames.append([seq_fname, seq_label])
        self.sequence_filenames = sequence_filenames
        return sequence_filenames
    
        
    def load_pipeline_samples(self):
        """
        Load samples.
        """
        print "Loading pipeline samples..."
        samples = []
        num_seq_files = len(self.sequence_filenames)
        print "  - Total of %d sequence files." %(num_seq_files)
        for seq_entry in self.sequence_filenames:
            seq_filename, sample_label = seq_entry
            # Ensure file exists
            if not os.path.isfile(seq_filename):
                print "Error: %s does not exist!" %(seq_filename)
                sys.exit(1)
            sample = Sample(sample_label,
                            seq_filename,
                            settings_info=self.settings_info)
            samples.append(sample)
        self.samples = samples


    def preprocess_reads(self, sample):
        """
        Pre-process reads.
        """
        print "Preprocessing: %s" %(sample)
        if sample.sample_type == "riboseq":
            # Preprocess riboseq samples by trimming trailing
            # As
            trimmed_filename = ribo_utils.trim_polyA_ends(sample.seq_filename,
                                                          self.pipeline_outdirs["rawdata"])
            # Adjust the trimmed file to be the "reads" sequence file for this
            # sample
            sample.reads_filename = trimmed_filename
        return sample
            
        
    def run(self):
        """
        Run pipeline. 
        """
        print "Running pipeline..."
        num_samples = len(self.samples)
        if num_samples == 0:
            print "Error: No samples to run on."
            sys.exit(1)
        else:
            print "Running on %d samples" %(num_samples)
        # For each sample
        for sample in self.samples:
            print "Processing %s" %(sample)
            # Pre-process the data if needed
            sample = self.preprocess_reads(sample)
            # Map the data
            sample = self.map_reads(sample)
            # Perform QC
            sample = self.run_qc(sample)
            # Run gene expression analysis
            sample = self.run_analysis(sample)

            
    def map_reads(self, sample):
        """
        Map reads.
        """
        print "Mapping reads..."
        mapper = self.settings_info["mapping"]["mapper"]
        mapping_cmd = None
        job_name = "%s_%s" %(sample.label, mapper)
        if mapper == "bowtie":
            print "Mapping sample: %s" %(sample)
            bowtie_path = self.settings_info["mapping"]["bowtie_path"]
            index_filename = self.settings_info["mapping"]["bowtie_index"]
            output_filename = "%s.bowtie" %(os.path.join(self.pipeline_outdirs["mapping"],
                                                         sample.label))
            bowtie_options = self.settings_info["mapping"]["bowtie_options"]
            # Number of mismatches to use in mapping
            # Optional bowtie arguments
            mapping_cmd, bowtie_output_filename = \
                  mapper_wrappers.get_bowtie_mapping_cmd(bowtie_path,
                                                         sample.reads_filename,
                                                         index_filename,
                                                         output_filename,
                                                         bowtie_options=bowtie_options)
            # Record the bowtie output filename for this sample
            sample.bowtie_filename = bowtie_output_filename
            print "Executing: %s" %(mapping_cmd)
            self.launch_and_wait(mapping_cmd, job_name,
                                 unless_exists=output_filename)
            job_id = self.my_cluster.launch_job(mapping_cmd, job_name,
                                                unless_exists=output_filename)
        elif mapper == "tophat":
            raise Exception, "Not implemented yet."
        else:
            print "Error: unsupported mapper %s" %(mapper)
            sys.exit(1)
        return sample
    

    def run_qc(self, sample):
        """
        Run QC.
        """
        return sample

    
    def run_analysis(self, sample):
        return sample
        
