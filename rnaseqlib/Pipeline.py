import os
import sys
import time
import settings

class Sample:
    """
    Sample to run on.
    """
    def __init__(self, label, seq_filename,
                 settings_info=None):
        self.label = label
        self.seq_filename = seq_filename
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
                 output_dir):
        """
        Initialize pipeline.
        """
        # Load settings file
        self.settings_filename = settings_filename
        # Load settings
        self.sequence_filenames = None
        self.parsed_settings = None
        self.settings_info = None
        self.data_type = None
        # Paired-end or not
        self.is_paired_end = None
        # Load samples
        self.sample_to_group = None
        self.group_to_samples = None
        self.samples = None
        # Check settings are correct
        self.load_pipeline_settings()
        self.check_settings()
        self.load_pipeline_samples()

        
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
        # Compile flags
        print "Loaded pipeline settings (source: %s)." \
            %(self.settings_filename)
        

    def load_samples(self):
        """
        Load samples. If paired-end, load their groups
        information.
        """
        if self.is_paired_end:
            # If paired-end, load the groups information
            self.load_groups()


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
        
        return sample
    

    def run_qc(self, sample):
        """
        Run QC.
        """
        return sample

    
    def run_analysis(self, sample):
        return sample
        
