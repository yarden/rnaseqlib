import os
import sys
import time
import glob
import settings

import pysam
import pandas

from collections import defaultdict

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.rpkm.rpkm_utils as rpkm_utils
import rnaseqlib.mapping.mapper_wrappers as mapper_wrappers
import rnaseqlib.mapping.bedtools_utils as bedtools_utils
import rnaseqlib.ribo.ribo_utils as ribo_utils
import rnaseqlib.QualityControl as qc
import rnaseqlib.RNABase as rna_base

# Import all paths
from rnaseqlib.paths import *

import cluster_utils.cluster as cluster

class Sample:
    """
    A sample to run on. For paired-end, represents a
    pair of samples. For single end, represents a
    single sample.
    """
    def __init__(self, label, rawdata):
        self.label = label
        # QC information for this sample
        self.qc = None
        # For paired-end samples, samples_info is a list
        # of samples.
        # For single-end samples, it is just one object.
        self.rawdata = rawdata
        self.paired = False
        self.sample_type = None
        # Bowtie mapping for sample
        self.bowtie_filename = None
        # Tophat mapping for sample
        self.tophat_filename = None
        # Directory where processed BAM dirs are
        self.processed_bam_dir = None
        # BAM filename
        self.bam_filename = None
        # Unique BAM filename
        self.unique_bam_filename = None
        # rRNA subtracted BAM filename
        self.ribosub_bam_filename = None
        # Sample's RPKM directory
        self.rpkm_dir = None
        # RPKM tables for the sample
        self.rpkm_tables = defaultdict(lambda: None)
        # Record if a sample is grouped
        if type(self.rawdata) == list:
            self.paired = True
            self.sample_type = self.rawdata[0].sample_type
        else:
            self.sample_type = self.rawdata.sample_type
            

    def __repr__(self):
        if self.paired:
            samples_str = ",".join([s.label for s in self.rawdata])
        else:
            samples_str = self.rawdata.label
        return "Sample(%s, samples=%s)" \
            %(self.label, samples_str)
    

    def __str__(self):
        return self.__repr__()
        

class SampleRawdata:
    """
    Representation a of a set of rawdata files related to a single
    sample. This represents the raw data: like fastq file.
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
        self.settings_info = settings_info
        self.sample_type = None
        if self.settings_info is not None:
            self.sample_type = self.settings_info["pipeline"]["data_type"]
            

    def __str__(self):
        return "SampleRawdata(%s, %s, %s)" %(self.label,
                                             self.sample_type,
                                             self.seq_filename)


class Pipeline:
    """
    Pipeline for RNA processing.
    """
    def __init__(self,
                 settings_filename,
                 log_output_dir,
                 curr_sample=None):
        """
        Initialize pipeline.
        """
        # If invoked to run on particular sample
        self.curr_sample = curr_sample
        self.genome = None
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
        # Directory where pipeline init files are
        self.init_dir = None
        # Paired-end or not
        self.is_paired_end = None
        self.sample_to_group = None
        self.group_to_samples = None
        self.samples = []
        # Cluster objects to use
        self.my_cluster = None
        # Check settings are correct
        self.load_pipeline_settings()
        # Pipeline output subdirectories
        self.pipeline_outdirs = {}
        # RPKM directory for teh pipeline
        self.rpkm_dir = None
        # QC objects for each sample in pipeline
        self.qc_objects = {}
        # Top-level output dirs
        self.toplevel_dirs = ["rawdata",
                              "mapping",
                              "qc",
                              "analysis",
                              "logs"]
        self.init_outdirs()
        pipeline_log_name = "Pipeline"
        if self.curr_sample is not None:
            pipeline_log_name = "Pipeline.%s" %(self.curr_sample)
        self.logger = utils.get_logger(pipeline_log_name,
                                       self.pipeline_outdirs["logs"])
        self.load_cluster()
        ## Load RNA Base: object storing all the relevant
        ## initialization information
        self.rna_base = None
        self.load_rna_base()
        ## Load samples
        self.load_pipeline_samples()
        ## Initialize QC for samples
        # QC header: order of QC fields to be outputted
        self.qc_header = []
        self.init_qc()
        self.na_val = "NA"


    def load_rna_base(self):
        """
        Load pipeline RNA base.
        """
        self.rna_base = rna_base.RNABase(self.genome,
                                         None,
                                         from_dir=self.init_dir)
        # Load genes information: tables only, without parsing
        # into gene objects
        self.rna_base.load_gene_tables(tables_only=True)
        

    def init_qc(self):
        """
        Make a QC object for each sample in pipeline.
        """
        # Create a QualityControl object for each sample in pipeline
        if len(self.samples) == 0:
            print "WARNING: No samples to create QC objects for."
            return
        self.load_qc()
        

    def load_qc(self):
        """
        Load QC information from samples.
        """
        for sample in self.samples:
            self.qc_objects[sample.label] = qc.QualityControl(sample,
                                                              self)
            sample.qc = self.qc_objects[sample.label]
        # Retrieve QC header: get the header from first sample
        self.qc_header = self.qc_objects[self.samples[0].label].qc_header


    def load_rpkms(self):
        """
        Load RPKM information from samples.
        """
        for sample_num, sample in enumerate(self.samples):
            # Get mapping from table names to loaded RPKM tables
            self.samples[sample_num].rpkm_tables = \
                rpkm_utils.load_sample_rpkms(sample,
                                             self.rna_base)
        # Load RPKMs into combined RPKM tables for all samples
        # Create a mapping from table name to RPKM DataFrame
        # that includes all samples
        self.rpkm_tables = defaultdict(lambda: None)
        for table_name in self.rna_base.rpkm_table_names:
            curr_sample = self.samples[0]
            # If the RPKM table is not available, skip it
            if curr_sample.rpkm_tables[table_name] is None: continue
            # If we have only one sample, make add it to the set of
            # RPKM tables by itself and continue to next table
            if len(self.samples) == 1:
                self.rpkm_tables[table_name] \
                    = curr_sample.rpkm_tables[table_name]
                continue
            # For multiple samples: collect each sample's RPKM table
            # for the current table (e.g. ensGene)
            combined_rpkm_table = curr_sample.rpkm_tables[table_name]
            for next_sample in self.samples[1:]:
                # Merge with next sample's RPKM table
                next_sample_rpkm = next_sample.rpkm_tables[table_name]
                combined_rpkm_table = \
                    pandas.merge(combined_rpkm_table,
                                 next_sample_rpkm,
                                 # Merge on common columns of
                                 # gene ID and exons
                                 on=["gene_id", "exons",
                                     "gene_symbol", "gene_desc"])
            # Record the combined RPKM table
            self.rpkm_tables[table_name] = combined_rpkm_table
        

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
        print "Initializing the pipeline output directories."
        utils.make_dir(self.output_dir)
        # Subdirectories of toplevel subdirs
        self.toplevel_subdirs = defaultdict(list)
        self.toplevel_subdirs["analysis"] = ["rpkm", "insert_lens"]
        for dirname in self.toplevel_dirs:
            dirpath = os.path.join(self.output_dir, dirname)
            print " - Creating: %s" %(dirpath)
            utils.make_dir(dirpath)
            self.pipeline_outdirs[dirname] = dirpath
            for subdir_name in self.toplevel_subdirs[dirname]:
                subdir_path = os.path.join(dirpath, subdir_name)
                utils.make_dir(subdir_path)
        # Variables storing commonly accessed directories
        self.rpkm_dir = os.path.join(self.pipeline_outdirs["analysis"],
                                     "rpkm")

            
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
        self.genome = self.settings_info["mapping"]["genome"]
        # Determine if we're in paired-end mode
        self.is_paired_end = False
        if self.settings_info["mapping"]["paired"]:
            self.is_paired_end = True
        # Load the sequence files
        self.load_sequence_files()
        # Load the directory where pipeline output should go
        self.output_dir = utils.pathify(self.settings_info["data"]["outdir"])
        print "Loaded pipeline settings (source: %s)." \
            %(self.settings_filename)
        # Pipeline init directory
        self.init_dir = \
            os.path.join(self.settings_info["pipeline-files"]["init_dir"])
        # Loading group information if there is any
        self.load_groups()
        

    def load_groups(self):
        """
        If paired-end, set sample groups.
        """
        if not self.settings_info["mapping"]["paired"]:
            return
        if "sample_groups" not in self.settings_info["data"]:
            print "Error: In paired-end mode, but cannot find \'sample_groups\' "\
                  "parameter. Please set \'sample_groups\'."
            sys.exit(1)
        sample_groups = self.settings_info["data"]["sample_groups"]
        sample_to_group = {}
        group_to_samples = defaultdict(list)
        # Map sample to its group
        for group, samples in sample_groups:
            group = str(group)
            samples = map(str, samples)
            for curr_sample in samples:
                sample_to_group[curr_sample] = group
                group_to_samples[group].append(curr_sample)
        self.sample_to_group = sample_to_group
        self.group_to_samples = group_to_samples
        # Tell each sample what group it is in
        for sample_obj in self.samples:
            sample_obj.group = self.sample_to_group[sample_obj.label]
        


    def load_cluster(self):
        """
        Load cluster submission object for the particular
        pipeline settings we were given.
        """
        self.my_cluster = \
            cluster.Cluster(self.settings_info["mapping"]["cluster_type"],
                            self.output_dir,
                            self.logger)
        

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
        input_dir = utils.pathify(self.settings_info["data"]["indir"])
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


    def get_samples_rawdata(self):
        """
        Load rawdata information related to each sample.
        """
        # Mapping from label to samples info
        all_samples_rawdata = []
        for seq_entry in self.sequence_filenames:
            seq_filename, sample_label = seq_entry
            # Ensure file exists
            if not os.path.isfile(seq_filename):
                self.logger.critical("Error: %s does not exist!" %(seq_filename))
                print "Error: %s does not exist!" %(seq_filename)
                sys.exit(1)
            sample_rawdata = SampleRawdata(sample_label,
                                           seq_filename,
                                           settings_info=self.settings_info)
            all_samples_rawdata.append(sample_rawdata)
        return all_samples_rawdata
    
        
    def load_pipeline_samples(self):
        """
        Load samples.
        """
        print "Loading pipeline samples..."
        samples = []
        # Get samples information
        all_samples_rawdata = self.get_samples_rawdata()
        # Mapping from labels to sample info
        samples_rawdata_by_label = dict([(s.label, s) for s in all_samples_rawdata])
        # If paired-end, also load sample groups information
        if self.is_paired_end:
            self.load_groups()
            for group_label, samples_in_group in self.group_to_samples.iteritems():
                # Get all the samples in the group
                group_samples = [samples_rawdata_by_label[label] \
                                 for label in samples_in_group]
                # Create a sample consisting of all the samples info
                # objects in this group
                sample = Sample(group_label, group_samples)
                samples.append(sample)
        else:
            # If single-end, then each sample consists of just one
            # sample rawdata object
            for sample_rawdata in all_samples_rawdata:
                sample = Sample(sample_rawdata.label, sample_rawdata)
                samples.append(sample)
        self.samples = samples
        # Tell each sample locations of its various output directories
        # (e.g. where its RPKM directory is)
        for sample in self.samples:
            sample.rpkm_dir = os.path.join(self.rpkm_dir, sample.label)
            

    def preprocess_reads(self, sample):
        """
        Pre-process reads.
        """
        print "Preprocessing: %s" %(sample)
        if sample.sample_type == "riboseq":
            # Preprocess riboseq samples by trimming trailing
            # As
            trimmed_filename = ribo_utils.trim_polyA_ends(sample.rawdata.seq_filename,
                                                          self.pipeline_outdirs["rawdata"])
            # Adjust the trimmed file to be the "reads" sequence file for this
            # sample
            sample.rawdata.reads_filename = trimmed_filename
        else:
            print "WARNING: Do not know how to pre-process type %s samples." \
                %(sample.sample_type)
            self.logger.info("Do not know how to pre-process type %s samples" \
                             %(sample.sample_type))
        return sample


    def get_sample_by_label(self, label):
        """
        Return a sample by its label.
        """
        for sample in self.samples:
            if sample.label == label:
                return sample
        return None
            
            
    def run_on_samples(self):
        samples_job_ids = []
        self.logger.info("Running on samples..")
        for sample in self.samples:
            print "Processing sample %s" %(sample)
            job_name = "pipeline_run_%s" %(sample.label)
            sample_cmd = "python %s --run-on-sample %s --settings %s --output-dir %s" \
                %(PIPELINE_RUN_SCRIPT,
                  sample.label,
                  self.settings_filename,
                  self.output_dir)
            self.logger.info("Executing: %s" %(sample_cmd))
            print "Executing: %s" %(sample_cmd)
            job_id = self.my_cluster.launch_job(sample_cmd, job_name)
            self.logger.info("Job launched with ID %s" %(job_id))
            samples_job_ids.append(job_id)
        return samples_job_ids
            
        
    def run(self, label=None):
        """
        Run pipeline. 
        """
        self.logger.info("Running pipeline.")
        print "Running pipeline..."
        num_samples = len(self.samples)
        if num_samples == 0:
            print "Error: No samples to run on."
            sys.exit(1)
        else:
            print "Running on %d samples" %(num_samples)
        # Job IDs for each sample
        job_ids = self.run_on_samples()
        # Wait until all jobs completed 
        self.my_cluster.wait_on_jobs(job_ids)
        # Compile all the QC results
        self.compile_qc_output()
        # Compile all the analysis results
        self.compile_analysis_output()
        # Signal completion
        self.logger.info("Run completed!")


    def run_on_sample(self, label):
        self.logger.info("Running on sample: %s" %(label))
        self.logger.info("Retrieving sample...")
        # Fetch the sample by its label
        sample = self.get_sample_by_label(label)
        if sample is None:
            self.logger.info("Cannot find sample %s! Exiting.." \
                             %(label))
            print "Error: Cannot find sample %s" %(label)
            sys.exit(1)
        # Pre-process the data if needed
        self.logger.info("Preprocessing reads")
        sample = self.preprocess_reads(sample)
        # Map the data
        self.logger.info("Mapping reads")
        sample = self.map_reads(sample)
        # Perform QC
        self.logger.info("Running QC")
        sample = self.run_qc(sample)
        # Run gene expression analysis
        self.logger.info("Running analysis")
        sample = self.run_analysis(sample)

            
    def map_reads(self, sample):
        """
        Map reads using a read mapper.

        Also creates a BAM file containing only the uniquely mapped
        reads, which is used for downstream analyses (e.g. in QC).
        """
        self.logger.info("Mapping reads for sample: %s" %(sample.label))
        print "Mapping reads..."
        mapper = self.settings_info["mapping"]["mapper"]
        mapping_cmd = None
        job_name = "%s_%s" %(sample.label, mapper)
        print "Mapping sample: %s" %(sample)
        print "  - mapper: %s" %(mapper)
        self.logger.info("Mapper: %s" %(mapper))
        if mapper == "bowtie":
            bowtie_path = self.settings_info["mapping"]["bowtie_path"]
            index_filename = self.settings_info["mapping"]["bowtie_index"]
            output_filename = "%s" %(os.path.join(self.pipeline_outdirs["mapping"],
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
            sample.bam_filename = sample.bowtie_filename
            self.my_cluster.launch_and_wait(mapping_cmd, job_name,
                                            unless_exists=output_filename)
        elif mapper == "tophat":
            tophat_path = self.settings_info["mapping"]["tophat_path"]
            sample_mapping_outdir = os.path.join(self.pipeline_outdirs["mapping"],
                                                 sample.label)
            print "Creating: %s" %(sample_mapping_outdir)
            utils.make_dir(sample_mapping_outdir)
            tophat_cmd, tophat_outfilename = \
                mapper_wrappers.get_tophat_mapping_cmd(tophat_path,
                                                       sample,
                                                       sample_mapping_outdir,
                                                       self.settings_info)
            print "Executing: %s" %(tophat_cmd)
            # Check that Tophat file does not exist
            self.my_cluster.launch_and_wait(tophat_cmd, job_name,
                                            unless_exists=tophat_outfilename)
            sample.bam_filename = tophat_outfilename
        else:
            print "Error: unsupported mapper %s" %(mapper)
            sys.exit(1)
        ##
        ## Post processing of BAM reads
        ##
        # Index the main BAM file
        self.index_bam(sample.bam_filename)
        # Create a directory for processed BAMs
        sample.processed_bam_dir = os.path.join(self.pipeline_outdirs["mapping"],
                                                sample.label,
                                                "processed_bams")
        utils.make_dir(sample.processed_bam_dir)
        # Get the uniquely mapping reads
        sample.unique_bam_filename = self.get_unique_reads(sample)
        # Get the ribo-subtracted mapping reads
        sample.ribosub_bam_filename = self.get_ribosub_bam_reads(sample)
        # Sort and index the ribo-subtracted mapped reads BAM
        sample.bam_filename = self.sort_and_index_bam(sample.bam_filename)
        # Sort and index the unique BAM reads
        sample.unique_bam_filename = self.sort_and_index_bam(sample.unique_bam_filename)
        # Sort and index the ribosubtracted BAM reads
        sample.ribosub_bam_filename = self.sort_and_index_bam(sample.ribosub_bam_filename)
        return sample


    def index_bam(self, bam_filename):
        """
        Index a BAM filename if it's not already indexed.
        """
        index_cmd = "samtools index %s" %(bam_filename)            
        index_filename = "%s.bai" %(bam_filename)
        print "Indexing %s" %(bam_filename)
        if not os.path.isfile(index_filename):
            os.system(index_cmd)


    def sort_and_index_bam(self, bam_filename):
        """
        Sort and index the BAM for the sample.

        Return the filename of the sorted, index filename.
        """
        self.logger.info("Sort and indexing BAM: %s" \
                         %(bam_filename))
        bam_basename = os.path.basename(bam_filename).split(".bam")[0]
        sorted_bam_filename = os.path.join(os.path.dirname(bam_filename),
                                           "%s.sorted" %(bam_basename))
        print "Sorting %s as %s" %(bam_filename,
                                   sorted_bam_filename)
        sort_cmd = "samtools sort %s %s" %(bam_filename,
                                           sorted_bam_filename)
        job_name = "sorted_bam_%s" %(bam_basename)
        expected_bam_filename = "%s.bam" %(sorted_bam_filename)
        if not os.path.isfile(expected_bam_filename):
            os.system(sort_cmd)
        # Index the sorted BAM
        self.index_bam(expected_bam_filename)
        return expected_bam_filename



    def get_ribosub_bam_reads(self, sample, chr_ribo="chrRibo"):
        """
        Subtract the rRNA-mapping reads away from the BAM file
        and create a new BAM file.
        """
        self.logger.info("Getting rRNA-subtracted BAM file for %s" \
                         %(sample.label))
        mapped_reads = pysam.Samfile(sample.bam_filename, "rb")
        if not sample.bam_filename.endswith(".bam"):
            self.logger.critical("BAM %s file does not end in .bam" \
                                 %(sample.bam_filename))
        bam_basename = os.path.basename(sample.bam_filename)
        ribosub_bam_filename = os.path.join(sample.processed_bam_dir,
                                           "%s.ribosub.bam" %(bam_basename[0:-4]))
        print "Getting rRNA-subtracted reads for %s" %(sample.label)
        print "  - Output file: %s" %(ribosub_bam_filename)
        if not os.path.isfile(ribosub_bam_filename):
            # Get the ribosomal rRNA mapping reads
            ribo_read_ids = {}
            ribo_reads = mapped_reads.fetch(reference=chr_ribo,
                                            start=None,
                                            end=None)
            for ribo_read in ribo_reads:
                ribo_read_ids[ribo_read.qname] = True
            ribosub_bam = pysam.Samfile(ribosub_bam_filename, "wb",
                                        # Use original file's headers
                                        template=mapped_reads)
            for read in mapped_reads:
                # If the read has any mapping to rRNA, then
                # skip it
                if read.qname in ribosub_bam:
                    continue
                ribosub_bam.write(read)
            ribosub_bam.close()
        else:
            print "Found %s. Skipping.." %(ribosub_bam_filename)
        mapped_reads.close()
        return ribosub_bam_filename
        

    def get_unique_reads(self, sample):
        """
        Get only the uniquely mapping reads from the reads BAM
        file and put them in a new file.
        """
        self.logger.info("Getting unique reads for %s" %(sample.label))
        mapped_reads = pysam.Samfile(sample.bam_filename, "rb")
        if not sample.bam_filename.endswith(".bam"):
            self.logger.critical("BAM %s file does not end in .bam" \
                                 %(sample.bam_filename))
        bam_basename = os.path.basename(sample.bam_filename)
        unique_bam_filename = os.path.join(sample.processed_bam_dir,
                                           "%s.unique.bam" %(bam_basename[0:-4]))
        print "Getting unique reads for %s" %(sample.label)
        print "  - Output file: %s" %(unique_bam_filename)
        if not os.path.isfile(unique_bam_filename):
            # Get unique reads if file for them does not already
            # exist
            unique_reads = pysam.Samfile(unique_bam_filename, "wb",
                                         # Use original file's headers
                                         template=mapped_reads)
            unique_bam_filename
            num_unique = 0
            for read in mapped_reads:
                # Keep only reads with 'NH' tag equal to 1
                if ("NH", 1) in read.tags:
                    unique_reads.write(read)
            unique_reads.close()
        else:
            print "Found %s. Skipping.." %(unique_bam_filename)
        mapped_reads.close()
        return unique_bam_filename
    

    def run_qc(self, sample):
        """
        Run QC for this sample.
        """
        self.logger.info("Running QC on %s" %(sample.label))
        print "Running QC on sample: %s" %(sample.label)
        # Retrieve QC object for sample
        qc_obj = self.qc_objects[sample.label]
        if qc_obj.qc_loaded:
            self.logger.info("QC objects already loaded.")
            # Don't load QC information if it already exists
            print "  - QC objects already loaded from file."
        else:
            # Run QC metrics
            qc_obj.compute_qc()
            # Output QC to file
            qc_obj.output_qc()
        sample.qc = qc_obj
        return sample
        

    def compile_qc_output(self):
        """
        Compile QC output for all samples.
        """
        self.logger.info("Compiling QC output for all samples...")
        print "Compiling QC output for all samples..."
        # Load QC information from existing files
        self.load_qc()
        # Get a compiled object representing the QC
        # for all samples in the pipeline
        qc_stats = qc.QCStats(self.samples,
                              self.qc_header,
                              self.qc_objects)
        qc_stats.compile_qc()
        qc_output_filename = os.path.join(self.pipeline_outdirs["qc"],
                                          "qc_stats.txt")
        print "  - Outputting QC to: %s" %(qc_output_filename)
        self.logger.info("Outputting QC to: %s" %(qc_output_filename))
        qc_stats.to_csv(qc_output_filename)


    def compile_analysis_output(self):
        """
        Compile analysis output for all samples.
        """
        self.logger.info("Compiling analysis output for all samples...")
        print "Compiling analysis output for all samples..."
        # Compile RPKM results
        self.compile_rpkms_output()


    def compile_rpkms_output(self):
        """
        Compile and output RPKMs for all samples.
        """
        # Load RPKMs for all samples
        self.load_rpkms()
        # Output RPKM tables
        # Order in which table columns should be serialized:
        # Gene ID first, followed by gene symbol, the counts for each sample,
        # followed by the exons used in the calculation and the
        # gene description
        fieldnames = ["gene_id", "gene_symbol"]
        fieldnames.extend(["rpkm_%s" %(sample.label) \
                           for sample in self.samples])
        fieldnames.extend(["counts_%s" %(sample.label) \
                           for sample in self.samples])
        fieldnames.extend(["gene_desc", "exons"])
        for table_name, rpkm_table in self.rpkm_tables.iteritems():
            if rpkm_table is None: continue
            rpkm_table_filename = os.path.join(self.rpkm_dir,
                                               "%s.rpkm.txt" %(table_name))
            rpkm_table.to_csv(rpkm_table_filename,
                              cols=fieldnames,
                              na_rep=self.na_val,
                              sep="\t",
                              index=False)


    def output_rpkms(self, sample):
        """
        Output RPKMs.
        """
        print "Outputting RPKMs for sample: %s" %(sample.label)
        sample_rpkm_outdir = os.path.join(self.rpkm_dir, sample.label)
        rpkm_tables = rpkm_utils.output_rpkm(sample,
                                             sample_rpkm_outdir,
                                             self.settings_info,
                                             self.rna_base,
                                             self.logger)
        return sample

    
    def run_analysis(self, sample):
        """
        Run analysis on a sample.
        """
        # Compute RPKMs
        self.logger.info("Computing RPKMs for sample: %s" %(sample.label))
        self.output_rpkms(sample)
        return sample
        
