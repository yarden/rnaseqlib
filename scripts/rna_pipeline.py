##
## Main driver for running RNA-Seq pipeline
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.settings as settings
import rnaseqlib.Pipeline as rna_pipeline
import rnaseqlib.RNABase as rna_base

def run_pipeline(settings_filename,
                 output_dir):
    """
    Run pipeline on all samples given settings file.
    """
    # Create pipeline instance
    pipeline = rna_pipeline.Pipeline(settings_filename,
                                     output_dir)
    # Run pipeline
    pipeline.run()

    
def run_on_sample(sample_label,
                  settings_filename,
                  output_dir):
    """
    Run pipeline on one particular sample.
    """
    pipeline = rna_pipeline.Pipeline(settings_filename,
                                     output_dir)
    pipeline.run_on_sample(sample_label)


def check_requirements():
    print "Checking that all required programs are available..."
    # Utilities that need to be on path for pipeline to run
    REQUIRED_PROGRAMS = [# UCSC utils
                         "genePredToGtf",
                         # Tophat/Bowtie
                         "bowtie-build",
                         "tophat",
                         # Bedtools
                         "intersectBed",
                         "tagBam",
                         # Related utils
                         "gtf2gff3.pl",
                         # Unix utils
                         "gunzip",
                         "wget",
                         "cat",
                         "zcat",
                         "cut"]
    found_all = True
    for program in REQUIRED_PROGRAMS:
        if utils.which(program) is None:
            print "WARNING: Cannot find required program \'%s\' " \
                  "on your path.  Please install it or add it. " \
                  "to your path if already installed." %(program)
            if program == "genePredToGtf":
                genePredToGtf_msg()
            print "  - Proceeding anyway..."
            found_all = False
    if found_all:
        print "Found all required programs."


def genePredToGtf_msg():
    print "To install genePredToGtf, download the executable for your OS " \
          "from: "
    print "http://hgdownload.cse.ucsc.edu/admin/exe/"
    print "And place the executable in your PATH."


def initialize_pipeline(genome,
                        output_dir):
    """
    Initialize the pipeline.
    """
    # Check for required programs
    check_requirements()
    base_obj = rna_base.RNABase(genome, output_dir)
    base_obj.initialize()


def greeting():
    print "rnaseqlib: a lightweight pipeline for RNA-Seq analysis"
    print "=" * 10 
    print "For help, pass --help\n"
    
    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run", dest="run", action="store_true",
                      default=False,
                      help="Run pipeline.")
    parser.add_option("--run-on-sample", dest="run_on_sample", nargs=1, default=None,
                      help="Run on a particular sample. Takes as input the sample label.")
    parser.add_option("--settings", dest="settings", nargs=1,
                      default=None,
                      help="Settings filename.")
    parser.add_option("--init", dest="initialize", nargs=1, default=None,
                      help="Initialize the pipeline. Takes as input a genome, "
                      "e.g. mm9 or hg18")
    parser.add_option("--output-dir", dest="output_dir", nargs=1,
                      default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    greeting()

    if options.output_dir == None:
        print "Error: need --output-dir"
        sys.exit(1)
        
    output_dir = utils.pathify(options.output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    settings_filename = None
    if options.run:
        if options.settings == None:
            # Running of pipeline requires settings filename
            print "Error: need --settings"
            sys.exit(1)
        settings_filename = utils.pathify(options.settings)
        run_pipeline(settings_filename,
                     output_dir)

    if options.run_on_sample is not None:
        if options.settings == None:
            # Running of pipeline requires settings filename
            print "Error: need --settings"
            sys.exit(1)
        settings_filename = utils.pathify(options.settings)
        sample_label = options.run_on_sample
        run_on_sample(sample_label, settings_filename,
                      output_dir)

    if options.initialize is not None:
        genome = options.initialize
        initialize_pipeline(genome,
                            output_dir)
    

if __name__ == '__main__':
    main()
