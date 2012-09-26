##
## Main driver for running pipeline
##

import os
import sys
import time

import utils
import settings
import Pipeline as rna_pipeline


def run_pipeline(settings_filename,
                 output_dir):
    """
    Run pipeline given settings file.
    """
    # Create pipeline instance
    pipeline = rna_pipeline.Pipeline(settings_filename,
                                     output_dir)
    # Run pipeline
    pipeline.run()
    # Summarize the results
    #rna_pipeline.summarize_run()

def run_on_sample(sample_label,
                  settings_filename,
                  output_dir):
    """
    Run pipeline on one particular sample.
    """
    pipeline = rna_pipeline.Pipeline(settings_filename,
                                     output_dir)
    pipeline.run_on_sample(sample_label)


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
#    parser.add_option("--summarize-run", dest="summarize_run", nargs=1,
#                      default=None,
#                      help="Summarize the run.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1,
                      default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.settings == None:
        print "Error: need --settings"
        sys.exit(1)
        
    settings_filename = utils.pathify(options.settings)
    if options.output_dir == None:
        print "Error: need --output-dir"
        sys.exit(1)
        
    output_dir = utils.pathify(options.output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.run:
        run_pipeline(settings_filename,
                     output_dir)

    if options.run_on_sample is not None:
        sample_label = options.run_on_sample
        run_on_sample(sample_label, settings_filename,
                      output_dir)
    

if __name__ == '__main__':
    main()
