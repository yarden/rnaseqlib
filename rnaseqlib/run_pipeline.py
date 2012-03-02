##
## Main driver for running pipeline
##

import os
import sys
import time

import utils
import settings

def run_pipeline(settings_filename, output_dir):
    """
    Run pipeline given settings file.
    """
    settings_info, parsed_settings = \
        settings.load_settings(settings_filename)


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run", dest="run", nargs=1, default=None,
                      help="Run pipeline.")
    parser.add_option("--settings", dest="settings", nargs=1, default=None,
                      help="Settings filename.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.settings_filename == None:
        print "Error: need --settings".
        sys.exit(1)
        
    settings_filename = utils.pathify(options.settings)
    if options.output_dir == None:
        print "Error: need --output-dir"
        sys.exit(1)
        
    output_dir = utils.pathify(output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if options.run != None:
        run_pipeline(settings_filename, output_dir)
    

if __name__ == '__main__':
    main()
