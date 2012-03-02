##
## Main driver for running pipeline
##

import os
import sys
import time

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

    output_dir = options.output_dir
    if output_dir != None:
        output_dir = pathify(output_dir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    

if __name__ == '__main__':
    main()
