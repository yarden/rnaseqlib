##
## Utilities for running MEME
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils


def get_meme_default_params():
    """
    Get default parameters for calling MEME.
    """
    params = {
        # Use -text output by default (not HTML)
        "-text": ""
    }
    return params
    

def run_meme(input_fasta_fname, output_dir,
             user_params=None):
    """
    Run MEME against an input FASTA file.
    """
    # Get default parameters for MEME
    params = get_meme_default_params()
    if user_params is not None:
        # Update parameters with user given parameters, if any
        params.update(user_params)
    print "Running MEME with parameters: ", params
    # Check if MEME is available
    meme_path = utils.which("meme")
    if meme_path is None:
        print "Error: Cannot find or execute \'meme\' program."
        sys.exit(1)
    params_str =  " ".join(["%s %s" %(p, params[p]) for p in params])
    meme_cmd = "%s %s %s" %(meme_path,
                            params_str,
                            input_fasta_fname)
    ret_val = os.system(meme_cmd)
    if ret_val != 0:
        print "Error: MEME call failed."
        sys.exit(1)
    return output_dir
