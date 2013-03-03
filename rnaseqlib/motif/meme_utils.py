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
        "-text": "",
        "-maxsize": "100000000"
    }
    return params
    

def run_meme(logger, input_fasta_fname, output_dir,
             meme_params=None):
    """
    Run MEME against an input FASTA file.
    """
    # Get default parameters for MEME
    params = get_meme_default_params()
    # Set output directory for MEME
    params.update({"-o": output_dir})
    if meme_params is not None:
        # Update parameters with user-given parameters, if any
        params.update(meme_params)
    # Check if MEME is available
    meme_path = utils.which("meme")
    if meme_path is None:
        logger.critical("Error: Cannot find or execute \'meme\' program.")
        sys.exit(1)
    params_str =  " ".join(["%s %s" %(p, params[p]) for p in params])
    fasta_basename = \
        os.path.basename(input_fasta_fname).rsplit(".", 1)[0]
    meme_output_fname = \
        os.path.join(output_dir, "%s.meme" %(fasta_basename))
    if os.path.isfile(meme_output_fname):
        logger.info("Found MEME file %s, skipping..." \
                    %(meme_output_fname))
        return meme_output_fname
    meme_cmd = "%s %s %s &> %s" %(meme_path,
                                  params_str,
                                  input_fasta_fname,
                                  meme_output_fname)
    logger.info("Calling MEME: ")
    logger.info("Executing: %s" %(meme_cmd))
    t1 = time.time()
    ret_val = os.system(meme_cmd)
    if ret_val != 0:
        logger.critical("Error: MEME call failed.")
        sys.exit(1)
    t2 = time.time()
    logger.info("MEME completed in %.2f minutes" %((t2 - t1)/60.))
    return meme_output_fname
