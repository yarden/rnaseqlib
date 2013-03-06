##
## Utilities for running Homer
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils

homer_path = utils.which("findMotifsGenome.pl")


def run_homer(logger, bed_fname, genome, output_dir,
              params):
    """
    Run Homer against an input BED file.

    findMotifsGenome.pl <pos file> <genome> <output directory> 
    """
    if homer_path is None:
        logger.critical("Error: Cannot find or execute Homer program.")
        sys.exit(1)
    params_str =  " ".join(["%s %s" %(p, params[p]) for p in params])
    utils.make_dir(output_dir)
    # If there's a Homer results directory in the target
    # directory, then don't rerun Homer
    if os.path.isdir(os.path.join(output_dir, "homerResults")):
        logger.info("Found Homer results, skipping..")
        return output_dir
    homer_cmd = "%s %s %s %s %s" %(homer_path,
                                   bed_fname,
                                   genome,
                                   output_dir,
                                   params_str)
    logger.info("Calling Homer: ")
    logger.info("Executing: %s" %(homer_cmd))
    t1 = time.time()
    ret_val = os.system(homer_cmd)
    if ret_val != 0:
        logger.critical("Error: Homer call failed.")
        sys.exit(1)
    t2 = time.time()
    logger.info("Homer completed in %.2f minutes" %((t2 - t1)/60.))
    return output_dir


