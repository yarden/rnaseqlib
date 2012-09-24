##
## Utilities for running on cluster
##
import os
import sys
import time

import rnaseqlib.cluster_utils as cluster_utils
from cluster_utils import Mybsub, Mypbm, Mysge

class Cluster:
    """
    Cluster submission.
    """
    def __init__(self, settings_info,
                 supported_types=["bsub", "qsub"]):
        self.cluster_type = settings_info["mapping"]["cluster_type"]
        if self.cluster_type not in supported_types:
            print "Error: unsupported cluster type %s" %(self.cluster_type)
            sys.exit(1)
            

    def launch_job(self, cmd, job_name,
                   unless_exists=None):
        """
        Launch job on cluster and return a job id.

        if unless_exists flag is given, do not execute command
        if the given filename path exists.
        
        Wrapper to Mysge/Mypbm/Mybsub
        """
        job_id = None
        script_options = {}
        if (unless_exists is not None) and \
            os.path.isfile(unless_exists):
            print "launch_job: SKIPPING %s since %s exists." \
                %(cmd, unless_exists)
            return job_id
        if self.cluster_type == "bsub":
            job_id = Mybsub.launchJob(cmd, job_name,
                                      script_options,
                                      queue_type="normal")
        if job_id is None:
            print "WARNING: Job %s not submitted." %(job_name)
        return job_id
        

    def wait_on_job(self, job_id):
        pass
