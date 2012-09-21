##
## Utilities for running on cluster
##
import os
import sys
import time

class Cluster:
    """
    Cluster submission.
    """
    def __init__(self, settings_info,
                 supported_types=["bsub", "qsub"]):
        self.cluster_type = settings_info["mapping"]["cluster_type"]
        if self.cluster_type not in supported_type:
            print "Error: unsupported cluster type %s" %(self.cluster_type)
            sys.exit(1)
            

    def launch_job(self, cmd, job_name):
        """
        Launch job on cluster and return a job id.

        Wrapper to Mysge/Mypbm/Mybsub
        """
        job_id = None
        if self.cluster_type == "bsub":
            self.launch_job(cmd, job_name, self.cluster_type)
        return job_id
        

    def wait_on_job(self, job_id):
        pass
