##
## Utilities for running on cluster
##
import os
import sys
import time

import rnaseqlib
from rnaseqlib.cluster_utils import Mybsub, Mypbm, Mysge

class Cluster:
    """
    Cluster submission.
    """
    def __init__(self,
                 cluster_type,
                 output_dir,
                 logger,
                 supported_types=["bsub", "qsub"]):
        self.logger = logger
        self.cluster_type = cluster_type
        self.output_dir = output_dir
        if self.cluster_type not in supported_types:
            self.logger.critical("Unsupported cluster type: %s" %(self.cluster_type))
            print "Error: unsupported cluster type %s" %(self.cluster_type)
            sys.exit(1)
            

    def launch_and_wait(self, cmd, job_name,
                        unless_exists=None,
                        extra_sleep=60):
        """
        Launch job and wait until it's done.
        """
        job_id = self.launch_job(cmd, job_name,
                                 unless_exists=unless_exists)
        if job_id is None:
            self.logger.critical("Job submission failed for %s, %s" %(cmd,
                                                                      job_name))
            # Job submission failed
            return None
        else:
            # Job is submitted (assigned an ID) so now
            # wait for it to finish
            self.wait_on_job(job_id)
        time.sleep(extra_sleep)
    

    def launch_job(self, cmd, job_name,
                   unless_exists=None):
        """
        Launch job on cluster and return a job id.

        if unless_exists flag is given, do not execute command
        if the given filename path exists.
        
        Wrapper to Mysge/Mypbm/Mybsub.
        """
        ##
        ## TODO: add wrapper for Python multiprocess for local running
        ##
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
                                      self.output_dir,
                                      queue_type="normal")
        if job_id is None:
            print "WARNING: Job %s not submitted." %(job_name)
        return job_id
        

    def wait_on_job(self, job_id):
        if self.cluster_type == "bsub":
            print "Waiting on %s.." %(job_id)
            Mybsub.waitUntilDone(job_id)
            print "  - Completed"
            return True
        else:
            raise Exception, "Not implemented yet."
        

    def wait_on_jobs(self, job_ids):
        num_jobs = len(job_ids)
        print "Starting to wait on a collection of %d jobs" \
            %(num_jobs)
        jobs_completed = {}
        for job_id in job_ids:
            if job_id in jobs_completed: continue
            if self.wait_on_job(job_id):
                jobs_completed[job_id] = True
        print "All jobs completed."
