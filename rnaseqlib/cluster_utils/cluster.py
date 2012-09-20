##
## Utilities for running on cluster
##
import os
import sys
import time


def launch_job(cmd, job_name, settings):
    """
    Launch job on cluster and return a job id.

    Wrapper to Mysge/Mypbm/Mybsub
    """
    # Have a bsub/qsub branch here
    # LAUNCH JOB HERE..
    job_id = None
    return job_id
