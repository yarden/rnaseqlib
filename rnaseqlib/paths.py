##
## Relevant paths
##
import os
import sys

PIPELINE_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_RUN_SCRIPT = os.path.join(PIPELINE_CODE_DIR, "run_pipeline.py")
SCRIPTS_DIR = os.path.join(PIPELINE_CODE_DIR, "..", "scripts")
