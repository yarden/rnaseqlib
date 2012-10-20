##
## Relevant paths
##
import os
import sys

PIPELINE_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.abspath(os.path.join(PIPELINE_CODE_DIR, "..", "scripts"))
PIPELINE_RUN_SCRIPT = os.path.join(SCRIPTS_DIR, "rna_pipeline.py")
