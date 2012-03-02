##
## Utilities
##
import os

def pathify(filename):
    return os.path.abspath(os.path.expanduser(filename))
