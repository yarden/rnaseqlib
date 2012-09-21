##
## Utilities
##
import os

def make_dir(dirpath):
    if os.path.isfile(dirpath):
        print "Error: %s is a file!" %(dirpath)
        sys.exit(1)
    # Try to make the directory
    try:
        os.makedirs(dirpath)
    except OSError:
        pass
    
def pathify(filename):
    return os.path.abspath(os.path.expanduser(filename))
