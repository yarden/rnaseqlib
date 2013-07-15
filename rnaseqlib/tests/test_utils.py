import os
import sys
import time

TESTDIR = os.path.dirname(os.path.abspath(__file__))

def load_test_data(name):
    """
    Return filename of test data.
    """
    test_fname = os.path.join(TESTDIR, "test_data", name)
    if not os.path.isfile(test_fname):
        raise Exception, "Cannot find %s" %(test_fname)
    return test_fname
    
    
    
    
