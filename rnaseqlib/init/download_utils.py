##
## Download utilities
##

import os
import sys
import time

import shutil
import urllib2
import posixpath

import rnaseqlib


def wget(url):
    """
    wget a URL.
    """
    t1 = time.time()
    wget_cmd = "wget \'%s\'" %(url)
    os.system(wget_cmd)
    t2 = time.time()
    print "  Downloading took %.2f minutes." %((t2 - t1)/60.)


def download_url(url, output_dir,
                 binary=True,
                 basename=None,
                 unless_exists=True):
    """
    Download a url and put it at the desired location.
    """
    print "Downloading: %s" %(url)
    t1 = time.time()
    url_in = None
    try:
        url_in = urllib2.urlopen(url)
    except urllib2.HTTPError:
        print "Error in urllib2: Could not fetch %s" %(url)
        sys.exit(1)
    url_name = posixpath.basename(url)
    if basename != None:
        url_name = basename
    output_filename = os.path.join(output_dir, url_name)
    if unless_exists and os.path.isfile(output_filename):
        print "  - File exists, skipping."
        return output_filename
    if binary:
        url_out = open(output_filename, "wb")
    else:
        url_out = open(output_filename, "w")
    shutil.copyfileobj(url_in, url_out)
    t2 = time.time()
    print "  Downloading took %.2f minutes." %((t2 - t1)/60.)
    return output_filename

