##
## Utilities
##
import os
import sys
import time

def make_dir(dirpath):
    if os.path.isfile(dirpath):
        print "Error: %s is a file!" %(dirpath)
        sys.exit(1)
    # Try to make the directory
    try:
        os.makedirs(dirpath)
    except OSError:
        pass

    
def gunzip_file(filename, output_dir):
    print "Unzipping: %s into directory" %(filename,
                                           output_dir)
    os.chdir(output_dir)
    os.system("gunzip %s" %(filename))
    
    
def pathify(filename):
    return os.path.abspath(os.path.expanduser(filename))

from os.path import basename
from urlparse import urlsplit

def url2name(url):
    return basename(urlsplit(url)[2])

# def download_url(url, localFileName = None):
#     localName = url2name(url)
#     req = urllib2.Request(url)
#     r = urllib2.urlopen(req)
#     if r.info().has_key('Content-Disposition'):
#         # If the response has Content-Disposition, we take file name from it
#         localName = r.info()['Content-Disposition'].split('filename=')[1]
#         if localName[0] == '"' or localName[0] == "'":
#         localName = localName[1:-1]
#     elif r.url != url: 
#         # if we were redirected, the real file name we take from the final URL
#         localName = url2name(r.url)
#     if localFileName: 
#         # we can force to save the file as specified name
#         localName = localFileName
#     f = open(localName, 'wb')
#     f.write(r.read())
#     f.close()



