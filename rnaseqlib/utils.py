##
## Utilities
##
import os
import sys
import time
import datetime
from time import gmtime, strftime
import glob
import re

import operator

from os.path import basename
from urlparse import urlsplit

import itertools
import logging


def invert_dict(d):
    return dict((v,k) for k, v in d.iteritems())


def parse_attributes(attributes_str):
    """
    Parse semi-colon separated key=val attributes
    (e.g. from GFF attributes fields).
    Returns a dictionary.
    """
    if attributes_str.endswith(";"):
        attributes_str = attributes_str[0:-1]
    attribute_fields = attributes_str.split(";")
    attributes = {}
    for attr in attribute_fields:
        if "=" not in attr: continue
        attr_val = attr.split("=")
        attributes[attr_val[0]] = attr_val[1]
    return attributes


def elts_from_list(l, indices):
    """
    Return elements from list by their indices.
    """
    return [l[ind] for ind in indices]


def sort_list_by_list(list_to_sort, sort_values,
                      reverse=False):
    """
    Sort list 'list_to_sort' by the values in
    'sort_values'
    """
    assert (len(list_to_sort) == len(sort_values)), \
           "Lists must be of same length."
    sorted_list = \
        [x for (y,x) in sorted(zip(sort_values, list_to_sort),
                               reverse=reverse)]
    return sorted_list


def get_curr_datetime():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def get_logger(logger_name, log_outdir,
               level=logging.INFO,
               include_stdout=True):
    """
    Return a logging object.
    """
    make_dir(log_outdir)
    logger = logging.getLogger(logger_name)
    formatter = \
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                          datefmt='%m/%d/%Y %I:%M:%S %p')
    log_filename = os.path.join(log_outdir, "%s.log" %(logger_name))
    fh = logging.FileHandler(log_filename)
    fh.setLevel(level)
    fh.setFormatter(formatter)    
    logger.addHandler(fh)
    logging.root.setLevel(level)
    # Optionally add handler that streams all logs
    # to stdout
    if include_stdout:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    logger.info("Created logger %s" %(logger_name))
    return logger


def trim_fastq_ext(fastq_filename):
    """
    Trim .fastq or .fastq.gz (case-insensitive)
    from basename of given fastq filename.
    """
    fastq_dirname = os.path.dirname(fastq_filename)
    fastq_basename = os.path.basename(fastq_filename)
    # Trim trailing .fastq (optionally followed by .gz)
    ext_pattern = re.compile("\.fastq(\.gz)?$", re.IGNORECASE)
    trimmed_basename = ext_pattern.sub("", fastq_basename)
    trimmed_fastq_fname = \
        os.path.join(fastq_dirname, trimmed_basename)
    return trimmed_fastq_fname


def trim_gff_ext(gff_filename):
    """
    Trim .gff or .gff3 (case-insensitive)
    """
    gff_dirname = os.path.dirname(gff_filename)
    gff_basename = os.path.basename(gff_filename)
    ext_pattern = re.compile("\.gff(3)?$", re.IGNORECASE)
    trimmed_basename = ext_pattern.sub("", gff_basename)
    trimmed_gff_fname = \
        os.path.join(gff_dirname, trimmed_basename)
    return trimmed_gff_fname


def endsin_gff(gff_filename):
    """
    Return true if filename ends in gff/gff3 (case-insensitive).
    """
    gff_basename = os.path.basename(gff_filename)
    ext_pattern = re.compile("\.gff(3)?$", re.IGNORECASE)
    result = ext_pattern.search(gff_filename)
    if result is not None:
        return True
    return False


def make_dir(dirpath):
    if os.path.isfile(dirpath):
        print "Error: %s is a file!" %(dirpath)
        sys.exit(1)
    # Try to make the directory
    try:
        os.makedirs(dirpath)
    except OSError:
        pass

    
def gunzip_file(filename, output_dir,
                ext=".txt",
                force=True,
                unless_exists=True):
    """
    Unzip the file into the given output directory.

    If file ends with .gz, strip the .gz. Otherwise,
    add a .txt at the end.
    """
    print "Unzipping: %s into directory %s" %(filename,
                                              output_dir)
    if filename.endswith(".gz"):
        unzipped_filename = filename[0:-3]
    else:
        unzipped_filename = "%s.%s" %(filename,
                                      ext)
    print "  - Unzipped filename: %s" %(unzipped_filename)
    if unless_exists and os.path.isfile(unzipped_filename):
        print "  - File exists, skipping.."
        return unzipped_filename
    os.chdir(output_dir)
    gunzip = "gunzip "
    if force:
        gunzip += "--force"
    os.system("%s %s > %s" %(gunzip,
                             filename,
                             unzipped_filename))
    return unzipped_filename
    
    
def pathify(filename):
    return os.path.abspath(os.path.expanduser(filename))


def url2name(url):
    return basename(urlsplit(url)[2])


def which(program):
    """
    Check if program exists on path.
    """
    def is_exe(fpath):
        if not os.path.isfile(fpath):
            return False
        elif not os.access(fpath, os.X_OK):
            # If the file exists but is not executable, warn
            # the user
            print "WARNING: Found %s but it is not executable." %(fpath)
            print "Please ensure %s is executable." %(fpath)
            print "On Unix, use something like: "
            print "  chmod +x %s" %(fpath)
            time.sleep(10)
            return False
        return True
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
            

def count_lines(fname, skipstart="#"):
    """
    Return number of lines in file.

    Skip liens starting with 'skipstart'.
    """
    num_lines = 0
    for line in open(fname):
        if line.startswith(skipstart): continue
        num_lines += 1
    return num_lines

##
## Misc. utilities for manipulating strings/lists
##
def iter_by_pair(l, step=2):
    """
    Return an iterator over pairs in the list 'l'.

    'step' determines the number of elements to advance
    the pair iterator by.
    """
    return itertools.izip(l, l[1::step])


def unique_list(l):
    """
    Return unique list, with consistent ordering.
    """
    seen = {}
    uniq_l = []
    for elt in l:
        if elt in seen: continue
        uniq_l.append(elt)
        seen[elt] = True
    return uniq_l


def min_item(values):
    """
    Return minimum value index, minimum value
    """
    min_index, min_value = min(enumerate(values),
                               key=operator.itemgetter(1))
    return min_index, min_value


def max_item(values):
    """
    Return maximum value index, maximum value
    """
    max_index, max_value = max(enumerate(values),
                               key=operator.itemgetter(1))
    return max_index, max_value

def min_item(values):
    """
    Return min value index, min value
    """
    min_index, min_value = min(enumerate(values),
                               key=operator.itemgetter(1))
    return min_index, min_value



def get_pairwise_from_sets(first_samples, second_samples):
    seen_pairs = []
    for sample_pair in itertools.product(first_samples,
                                         second_samples):
        sample1, sample2 = sample_pair
        if (sample_pair in seen_pairs) or (sample1 == sample2) or \
            ((sample2, sample1) in seen_pairs):
            continue
        seen_pairs.append(sample_pair)
    return seen_pairs
    
        
def get_pairwise_comparisons(samples):
    """
    Return pairwise comparisons between samples.
    """
    return get_pairwise_from_sets(samples, samples)


def abs_fold_change(fold_changes):
    """
    Get the absolute fold change.
    """
    abs_fc = map(lambda v: \
                 1/float(v) if v <= 1 else v,
                 fold_changes)
    return abs_fc


def flatten(l):
    """
    Flatten a list.
    """
    return [item for sublist in l for item in sublist]


def slices_of_list(l, slice_size):
    list_of_slices = zip(*(iter(l),) * slice_size)
    return list_of_slices


def get_sliding_pairs(l):
    """
    Get sliding pairs from list.

    [1, 2, 3, 4] => [(1,2), (3,4)]

    [1, 2, 3] => [(1,2), (2,3)]
    """
    if len(l) % 2 == 0:
        # Even number of elements lists can be sliced
        # by 2
        pairs = slices_of_list(l, 2)
    else:
        # If it's an odd number of elements list, slide by two
        # up until last element and then append last pair
        pairs = slices_of_list(l[0:-1], 2)
        pairs.append((l[-2], l[-1]))
    return pairs
        

def chunk_list(l, delim):
    record = []
    for r in l:
        if delim in r:
            if record:
                yield record
            record = [r]
        else:
            record.append(r)
    if record:
        yield record         
    

def pathify(filename):
    return os.path.abspath(os.path.expanduser(filename))


def parse_coords(coords):
    """
    Parse coordinates of form:

      chrom:start:end:strand

    Ensure that start < end.

    where ':strand' is optional.
    """
    fields = coords.split(":")
    if len(fields) < 2:
        return None
    strand = None
    chrom = fields[0]
    if len(fields) == 3:
        strand = fields[2]
    start, end = int(fields[1]), int(fields[2])
    if start > end:
        # Reverse coordinates of start > end
        start, end = end, start
    return chrom, start, end, strand


def get_gff_filenames_in_dir(dirname):
    """
    Get all file names that end in the given extensions
    in case-insensitive manner.
    """
    dirname = pathify(dirname)
    #ext_pattern = re.compile("\.gff(3)?$", re.IGNORECASE)
    filenames = glob.glob(os.path.join(dirname, "*.*"))
    matching_fnames = \
        [f for f in filenames \
         if f.lower().endswith((".gff", ".gff3"))]
    return matching_fnames


def in_any(s, l):
    """
    Return True if the string 's' is
    in any element of 'l'.
    """
    return any([s in elt for elt in l])


def is_rpy2_available():
    """
    Return True if rpy2 is available, False otherwise.
    """
    try:
        import rpy2
        return True
    except:
        return False


def get_strtime():
    return strftime("%Y-%m-%d %H:%M:%S", gmtime())




    



