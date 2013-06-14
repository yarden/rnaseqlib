##
## Utilities for FASTQ
##

import os
import time
from itertools import ifilter, islice

import gzip

def read_open_fastq(fastq_filename):
    fastq_file = None
    if fastq_filename.endswith(".gz"):
        fastq_file = gzip.open(fastq_filename, "rb")
    else:
        fastq_file = open(fastq_filename, "r")
    return fastq_file

def get_fastq_entries(fastq_filename):
    fastq_handle = read_open_fastq(fastq_filename)
    fastq_entries = read_fastq(fastq_handle)
    return fastq_entries

    
def write_open_fastq(fastq_filename):
    fastq_file = None
    if fastq_filename.endswith(".gz"):
        fastq_file = gzip.open(fastq_filename, "wb")
    else:
        fastq_file = open(fastq_filename, "w")
    return fastq_file
    

def read_fastq(fastq_in):
    """
    parse a fastq-formatted file, yielding a
    (header, sequence, header2, quality) tuple
    """
    fastqfile = None
    if type(fastq_in) == str:
        # Parse as filename
        fastqfile = read_open_fastq(fastq_in)
    else:
        # Otherwise assume it's a file handle
        fastqfile = fastq_in
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines
    fastqiter = ifilter(lambda l: l, fastqiter)  # skip blank lines
    line_num = 0
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1,seq,header2,qual = fqlines
            # Allow header1 to have '@' somewhere in it, not just in the
            # first line, mainly to header FASTQ files with odd headers
            if (header1.startswith('@') or ("@" in header1)) and \
                header2.startswith('+'):
                yield header1[1:], seq, header2, qual
            else:
                print "Problem with formatting of FASTQ file detected."
                print "header1: ", header1, " header2: ", header2
                raise ValueError("Invalid header lines in FASTQ: %s and %s (line %d)" \
                                 %(header1, header2, line_num))
        elif len(fqlines) == 0:
            raise StopIteration
        else:
            raise EOFError("Failed to parse four lines from fastq file!")
        line_num += 1


def write_fastq(fastq_file, fastq_rec):
    header, seq, header2, quality = fastq_rec
    if not header.startswith("@"):
        header = "@%s" %(header)
    fastq_file.write("%s\n" %(header))
    fastq_file.write("%s\n" %(seq))
    fastq_file.write("%s\n" %(header2))
    fastq_file.write("%s\n" %(quality))
    
                                                                                
def fastq_get_rec_id(line):
    if line.startswith("@"):
        return line[1:]
    else:
        return line


def fastq2fasta(settings_filename,
                output_dir,
                fieldname="fastq_filenames"):
    """
    Convert FASTQ to FASTA.
    """
    output_dir = os.path.join(output_dir, "fasta")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    settings_info, parsed_settings = \
                   settings.load_settings(settings_filename)

    fastq_filenames = settings_info["data"][fieldname]

    queue_name = "quick"

    for fastq in fastq_filenames:
        sample_id, fastq_filename = fastq
        print "Converting %s to FASTA" %(sample_id)
        print "  FASTQ: %s" %(fastq_filename)
        job_id = "fastq2fasta_%s" %(sample_id)
        fastq_basename = os.path.basename(fastq_filename)
        output_filename = os.path.join(output_dir,
                                       "%s.fa" %(fastq_basename))

        if os.path.isfile(output_filename):
            print "WARNING: %s exists. Aborting..." \
                  %(output_filename)
            continue
        fastq2fasta_cmd = "perl %s fq2fa %s > %s" %(FQCONVERT,
                                                    fastq_filename,
                                                    output_filename)
        print "Running: %s" %(fastq2fasta_cmd)
        cluster.run_on_cluster(fastq2fasta_cmd, job_id, output_dir,
                               queue_name=queue_name,
                               other_args="-l nodes=1:E5620")

def file_type(input):
    """given an input file, determine the type and return both type and record delimiter (> or @)"""
    name, extension = os.path.splitext(os.path.basename(input))
    fastas = set(['.fsa','.fasta','.fa'])
    fastqs = set(['.fastq','.fq'])
    gffs = set(['.gff'])
    if extension in fastas:
        ft = 'fasta'
        delim = '>'
    elif extension in fastqs:
        ft = 'fastq'
        delim = '@'
    # TODO:  sff ???
    #elif extension in gffs:
    #    ft = 'sff'
    #    delim = None
    else:
        raise IOError, "Input file not of correct extension"
    return ft, delim
    
def _get_file_chunks(input, delim, size):
    """given input, record delimiter, and chunk size, yield an iterator contains file seek (start)
    and file read (stop) positions.  Return final position as (6365605, None)."""
    f = open(input)
    while 1:
        start = f.tell()
        f.seek(size, 1)
        line = f.readline()
        if not line:
            break
        # if this isn't a fasta header line, read forward until
        # we get to one
        while not line.startswith(delim):
            line = f.readline()
        else:
            # now that we got to a fasta header, we're at the end.
            # back up the length of the fasta header.
            f.seek(-len(line), 1)
            # tuple up
            yield start, f.tell() - start, input
    # make sure we catch the (start, distance) for the end of the file, too
    yield start, None, input
    f.close()
    
def get_chunks(input, delim, split_type, mb=1, splits = None):
    """return a tuple of file seek (start, distance) positions covering chunks of a file"""
    if split_type == 'size':
        size = mb * (1024**2)
#    if split_type == 'pieces':
#        if not splits:
#            splits = multiprocessing.cpu_count() - 1
#        size = int(round((os.path.getsize(input)/float(splits)), 0))
    return _get_file_chunks(input, delim, size)


def split_file(chunk, output_filename):
    """function to split a file into appropriate pieces given (start, stop) file seek coords"""
    f = open(chunk[2])
    f.seek(chunk[0])
    if chunk[1]:
        d = f.read(chunk[1])
    else:
        d = f.read()
    otf = open(output_filename, 'w')
    otf.write(d)
    otf.close()
    f.close()


def chunk_fasta(fasta_filename, output_dir,
                mb=1):
    """
    Chunk FASTA files.
    """
#    output_dir = os.path.join(output_dir)
#                              "chunked_%d_mb" %(mb))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    print "Chunking %s..." %(fasta_filename)

    t1 = time.time()

    ft, delim = file_type(fasta_filename)

    fasta_chunks = get_chunks(fasta_filename, delim,
                              "size", mb=mb)
    
    chunk_num = 1
    fasta_basename = os.path.basename(fasta_filename)

    if fasta_basename.endswith(".fa"):
        fasta_basename = fasta_basename[0:-3]
    elif fasta_basename.endswith(".fasta"):
        fasta_basename = fasta_basename[0:-6]
    
    for chunk in fasta_chunks:
        chunk_filename = os.path.join(output_dir,
                                      "%s.chunk_%d.fasta" \
                                      %(fasta_basename,
                                        chunk_num))
        split_file(chunk, chunk_filename)
        chunk_num += 1

    t2 = time.time()
    print "Chunking into %d chunks took %.2f seconds" \
          %(chunk_num - 1, (t2 - t1))


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--fastq2fasta", dest="fastq2fasta", nargs=1,
                      help="Convert FASTQ to FASTA. Takes settings file with "
                      "FASTQ filenames.")
    parser.add_option("--fastq-fieldname", dest="fastq_fieldname", default="fastq_files",
                      type="str", nargs=1)
    parser.add_option("--chunk-fasta", dest="chunk_fasta", nargs=2, default=None,
                      help="Chunk FASTA filename. Takes FASTA file and size (in Megabytes) to "
                      "chunk to.")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None,
                      help="Output directory.")
    (options, args) = parser.parse_args()

    if options.output_dir == None:
        print "Error: need --output-dir"
        return

    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    if options.fastq2fasta != None:
        settings_filename = os.path.abspath(os.path.expanduser(options.fastq2fasta))
        fastq2fasta(settings_filename, output_dir, fieldname=options.fastq_fieldname)
    elif options.chunk_fasta != None:
        fasta_filename = os.path.abspath(os.path.expanduser(options.chunk_fasta[0]))
        chunk_size = int(options.chunk_fasta[1])
        chunk_fasta(fasta_filename, output_dir,
                    mb=chunk_size)

if __name__ == '__main__':
    main()
