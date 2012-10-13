import os
import time
from itertools import ifilter, islice

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

class fasta_sequence:
    """
        fasta sequence with a header

        two members: header and seq
    """
    __slots__ = ['seq', 'header']
    def __init__(self, h, s):
        self.header=h
        self.seq=s

def fasta_read(input):
    """
        fasta_read(input):

        @param input can be either a file or the name of a file.

        returns a list of fasta_sequence objects with all the sequences in the file.
        comments (lines starting with ';') are ignored.
    """
    if type(input) == str:
        if input.endswith('.gz'):
            import gzip
            input=gzip.gzipfile(input)
        else:
            input=file(input)
    results = []
    header = ''
    seq_items = []
    first = True
    for line in input:
        if line[0] == ';':
            continue # comment
        elif line[0] == '>':
            if not first:
                seq= "".join(seq_items)
                results.append(fasta_sequence(header,seq))
                seq_items = []
            header = line[1:-1] # eat '>' and '\n'
            first = False
        else:
            seq_items.append(line[:-1])
    if len(seq_items) > 0:
        seq = "".join(seq_items)
        results.append(fasta_sequence(header,seq))
    return results
