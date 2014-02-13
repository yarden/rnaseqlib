##
## Utilities for dealing with the atrocious RES format
## http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/res
##
import os
import time
import pandas
import re

def output_res_as_tsv(res_fname, output_fname,
                      delim_in="\s+",
                      delim_out="\t"):
    """
    Output RES as TSV file.

    - delim_in: Delimiter in RES format
    - delim_out: Delimiter to use in conversion
    """
    res_in = open(res_fname)
    header = res_in.readline().strip()
    header_fields = re.split(delim_in, header)
    print "Converting RES to TSV..."
    print "  - Input file: %s" %(res_fname)
    print "  - Output file: %s" %(output_fname)
    t1 = time.time()
    with open(output_fname, "w") as output_file:
        # Write new header
        output_file.write(delim_out.join(header_fields) + "\n")
        for line in res_in:
            line = line.strip()
            fields = re.split(delim_in, line)
            # Skip lines with one or fewer fields
            if len(fields) <= 1:
                continue
            curr_fields = re.split(delim_in, line)
            tsv_line = "\t".join(curr_fields) + "\n"
            output_file.write(tsv_line)
    t2 = time.time()
    print "Conversion took %.2f seconds" %(t2 - t1)
        
        
