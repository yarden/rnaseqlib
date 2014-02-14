##
## Utilities for dealing with the atrocious RES format
## http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/res
##
import os
import time
import pandas
import re

def output_res_as_tsv(res_fname, output_fname,
                      delim_header="\t+",
                      delim_row="\t\s?",
                      delim_out="\t",
                      na_val="NA"):
    """
    Output RES as TSV file.

    - delim_in: Delimiter in RES format
    - delim_out: Delimiter to use in conversion
    """
    res_in = open(res_fname)
    header = res_in.readline().strip()
    header_fields = re.split(delim_header, header)
    # Ignore description header field
    # and treat first column as Gene
    header_fields = ["Gene"] + header_fields[1:]
    print "Converting RES to TSV..."
    print "  - Input file: %s" %(res_fname)
    print "  - Output file: %s" %(output_fname)
    if os.path.isfile(output_fname):
        print "Already found file %s" %(output_fname)
        print "Skipping..."
        return
    t1 = time.time()
    line_num = 0
    with open(output_fname, "w") as output_file:
        # Write new header
        output_file.write(delim_out.join(header_fields) + "\n")
        for line in res_in:
            line_num += 1
            line = line.rstrip()
            curr_fields = re.split(delim_row, line)
            # Skip lines with one or fewer fields
            if len(curr_fields) <= 1:
                print "Skipping line number %d" %(line_num)
                continue
            # Get every other element (ignore A/P/M calls
            # of microarray)
            filtered_fields = curr_fields[0:2] + curr_fields[2::2]
            if filtered_fields[0] == "":
                # Record missing gene names as NA
                filtered_fields[0] = na_val
            n = 0
            for x in filtered_fields:
                if x == "A":
                    print filtered_fields
                    raise Exception, "%s" %(line)
                n += 1
            tsv_line = delim_out.join(filtered_fields)
            output_file.write(tsv_line + "\n")
    t2 = time.time()
    print "Conversion took %.2f seconds" %(t2 - t1)
        
        
