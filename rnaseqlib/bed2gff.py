##
## BED2GFF
##
# import pybedtools
# import sys
# import os

# if len(sys.argv) < 3:
#     raise Exception, "Need arguments: input.bed output.gff"

# if not os.path.isfile(sys.argv[1]):
#     raise Exception, "%s does not exist" %(sys.argv[1])

# if not os.path.isfile(sys.argv[2]):
#     raise Exception, "%s does not exist" %(sys.argv[2])

# input_bed = pybedtools.BedTool(sys.argv[1])
# output_bed = input_bed.saveas(sys.argv[2])
