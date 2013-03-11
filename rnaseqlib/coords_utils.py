##
## Dealing with genomic coordinates. Try to stick to
## bedtools-style coordinates and conventions internally
## whenever possible.
##
import pybedtools


def overlap(first_coords, second_coords):
    """
    Return overlap between two intervals.
    """
    overlap = pybedtools.overlap(first_coords[0], second_coords[0],
                                 first_coords[1], second_coords[1])
    return overlap
