#! /usr/bin/env python

import pysam
threshold = 1000

infile = pysam.AlignmentFile("-", "r")
for alignment in infile:
    if alignment.reference_length > threshold and alignment.mapping_quality >= 1:
        print ("%s\t%s\t%d\t%d" % (alignment.query_name, alignment.reference_name, alignment.reference_start, alignment.reference_end))
