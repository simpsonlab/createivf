#!/usr/bin/env python

#
# Code extracted from https://github.com/mike-molnar/nanopore-workflow
# by Dr. Mike Molnar
#

import sys
import os
import re
import csv
import argparse
import pybedtools
from pybedtools import BedTool
from pybedtools import genome_registry
from pybedtools.featurefuncs import extend_fields

def init_args():
    """
    Initialize command line arguments.
    """
    description = 'Filter copy number calls from BED files.'
    # parse the arguements from the user
    parser = argparse.ArgumentParser(description=description)
    # define input BED files
    parser.add_argument('-depth', '--coverage-depth', type=str, required=False,
            help='samtools depth output')
    parser.add_argument('-gaps', '--genome-gaps', type=str, required=True,
            help='full path to the gaps BED file')
    parser.add_argument('-cent', '--centromeres', type=str, required=True,
            help='full path to the centromeres BED file')
    # define output BED files
    parser.add_argument('-cnv_out', '--copy-number-variant-output', type=str, required=False,
            help='name of the copy number output file to write data to')
    # define variables for filtering
    parser.add_argument('-cov', '--mean-coverage', type=float, required=True,
            help='mean coverage of the sample')
    parser.add_argument('-p', '--ploidy', type=int, required=False, default=2,
            help='ploidy of the sample (default: 2)')
    parser.add_argument('-hc', '--high-copy-threshold', type=float, required=False, default=2.7,
            help='upper copy threshold (default: 2.7)')
    parser.add_argument('-lc', '--low-copy-threshold', type=float, required=False, default=1.3,
            help='lower copy threshold (default: 1.3)')
    parser.add_argument('-slop', '--slop-pct', type=float, required=False, default=0.5,
            help='proportion of size to increase a feature by (default: 0.5)')
    parser.add_argument('-size', '--genome-size', type=str, required=True,
            help='genome size file')
    return parser.parse_args()


# define a function to add the copy number to the end of the BED file
def process_cnv(feature):
    # get the copy number information
    coverage = float(feature[3])
    copy_number = (coverage / args.mean_coverage) * args.ploidy
    length = int(feature[2]) - int(feature[1])
    
    # check for gains and/or losses
    gain = 0
    loss = 0
    cnv_list = feature[4].split(",")
    for current_cnv in cnv_list:
        if float(current_cnv) >= args.high_copy_threshold:
            gain += 1
        else:
            loss += 1
    if gain > 0 and loss == 0:
        gain_loss = "gain"
    if gain == 0 and loss > 0:
        gain_loss = "loss"
    if gain > 0 and loss > 0:
        gain_loss = "gain and loss"
    
    # the feature has to be expanded before adding the copy number
    feature = extend_fields(feature, 6)
    feature[3] = str(length)
    feature[4] = str(copy_number)
    feature[5] = str(gain_loss)
    
    return feature

args = init_args()
# store BED files into a BedTool object
if args.coverage_depth:
    in_depth = BedTool(args.coverage_depth)
if os.path.exists(args.genome_gaps):
    in_gaps = BedTool(args.genome_gaps)
else:
    print(f'Gap file {args.genome_gaps} does not exist, exiting...\n')
    sys.exit(1)
if os.path.exists(args.centromeres):
    in_centromeres = BedTool(args.centromeres)
else:
    print(f'Centromere file {args.centromeres} does not exist, exiting...\n')
    sys.exit(1)


# create a BED file of filtered copy number variants
if args.coverage_depth and args.copy_number_variant_output:
    filter_cnv = in_depth\
    .filter(lambda x: ((float(x[3])/args.mean_coverage)*args.ploidy) <= args.low_copy_threshold or ((float(x[3])/args.mean_coverage)*args.ploidy) >= args.high_copy_threshold)\
    .intersect(in_gaps.slop(b=args.slop_pct, pct=True, g=args.genome_size), v=True)\
    .merge(d=1000000, c=(4,4), o=("mean","collapse"))
    #.intersect(in_gaps.slop(b=args.slop_pct, pct=True, genome="hg38"), v=True)\

    # process the estimated copy number variants
    cnv_final = filter_cnv.each(process_cnv)
    cnv_final.saveas(args.copy_number_variant_output)

