#! /usr/bin/env python

import sys
import argparse


def init_args():
    """
    Initialize command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--file1', required=True,
            help='first breakpoint file')
    parser.add_argument('--file2', required=True,
            help='second breakpoint file')
    return parser.parse_args()

args = init_args()
#f1 = sys.argv[1]
f1 = args.file1
#f2 = sys.argv[2]
f2 = args.file2

d = dict()

fh1 = open(f1)
for line in fh1:
    fields = line.rstrip().split()
    d[fields[0]] = fields

fh2 = open(f2)
for line in fh2:
    fields = line.rstrip().split()
    if fields[0] in d:
        print (fields[0], " ".join(d[fields[0]][1:4]), " ".join(fields[1:4]))
