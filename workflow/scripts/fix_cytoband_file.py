#!/usr/bin/env python

"""
Append a column to the UCSC cytoband file
"""

import os
import sys
import argparse
import csv


def init_args():
    """
    Initialize command line arguments
    """
    description = 'Append chromosome identifier column to UCSC cytoband file'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--file', required=True,
        help='full path to UCSC cytoband file')
    return parser.parse_args()


def main():
    """
    Main program
    """
    args = init_args()
    fieldnames = ['chrom', 'start', 'end', 'band', 'desc']
    with open(args.file, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t', fieldnames=fieldnames)
        for line in reader:
            chr_cyto = ''.join(line['chrom'], line['band']])
            print('\t'.join([
                line['chrom'],
                line['start'],
                line['end'],
                line['band'],
                line['desc'],
                chr_cyto
            ]))
    ifh.close()


if __name__ == '__main__':
    main()


#__END__
