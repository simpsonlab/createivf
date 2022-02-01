#!/usr/bin/env python

"""
Python script to generate read support files for breakpoints.
"""

import os
import sys
import argparse
import csv
import pysam
from operator import itemgetter

def init_args():
    """
    Initialize command line arguments
    """
    description = 'CReATe IVF breakpoint generator'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--sample', required=True,
            help='name of the sample processed')
    parser.add_argument('-m', '--metadata', required=True,
            help='path to the sample metadata file')
    parser.add_argument('-c', '--cytobands', required=True,
            help='path to the cytoband composite file')
    parser.add_argument('-b', '--bam', required=True,
            help='full path to the merged sorted BAM file to process')
    parser.add_argument('-o', '--outprefix', default='bp',
            help='read file prefix (default: bp)')
    parser.add_argument('-t', '--threshold', required=False, default=1000, type=int,
            help='reference length threshold (default: 1000)')
    parser.add_argument('--mapping_quality', required=False, default=1, type=int,
            help='mapping quality lower threshold (default: 1)')
    return parser.parse_args()


def add_chr_prefix(band):
    """
    Return the band string with chr prefixed
    """
    return ''.join(['chr', band])


def import_metadata(file, delimiter='\t'):
    """
    Import the metadata file and return a
    dictionary containing the contents.
    """
    target = dict()
    with open(file, 'r') as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for line in reader:
            if line['sample'] not in target:
                target.update({line['sample'] :
                    {'region1': add_chr_prefix(band=line['region1']),
                     'region2': add_chr_prefix(band=line['region2'])}})
            else:
                print(f'{line["sample"]} found in dictionary, skipping...')
                continue
            #(sampleid, region1, region2) = line.rstrip().split()
    return target


def get_region_for_sample(sample, bpid, target):
    """
    Parse a dictionary and obtain the region for a given sample
    """
    for sample_id in target:
        if sample.find(sample_id) != -1:
            if bpid == 1:
                return target[sample_id]['region1']
            else:
                return target[sample_id]['region2']
            break


def import_cytoband(file, target, delimiter='\t'):
    """
    Returns the genomic region from the cytoband composite file
    for a given target
    """
    bp_chr = list()
    bp_coords = list()
    with open(file, 'r') as fh:
        for line in fh:
            fields = line.rstrip().split()
            if len(fields) != 6:
                continue
            (chrom, start, end, _, _, key) = fields
            start = int(start)
            end = int(end)

            # handle ambiguous regions
            if key.find(target) != -1:
                bp_chr.append(chrom)
                bp_coords.append(start)
                bp_coords.append(end)
    bp_coords = sorted(bp_coords)
    chrom = bp_chr[0]
    start = bp_coords[0]
    end = bp_coords[-1]
    region = {"chrom" : chrom, "start": start, "end": end}
    return region
    #return f'{chrom}:{start}-{end}'


def get_reads_in_region(bam, chrom, start, end, threshold=1000, quality=1):
    """
    Print
    """
    region_reads = list()
    bamfile = pysam.AlignmentFile(bam, 'rb')
    for read in bamfile.fetch(chrom, start, end):
        if read.reference_length > threshold and read.mapping_quality >= 1:
            region_reads.append(read)
    return region_reads


def print_region_reads(outfile, reads):
    """
    Print the reads to a file
    """
    with open(outfile, 'w') as ofh:
        for read in reads:
            ofh.write('\t'.join(read))
            ofh.write('\n')
    ofh.close()


def create_sorted_reads(reads):
    """
    Create a list of sorted reads by query name (read ID)
    """
    readids = list()
    read_list = list()
    for read in reads:
        read_list.append([read.query_name,
                          read.reference_name,
                          str(read.reference_start),
                          str(read.reference_end)])
    return sorted(read_list, key=itemgetter(0))


def get_intersecting_reads(reads1, reads2):
    """
    Return reads that intersect breakpoint 1 and breakpoint
    2 files.
    """
    read1_dict = dict()
    read2_dict = dict()
    intersecting_reads = list()
    for read1 in reads1:
        if read1.query_name not in read1_dict:
            read1_dict.update({read1.query_name:
                {'chrom': read1.reference_name,
                 'start': str(read1.reference_start),
                 'end': str(read1.reference_end)}})
        else:
            continue
    for read2 in reads2:
        if read2.query_name in read1_dict:
            tmp_read1 = ' '.join([
                read1_dict[read2.query_name]['chrom'],
                read1_dict[read2.query_name]['start'],
                read1_dict[read2.query_name]['end']
                ])
            tmp_read2 = ' '.join([
                read2.reference_name,
                str(read2.reference_start),
                str(read2.reference_end)
                ])
            intersecting_reads.append([read2.query_name, tmp_read1, tmp_read2])
    return intersecting_reads


def print_intersecting_reads(file, reads):
    """
    Print reads that intersect breakpoint 1 and 2 to a file
    """
    with open(outfile, 'w') as ofh:
        for read in reads:
            ofh.write('\t'.join(read))
            ofh.write('\n')
    ofh.close()


def main():
    """
    Main program
    """
    args = init_args()
    target_dict = import_metadata(file=args.metadata)
    bp1_cyto = get_region_for_sample(sample=args.sample, bpid=1, target=target_dict)
    bp2_cyto = get_region_for_sample(sample=args.sample, bpid=2, target=target_dict)
    bp1_region = import_cytoband(file=args.cytobands, target=bp1_cyto)
    bp2_region = import_cytoband(file=args.cytobands, target=bp2_cyto)
    print('Processing read 1 region')
    bp1_reads = get_reads_in_region(bam=args.bam,
            chrom=bp1_region['chrom'],
            start=bp1_region['start'],
            end=bp1_region['end'],
            threshold=args.threshold,
            quality=args.mapping_quality)
    bp1_reads_list = create_sorted_reads(reads=bp1_reads)
    print('Processing read 2 region\n')
    bp2_reads = get_reads_in_region(bam=args.bam,
            chrom=bp2_region['chrom'],
            start=bp2_region['start'],
            end=bp2_region['end'],
            threshold=args.threshold,
            quality=args.mapping_quality)
    bp2_reads_list = create_sorted_reads(reads=bp2_reads)

    # write the reads to a tab separate file
    outfile1 = ''.join([args.outprefix, '1.reads'])
    print(f'Writing reads for breakpoint 1 in {outfile1}...')
    print_region_reads(outfile=outfile1, reads=bp1_reads_list)
    outfile2 = ''.join([args.outprefix, '2.reads'])
    print(f'Writing reads for breakpoint 2 in {outfile2}...\n')
    print_region_reads(outfile=outfile2, reads=bp2_reads_list)

    # get list of reads that occur in both breakpoint 1 and breakpoint 2
    int_reads = get_intersecting_reads(reads1=bp1_reads, reads2=bp2_reads)

    print(f'Detecting intersecting breakpoint reads...')
    int_file = ''.join([args.outprefix, '.intersecting.reads'])
    if len(int_reads) > 0:
        print_intersecting_reads(file=int_file, reads=int_reads)
    else:
        print('No overlapping breakpoints detected...\n')
    print('Breakpoint detection complete!\n')
    

if __name__ == '__main__':
    main()


#__END__

