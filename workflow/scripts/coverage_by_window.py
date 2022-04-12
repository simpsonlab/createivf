import os
import sys
import csv
import argparse
import re

def init_args():
    """
    Initialize command line arguments
    """
    parser = argparse.ArgumentParser( description='Calculate mean coverage levels in a sorted BED file by window size.')
    parser.add_argument('-i', '--input', type=str, required=True,
            help='the sample coverage file from "samtools depth"')
    parser.add_argument('-o', '--output', type=str, required=True,
            help='name of output file to write data to')
    parser.add_argument('-w', '--window', type=int, required=False, default=500000,
            help='window size to calculate segment coverage (default: 500000)')
    #parser.add_argument('-c', '--coverage', type=float, required=True)
    parser.add_argument('-p', '--ploidy', type=float, required=False, default=2.0,
            help='ploidy (default: 2.0)')
    parser.add_argument('-n', '--nanostats', type=str, required=True,
            help='full path to the NanoStats.tsv file')
    return parser.parse_args()


def get_coverage(file, pattern='Total bases aligned', genome_size=3100000000):
    """
    Get the average coverage from the NanoStats.txt file
    """
    total_bases = 0.0
    if os.path.exists(file):
        with open(file, 'r') as ifh:
            for line in ifh:
                if re.search(pattern, line):
                    total_bases = float(line.split(':')[1].strip().replace(',', ''))
    else:
        print(f'File {file} does not exist, exiting...\n')
        sys.exit(1)
    return float(total_bases)/genome_size


def main():
    """
    Main program
    """
    args = init_args()
    in_fh = open(args.input, 'r')
    out_fh = open(args.output, 'w')
    
    csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'locus', 'depth'], delimiter='\t')
    
    region_start = 0
    region_end = args.window
    total_coverage = 0
    total_locus = 0
    
    for record in csv_reader:
        if int(record['locus']) > region_end:
            while int(record['locus']) > region_end:
                if total_locus == 0:
                    out_fh.write("%s\t%d\t%d\t%d\t%.2f\n" % (record['chromosome'], region_start, region_end, 0, 0.0))
                else:
                    mean_coverage = total_coverage/total_locus
                    normalized = float(mean_coverage/(get_coverage(file=args.nanostats) * args.ploidy))
                    if normalized > 1.0:
                        normalized = 1.0
                    out_fh.write("%s\t%d\t%d\t%d\t%.2f\n" % (record['chromosome'], region_start, region_end, mean_coverage, normalized))
                total_coverage=0
                total_locus = 0
                region_start = region_end + 1
                region_end = region_end + args.window
            total_coverage = int(record['depth'])
            total_locus = 1
        else:
            total_coverage = total_coverage + int(record['depth'])
            total_locus = total_locus + 1


if __name__ == '__main__':
    main()


#__END__
