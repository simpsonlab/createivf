#!/usr/bin/env python


import os
import sys
from glob import glob
import argparse



def init_args():
    """
    Initialize command line arguments
    """
    description = ''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p', '--path', required=True,
            help='')
    parser.add_argument('-o', '--outdir', required=True,
            help='')
    return parser.parse_args()


def is_fast5(file, pattern='.fast5'):
    """
    Check if a file ends in .fast5
    """
    if os.path.basename(file).endswith(pattern):
        return True
    else:
        return False


def get_full_outdir(path):
    """
    Get the absolute path of the output directory to link all FAST5
    files to (resolves all symlinks).
    """
    return os.path.realpath(os.path.abspath(path))


def main():
    """
    Main program
    """
    args = init_args()
    files = glob(f'{args.path}/*/*.fast5')
    if not os.path.isdir(args.outdir):
        print(f'Creating directory {args.outdir} ...')
        os.mkdir(args.outdir)
    outdir = get_full_outdir(args.outdir)
    print(f'Linking files in {args.path} to {outdir}...')
    for file in files:
        if is_fast5(file):
            real_file = os.path.realpath(os.path.abspath(file))
            base_file = os.path.basename(file)
            outlink = f'{outdir}/{base_file}'
            if not os.path.exists(outlink):
                os.symlink(real_file, outlink)
            else:
                continue
        else:
            continue





if __name__ == '__main__':
    main()


#__END__
