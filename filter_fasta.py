#!/usr/bin/env python3

import argparse
import os
import gzip

def parse_args():
    parser = argparse.ArgumentParser(
    description='Remove contigs with fewer than N bases from fasta file.',
    usage='%(prog)s -i INPUT_FILE -o OUTPUT_FILE [-h] [-n MIN_BASES] [-z]'
    )
    parser.add_argument('-i', '--input_file', help='Input fasta file', required=True)
    parser.add_argument('-o', '--output_file', help='Output fasta file', required=True)
    parser.add_argument('-n', '--min_bases', type=int, default=0, help='Minimum number of bases for contigs to keep (default: 0)')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress output file using gzip')
    return parser.parse_args()

def main():
    """
    Remove contigs with fewer than N bases from fasta file.

    Usage: remove_short_contigs.py [-h] [-n MIN_BASES] [-z] -i INPUT_FILE -o OUTPUT_FILE

    required arguments:
      -i INPUT_FILE, --input_file INPUT_FILE
                            Input fasta file
      -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            Output fasta file

    optional arguments:
      -h, --help            show this help message and exit
      -n MIN_BASES, --min_bases MIN_BASES
                            Minimum number of bases for contigs to keep (default: 0)
      -z, --gzip            Compress output file using gzip
    """
    args = parse_args()

    if args.gzip:
        output_file = gzip.open(args.output_file, 'wt')
    else:
        output_file = open(args.output_file, 'w')

    with gzip.open(args.input_file, 'rt') if args.input_file.endswith('.gz') else open(args.input_file, 'r') as f:
        header, seq = None, []
        for line in f:
            if line.startswith('>'):
                if header:
                    if len(''.join(seq)) >= args.min_bases:
                        output_file.write(header + ''.join(seq) + '\n')
                    header, seq = line, []
                else:
                    header = line
            else:
                seq.append(line.strip())

        if header:
            if len(''.join(seq)) >= args.min_bases:
                output_file.write(header + ''.join(seq) + '\n')

    output_file.close()

if __name__ == '__main__':
    main()
