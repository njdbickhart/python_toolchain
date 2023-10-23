#!/usr/bin/python3
# This is a script designed to remove Freebayes gVCF non-variant information for use in standard calls

import os
import sys
import argparse
from Bio import bgzf

version = '0.0.1'

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Modify a gVCF file to create a preformatted VCF file. Version:" + version
            )
    parser.add_argument('-f', '--file',
                        help="A single gVCF file for modification",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--chr',
                        help="Modify a chromosome or scaffold name using comma delimiters (ie. X,39; may be entered more than once)",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--sample',
                        help="Replace the sample name with the following string (default: None)",
                        type=str, default='None'
                        )
    parser.add_argument('-o', '--output',
                        help="Output vcf file name",
                        type=str, required=True
                        )

    return parser.parse_args()

def main(args):
    # Save chromosome modification text
    modchr = False
    conversion = dict()
    if len(args.chr) > 0:
        modchr = True
        for c in args.chr:
            for o, n in c.split(','):
                conversion[o] = n

    # open the VCF and start processing
    fh = smartFile(args.file)
    with open(args.output, 'w') as out:
        for l in fh:
            if l.startswith('##'):
                out.write(l)
            elif l.startswith('#CHROM'):
                if args.sample != 'None':
                    s = l.rstrip().split()
                    s[-1] = args.sample
                    out.write('\t'.join(s) + '\n')
                else:
                    out.write(l)
            else:
                s = l.rstrip().split()
                if s[4] == '<*>':
                    continue
                s[0] = conversion.get(s[0], s[0])
                out.write('\t'.join(s) + '\n')

    fh.close()

def smartFile(filename : str, mode : str = 'r'):
    fh = None
    if filename.endswith('.gz') and mode == 'r':
        fh = bgzf.open(filename, mode='r')
    elif filename.endswith('.gz') and mode == 'w':
        fh = bgzf.open(filename, mode='w')
    else:
        fh = open(filename, mode)
    return fh

if __name__ == "__main__":
    args = parse_user_input()
    main(args)
