# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:00:00 2020

@author: derek.bickhart-adm
"""

import argparse
import re

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Rapidly reformat fasta files to conform to standards"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input fasta file to be reformatted",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output reformatted fasta file",
                        required=True, type=str,
                        )
    parser.add_argument('-l', '--lines',
                        help="Reformatted line number",
                        default=60, type=int,
                        )
    parser.add_argument('-v', '--verbose',
                        help="Verbose STDOUT details",
                        action='store_true'
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    lnum = args.lines
    with open(args.fasta, 'r') as fin, open(args.output, 'w') as fout:
        pattern = re.compile(rf'(.{{{lnum}}})')
        for name, seq in fastaReader(fin):
            fout.write(f'>{name}\n')
            if args.verbose:
                num = len(seq)
                (seq, nsubs) = re.subn(pattern, r'\1\n', seq)
                print(f'For {name} of length {num} implemented {nsubs} new lines')
            else:
                seq = re.sub(pattern, r'\1\n', seq)
            fout.write(f'{seq}\n')
    
def fastaReader(infile):
    header, sequence = '', ''
    for line in infile:
        if line[0] == '>':
            if sequence: yield header, sequence 
            header = line.rstrip().replace('>','').split()[0]
            sequence =''
        else:
            sequence += line.rstrip()
    yield header, sequence

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
