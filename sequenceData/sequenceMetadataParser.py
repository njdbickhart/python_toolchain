# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 08:41:05 2022

@author: Derek.Bickhart
"""

import argparse
import gzip
from collections import defaultdict
import contextlib

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A wrapper script for the FlashFry CRISPR design tool"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input reference fasta file. If a database has not been created before, the wrapper will make one",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are (basename).output and (basename).output.scored",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser

def main(args):
    print("hey")
    

def determineFile(filename):
    token = ''
    fsegs = filename.rstrip().split('.')
    if filename.endswith('gz'):
        token = fsegs[-2]
    else:
        token = fsegs[-1]
        
    if token == 'bam':
        return 'bam'
    elif token == 'vcf':
        return 'vcf'
    elif token == ''

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'r'):
    if filename.endswith('gz'):
        fh = gzip.open(filename, 'rt', encoding='utf-8')
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        fh.close()

if __name__ == "__main__":
    args = arg_parse()
    main(args)