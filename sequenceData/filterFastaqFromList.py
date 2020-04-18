# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:28:08 2019

@author: dbickhart
"""

import argparse
import sys
from fastqStats import fastq_reader_fh
from calcGCcontentFasta import read_fasta
import contextlib
import gzip

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Filter fastq/a file using a list of read names"
            )
    parser.add_argument('-f', '--fasta', 
                        help="The input sequence fasta file",
                        type=str, required=False, default = "None"
                        )
    parser.add_argument('-q', '--fastq', 
                        help="The input sequence fastq file",
                        type=str, required=False, default = "None"
                        )
    parser.add_argument('-o', '--output',
                        help="Output fastq/fasta file",
                        type=str, required=True
                        )
    parser.add_argument('-l', '--list',
                        help="tab delimited file where the first column is the read name",
                        type=str, required=True
                        )
    parser.add_argument('-v', '--reverse',
                        help="Exclude reads in this list",
                        action='store_true'
                        )
    return parser.parse_args()

def main(args):
    fastq = True if args.fastq != "None" else False
    fasta = True if args.fasta != "None" else False
    exclude = True if args.reverse else False
    
    if (fastq and fasta) or (not fastq and not fasta):
        print("Error! Must specify exactly one fastq or fasta file!")
        sys.exit()
        
    # Read in list of reads to filter
    names = set()
    with open(args.list, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            names.add(s[0])
    
    # It's ugly, but we're going to have to use different rules for the file formats
    file = args.fastq if fastq else args.fasta
    with smartFile(file) as input, open(args.output, 'w') as out:
        if fastq:
            for head, seq, qual in fastq_reader_fh(input):
                inlist = True if head[1:] in names else False
                if (inlist and not exclude) or (not inlist and exclude):
                    out.write(f'{head}\n{seq}\n+\n{qual}\n')
        else:
            for head, seq in read_fasta(input):
                inlist = True if head in names else False
                if (inlist and not exclude) or (not inlist and exclude):
                    out.write(f'>{head}\n{seq}\n') # This is without base formating! 
        
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
    args = parse_user_input()
    main(args)
