# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:01:20 2020

@author: derek.bickhart-adm
"""

import argparse
import gzip

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to extract subreads from PacBio CCS read datasets"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input gzipped CCS fasta file.",
                        required=True, type=str
                        )
    parser.add_argument('-q', '--fastq', 
                        help="Input gzipped subread fastq file.",
                        required=True, type=str
                        )
    parser.add_argument('-i', '--interval', 
                        help="The number of reads to sample after the first (inclusive; zerobased)",
                        default=3, type=int
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are (basename).(count).fastq",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Extract CCS read names from final CCS file
    ccsreads = dict()
    with gzip.open(args.fasta, 'rt') as ccs:
        for l in ccs:
            if l.startswith('>'):
                s = l.rstrip().split('/')
                ccsreads[s[1]] = 0
                
    # Now to read the subread fastq and process reads that match the dictionary
    outs = dict()
    for i in range(1, args.interval + 1):
        outs[i] = open(f'{args.output}.{i}.fastq', 'w')
        
    with gzip.open(args.fastq, 'rt') as fqfh:
        for name, seq, qual in fastq_reader_fh(fqfh):
            s = name.rstrip().split('/')
            if s[1] in ccsreads:
                if ccsreads[s[1]] != 0 and ccsreads[s[1]] <= args.interval:
                    outs[ccsreads[s[1]]].write(f'{name}\n{seq}\n+\n{qual}\n')
                ccsreads[s[1]] += 1
                
    for k, v in outs.items():
        v.close()
    
def fastq_reader_fh(infile):
  name = infile.readline().rstrip()
  while True:
    seq = ""
    for s in infile:
      if s[0] == '+':
        break
      else:
        seq += s.rstrip()
    qual = ""
    for q in infile:
      if len(qual) > 0 and  q[0] == '@':
        yield name, seq, qual
        name = q.rstrip()
        break
      else:
        qual += q.rstrip()
    else:
      yield name, seq, qual
      return

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
