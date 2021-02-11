# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:58:31 2021

@author: derek.bickhart-adm
"""

import argparse
import os
import glob
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Generate consensus bin taxonomy using blobtools taxify results"
            )
    parser.add_argument('-c', '--conversion', 
                        help="Input taxid conversion table for NCBI taxIDs",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output table file",
                        required=True, type=str,
                        )
    parser.add_argument('-t', '--table', 
                        help="Input table of taxified contigs",
                        required=True, type=str
                        )
    parser.add_argument('-d', '--directory', 
                        help="Input directory of bin fasta files",
                        required=True, type=str
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Create conversion table for Tax IDs
    taxIDConv = dict()
    with open(args.conversion, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            taxIDConv[s[0]] = s[2]
            
    # Create contig conversion table for Tax IDs
    NAS = 0
    ctgConv = dict()
    with open(args.table, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            ctgConv[s[0]] = taxIDConv.get(s[3], 'NA')
            if ctgConv[s[0]] == 'NA':
                NAS += 1
                
    if NAS > 0:
        print(f'Found {NAS} NA values')
            
    # Grab bin fasta files, prepare list for output
    bins = list()
    for file in glob.glob(f'{args.directory}/*.fa'):
        tbin = binFile(file)
        print(f'Working on {file}')
        with open(file, 'r') as input:
            for l in input:
                l = l.rstrip()
                if l.startswith('>'):
                    tid = ctgConv.get(l[1:], 'NA')
                    tbin.update(tid)
        bins.append(tbin)
        
    # Now print it all out
    with open(args.output, 'w') as out:
        for b in bins:
            out.write(b.createOutLine())
    print("Done!")

class binFile:

    def __init__(self, filename):
        self.filename = os.path.basename(filename)
        self.totctgs = 0
        self.taxid = defaultdict(int)
        
    def update(self, taxid):
        self.taxid[taxid] += 1
        self.totctgs += 1
        
    def createOutLine(self):
        sortVs = [f'{k}:{v}' for k, v in sorted(self.taxid.items(), key=lambda item: item[1], reverse=True)]
        return "{}\t{}\t{}\t{}\n".format(self.filename, self.totctgs, len(sortVs), "\t".join(sortVs))
        

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
