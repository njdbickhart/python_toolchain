# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:40:09 2021

@author: derek.bickhart-adm
"""

import argparse
import os
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to parse IsoPhase magphase humanreadable text data"
            )
    parser.add_argument('-f', '--folder', 
                        help="Folder with the haplotype files",
                        required=True, type=str
                        )
    parser.add_argument('-p', '--prefix', 
                        help="File prefix for the haplotype",
                        required=True, type=str
                        )
    parser.add_argument('-d', '--dastool', 
                        help="dastool completeness estimate",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name prefix. output = {output}.long and {output}.short",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    workhorse = Strain(args.prefix)
    
    nosnps = os.path.exists(f'{args.folder}/{args.prefix}.strain.NO_SNPS_FOUND')
    
    if not nosnps:
        with open(f'{args.folder}/{args.prefix}.strain.human_readable_by_hap.txt', 'r') as input:
            input.readline()
            for l in input:
                s = l.rstrip().split()
                workhorse.add(s)
                
        


class Strain:
    
    def __init__(self, binid):
        # Keys -> contig -> hap = number of reads
        self.binid = binid
        self.contigStrains = defaultdict(lambda : defaultdict(int))
        self.maxHapcount = 0
        self.comp = 0.0
        self.cont = 0.0
        
    def add(self, segs):
        # Dodging any haplotypes that contain "?" SNPs for now
        if "?" in segs[0]:
            return
            
        self.contigStrains[segs[2]][segs[1]] = int(segs[3])
        
    def update(self, comp, cont):
        self.comp = comp
        self.cont = cont
        
    def produceLongOut(self):
        tlist = list()
        for contig, v in self.contigStrains.items():
            for hap, count in v.items():
                tlist.append(f'{binid}\t{hap}\t{len(hap)}\t{count}\t{contig}\n')
                
        return tlist
    
    def produceShortOut(self):
        for contig, v in self.contigStrains.items():
            count = len(v.keys())
            if count > self.maxHapcount:
                self.maxHapcount = count 
                
        return f'{binid}\t{self.maxHapcount}\t{self.comp}\t{self.cont}\n'

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
