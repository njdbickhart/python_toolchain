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
                        required=False, type=str, default = "NO"
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
        
        if args.dastool != "NO":
            with open(args.dastool, 'r') as das:
                das.readline()
                for l in das:
                    s = l.rstrip().split()
                    if s[0] == args.prefix:
                        workhorse.update(s[-2], s[-1])
                        break
                
        with open(args.output + '.long', 'w') as long, open(args.output + '.short', 'w') as short:
            tlist = workhorse.produceLongOut()
            for t in tlist:
                long.write(t)
            short.write(workhorse.produceShortOut())
    else:
        with open(args.output + '.long', 'w') as long, open(args.output + '.short', 'w') as short:
            long.write(f'{args.prefix}\t0\t0\t0\t0\t0\t0\n')
            short.write(f'{args.prefix}\t0\t0\t0\n')
        


class Strain:
    
    def __init__(self, binid):        
        self.binid = binid
        # Keys -> contig -> hap = number of reads
        self.contigStrains = defaultdict(lambda : defaultdict(int))
        # Keys -> contig -> hap -> list of positions
        self.contigPositions = defaultdict(lambda : defaultdict(list))
        self.maxHapcount = 0
        self.comp = 0.0
        self.cont = 0.0
        
    def add(self, segs):
        # Dodging any haplotypes that contain "?" SNPs for now
        if "?" in segs[0]:
            return
            
        self.contigStrains[segs[2]][segs[0]] = int(segs[3])
        
    def position(self, segs):
        if "?" in segs[0]:
            return
        
        self.contigPositions[segs[2]][segs[0]].append(segs[3])
        
        
    def update(self, comp, cont):
        self.comp = float(comp)
        self.cont = float(cont)
        
    def produceLongOut(self):
        tlist = list()
        for contig, v in self.contigStrains.items():
            for hap, count in v.items():
                pos = [-1, -1]
                if hap in self.contigPositions[contig]:
                    tlist = self.contigPositions[contig][hap]
                    minp = min(tlist)
                    maxp = max(tlist)
                    pos = [minp, maxp]

                tlist.append(f'{self.binid}\t{hap}\t{len(hap)}\t{count}\t{contig}\t{pos[0]}\t{pos[1]}\n')
                
        return tlist
    
    def produceShortOut(self):
        for contig, v in self.contigStrains.items():
            count = len(v.keys())
            if count > self.maxHapcount:
                self.maxHapcount = count 
                
        return f'{self.binid}\t{self.maxHapcount}\t{self.comp}\t{self.cont}\n'

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
