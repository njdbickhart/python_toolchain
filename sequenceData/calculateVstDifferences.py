# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:12:20 2019

@author: dbickhart
"""

import argparse
from collections import defaultdict
import numpy as np

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate restriction enzyme maps for an assembly fasta and select fragments of a specific size range"
            )
    parser.add_argument('-c', '--gene', 
                        help="Input gene copy number list from annotation program",
                        type=str, required=True
                        )
    parser.add_argument('-p', '--pop', 
                        help="Population structure (sample(tab)popnumber) file",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output', 
                        help="Output file name",
                        type=str, required=True
                        )
    parser.add_argument('-p', '--skip', 
                        help="Skip this sample (can be specified more than once)",
                        action="append", default=[]
                        )
    return parser.parse_args()

def main(args):
    # Read in population structure and skip list
    skip = set()
    if len(args.skip) > 0:
        skip = {x for x in args.skip}
    
    popstruct = defaultdict(list)
    poplookup = {}
    with open(args.pop, 'r') as fh:
        for l in fh:
            l = l.rstrip()
            s = l.split()
            if s[0] in skip: continue
            poplookup[s[0]] = s[1]
            popstruct[s[1]].append(s[0])
            
    # Now load up data structures
    

class totalHolder:
    
    def __init__(self, markername: str, values: list, gene: str):
        self.markername = markername
        self.values = []
        self.values.append(values)
        self.genes = gene
    
    # Routines for overlapping gene regions
    def addValues(self, values: list):
        self.values.append(values)
        
    def addGene(self, gene: str):
        self.genes += ';' + gene
   
    def Vst(self, genoLookup, markers, skip):
        # Total variance calculation
        total = []
        numerator = 0
        denominator = 0
        # Data structures
        popnums = {}
        popvar = defaultdict(list)
        for i in self.values:
            for j in range(0, len(self.values[i])):                
                sample, pop = genoLookup(j)
                if sample == None:
                    continue
                
                total.append(self.values[i][j])
                popnums{pop} += 1
                
        
        Vt = np.var(total)
        
        
        
class genoLookup:
    
    def __init__(self, poplookup, header, skip):
        self.pops = {}
        self.samples = {}
        
        for i in range(6, len(header)):
            if header[i] in skip:
                self.samples[i - 6] = None
                self.pops[i - 6] = None
            else: 
                self.samples[i - 6] = header[i]
                self.pops[i - 6] = poplookup[header[i]]
            
    def getPopSamp(self, idx):
        return self.samples[i], self.pops[i]

    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
