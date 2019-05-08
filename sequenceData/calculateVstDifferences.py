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
    parser.add_argument('-s', '--skip', 
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
            
    # Now to process the data on a line-by-line basis
    numericErrs = 0
    with open(args.gene, 'r') as fh, open(args.output, 'w') as out:
        # Process header to get the column indicies for the samples
        head = fh.readline()
        geno = genoLookup(poplookup, head, skip)
        for l in fh:
            l = l.rstrip()
            s = l.split()
            if not isNumber(s[6]): 
                # Some of the entries are "null" because of the difference in chromosome names
                numericErrs += 1
                continue
            marker = s[1] + ':' + s[2] + '-' + s[3]
            worker = totalHolder(marker, s[6:], s[0])
            
            Vst = worker.Vst(geno, marker, skip)
            out.write(f'{s[1]}\t{s[2]}\t{s[3]}\t{Vst}\t{s[0]}\n')
            
    print(f'Dealt with {numericErrs} null fields')
    
def isNumber(value : str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False

class totalHolder:
    
    def __init__(self, markername: str, values: list, gene: str):
        self.markername = markername
        self.values = []
        self.values.append([float(x) for x in values])
        self.genes = gene
    
    # Routines for overlapping gene regions
    def addValues(self, values: list):
        self.values.append(values)
        
    def addGene(self, gene: str):
        self.genes += ';' + gene
   
    def Vst(self, genoLookup, markers, skip) -> float:
        # Total variance calculation
        total = []
        numerator = 0
        denominator = 0
        # Data structures
        popnums = {}
        popvar = defaultdict(list)
        for i in range(0, len(self.values)):
            for j in range(0, len(self.values[i])):                
                sample, pop = genoLookup.getPopSamp(j)
                if sample == None:
                    continue
                
                total.append(self.values[i][j])
                if not pop in popnums:
                    popnums[pop] = 0
                popnums[pop] += 1
                popvar[pop].append(self.values[i][j])
                
        for pop, num in popnums.items():
            pvar = np.var(popvar[pop])
            numerator += pvar * num
            denominator += num
            
        Vt = np.var(total)
        Vs = numerator / denominator
        if Vt == 0: return 0.0
        
        return ((Vt - Vs) / Vt)      
        
        
class genoLookup:
    
    def __init__(self, poplookup, head, skip):
        self.pops = {}
        self.samples = {}
        head = head.rstrip()
        header = head.split()
        
        for i in range(6, len(header)):
            if header[i] in skip:
                self.samples[i - 6] = None
                self.pops[i - 6] = None
            else: 
                self.samples[i - 6] = header[i]
                self.pops[i - 6] = poplookup[header[i]]
            
    def getPopSamp(self, idx: int):
        return self.samples[idx], self.pops[idx]

    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
