# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:55:43 2019

@author: dbickhart
"""

import argparse;
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate an upset plot using coordinate positions from alignment files"
            )
    parser.add_argument('-f', '--fai', 
                        help="The input reference fasta index file (Required; must have same names as in alignments)",
                        type=str, required=True
                        )
    parser.add_argument('-p', '--paf', 
                        help="An input paf alignment file (if used, must be specified more than once)",
                        action="append", default=[]
                        )
    parser.add_argument('-o', '--output', 
                        help="The output file name for the insersection list (must be same length as alignment file input!)",
                        action="append", default=[]
                        )
    parser.add_argument('-w', '--window', 
                        help="The window size (in bases) for the segments [default = 1000]",
                        type=int, default=1000
                        )
    return parser.parse_args()

def main(args):
    pafs = args.paf
    outs = args.output
    
    workHorse = flatIntersection(args.window)
    workHorse.generateWindows(args.fai)
   
    for p in pafs:
        workHorse.addFile(p, 'paf')
 
    workHorse.getOutput(outs)
 
class flatIntersection:
    
    def __init__(self, window : int):
        self.data = defaultdict(dict)
        self.fileIdx = dict()
        self.catIdx = dict()
        self.fileItr = 1
        self.chrOrder = []
        self.window = window
        
    def generateWindows(self, fai : str) -> None:
        chrlens = dict()
        chrorder = list()
        with open(fai, 'r') as input:
            for l in input:
                l = l.rstrip()
                s = l.split()
                chrlens[s[0]] = int(s[1])
                chrorder.append(s[0])
        
        print("Generating windows of size: " + str(self.window))
        self.chrOrder = chrorder
        
        for idx in range(len(chrorder)):
            chrom = chrorder[idx]
            chrlen = chrlens[chrom]
            
            for i in range(int(chrlen / self.window) + 1):
                self.data[chrom][i] = "1"
    
    def getOutput(self, outputs : list):
        outhandles = {x : open(outputs[x], 'w') for x in range(len(outputs))}
        
        for chrom in self.chrOrder:
            for i in sorted(self.data[chrom].keys()):
                matches = self.data[chrom][i].split(';')
                # Convert to unique list
                unique = set(matches)
                
                for f in unique:
                    outhandles[int(f) - 1].write(f'{chrom}:{i * self.window}\n')
        
        for handle in outhandles.values():
            handle.close()
            
    def addFile(self, file : str, ftype : str) -> None:
        self.fileItr += 1
        self.fileIdx[file] = self.fileItr
        self.catIdx[self.fileItr] = file
        
        if ftype == "paf":
            print(f'Adding paf file {file}')
            self._processPaf(file)
        else:
            print("Shouldn't be here!")
        
    def _processPaf(self, file : str) -> None:
        with open(file, 'r') as input:
            for l in input:
                l = l.rstrip()
                s = l.split()
                
                start = self._convertCoord(int(s[7]))
                end = self._convertCoord(int(s[8]))
                
                for j in range(start, end + 1):
                    self.data[s[5]][j] += f';{self.fileItr}'
        
    def _convertCoord(self, pos : int) -> int:
        return int(pos / self.window)
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
