# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:30:18 2019

@author: dbickhart
"""

import argparse
import defaultdict
import re

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Identify regions of potential inverted contigs in assemblies"
            )
    parser.add_argument('-f', '--file', 
                        help="The input chimera sam file",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output', 
                        help="The output file basename",
                        type=str, required=True
                        )
    return parser.parse_args()

def main(args):
    workhorse = coordIntersection()
    
    with open(args.file, 'r') as input:
        for l in input:
            l = l.rstrip()
            s = l.split()
            workhorse.loadCoordinate(s)
            
    workhorse.condenseCoords()
    
    workhorse.printOutList(args.output + ".allsites.bed")
    
    workhorse.slidingWindow(args.output + ".windowsites.bed")
   

class coordIntersection:

    def __init__(self):
        self.cigarRE = re.compile(r'\d+[MIDNSHP=X]')
        self.coords = defaultdict(list)
        
    def loadCoordinate(self, samStrs : list):
        # Only append coordinates where the two mates are on the same chromosome
        # Also, only if the samflags meet the "same" orientation requirement
        f = int(samStrs[1])
        if samStrs[2] == samStrs[6] and ((f & 16 == 16 and f & 32 == 32) or (f & 16 != 16 and f & 32 != 32)):
            self.coords[samStrs[2]].append(samCoord(samStrs[2], int(samStrs[3]), self.cigarRE, samStrs[5]))
   
    def condenseCoords(self) -> None:
        # do an overlap and create new coord class from the data
        for c, g in self.coords.items():
            g.sort()
            
            # Load the first element into our temporary list
            temp = [samCoord(c, g[0].start, self.cigarRE, "", g[0].end)]
            # normal consensus logic -- if the end doesn't overlap discard it
            for i in range(1, len(g)):
                cur = temp[-1]
                if cur.end > g[i].start and cur.start < g[i].end:
                    #overlap
                    cur.end = g[i].end
                    cur.counter += 1
                else:
                    # make new entry
                    temp.append(samCoord(c, g[i].start, self.cigarRE, "", g[i].end))
            
            self.coords[c] = temp # Load the condensed list back in the attribute
                    
    def printOutList(self, filename : str):
        with open(filename, 'w') as out:
            for c in sorted(self.coords.keys()):
                for g in sorted(self.coords[c]):
                    t = self.coords[c][g]
                    out.write(f'{c}\t{t.start}\t{t.end}\t{t.counter}\n')
                    
    def slidingWindow(self, filename : str, win = 4, cisd = 1000, transd = 5000):
        with open(filename, 'w') as out:
            # work with data in set windows
            # Assumption: there should be a small distance between read pileups on one side (cisd)
            # and there should be a large distance between two sets of closely spaced read pileups (transd)
            for c in sorted(self.coords.keys()):
                counts = {'Conforms' : 0, 'Close' : 0}
                for i in range(len(self.coords[c]) - win):
                    temp = []
                    for j in range(win):
                        temp.append(self.coords[c][i + j])
                    
                    front = True if temp[1].start - temp[0].end <= cisd else False
                    back = True if temp[-2].end - temp[-1].start <= cisd else False
                    middle = True if temp[-1].start - temp[0].end <= transd else False
                    
                    if front and back and middle:
                        out.write(f'{c}\t{temp[0].end}\t{temp[-1].start}\tConforms\n')
                        counts['Conforms'] += 1
                    elif sum([front, back, middle]) > 1:
                        out.write(f'{c}\t{temp[0].end}\t{temp[-1].start}\tClose\t{front}\t{middle}\t{back}\n')
                        counts['Close'] += 1
                
                print(f'Stats for {c}: Conforms: {counts["Conforms"]} Close: {counts["Close"]}')
            
    
class samCoord:
    
    def __init__(self, chrom : str, pos : int, cigarRe, cigar : str, end = 0):
        self.chrom = chrom
        self.start = pos
        self.counter = 1
        self.cigarRE = cigarRe
        if end == 0:
            self.end = self._getEndCoord(pos, cigar)
        else:
            self.end = end
        
    def _getEndCoord(self, start : int, cigar : str) -> int:
        end = start
        for cstr in self.cigarRE.findall(cigar):
            val = int(cstr[:-2])
            # We add the value of the cigar if it is part of MIS=X
            end += self._cigStrToVal(cstr) * val
        return end
            
    def _cigStrToVal(self, cstr : str) -> int:
        switcher = {
            'M' : 1,
            'I' : 1,
            'S' : 1,
            '=' : 1,
            'X' : 1,
        }
        return switcher.get(cstr, 0)
    
    def __lt__(self, other):
        # special method for sorting as a list
        return self.start < other.start
                 
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
