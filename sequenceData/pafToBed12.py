# -*- coding: utf-8 -*-
"""
Created on Wed May  4 13:29:19 2022

@author: derek.bickhart-adm
"""

import argparse
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to convert RNA alignments into UCSC bed file conventions"
            )
    parser.add_argument('-p', '--paf', 
                        help="Input paf file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    
    data = defaultdict(dict)
    with open(args.paf, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            name = s[0]
            chrom = s[5]
            start = int(s[7])
            end = int(s[8])
            orient = s[4]
            
            if name not in data[chrom]:
                data[chrom][name] = query(name, chrom)
            else:
                if abs(data[chrom][name].minstart - start) > 2000000:
                    continue
            data[chrom][name].addCoord(start, end, orient)
            
    # Now to go through each chromosome and sort by coordinates
    with open(args.output, 'w') as out:
        for c in data.keys():
            temp = sorted(data[c].items(), key = lambda kv: kv[1].minstart )
            for n, q in temp:
                out.write("\t".join(q.generateLine()) + "\n")
                

class query:
    
    def __init__(self, name, chromosome):
        self.name = name
        self.chromosome = chromosome
        self.coords = []
        
        self.minstart = 10000000000
        self.maxend = 0
        
    def addCoord(self, start, end, orientation):
        self.coords.append(coords(start, end, orientation))
        if start < self.minstart:
            self.minstart = start
        if end > self.maxend:
            self.maxend = end
        
    def generateLine(self):
        fixedCoords = []
        self.coords.sort()
        for c in self.coords:
            if len(fixedCoords) == 0:
                fixedCoords.append(c)
            elif c.start < fixedCoords[-1].end and c.end > fixedCoords[-1].start:
                fixedCoords[-1].start = min([c.start, fixedCoords[-1].start, c.end, fixedCoords[-1].end])
                fixedCoords[-1].end = max([c.start, fixedCoords[-1].start, c.end, fixedCoords[-1].end])
            else:
                fixedCoords.append(c)
        self.coords = fixedCoords
        for c in self.coords:
            c.makeRelative(self.minstart)

        blockCount = str(len(self.coords))
        blockSizes = [str(x.end - x.start) for x in self.coords]
        blockStarts = [str(x.start) for x in self.coords]

        orient = self.getOrient()
        return [self.chromosome, str(self.minstart), str(self.maxend), self.name, '1000', 
                          orient, str(self.minstart), str(self.maxend), '0,0,255', blockCount,
                          ",".join(blockSizes), ",".join(blockStarts)]
    
    def getOrient(self):
        oKeys = defaultdict(int)
        for o in self.coords:
            oKeys[o.orientation] += 1
        if len(oKeys.keys()) == 1:
            # If there is a consensus, return that
            return list(oKeys.keys())[0]
        else:
            # otherwise, return the most frequent key
            tlist = [k for k, v in sorted(oKeys.items(), key=lambda item: item[1], reverse = True)]
            return tlist[0]
        
    def __lt__(self, other):
        return self.minstart < other.minstart

class coords:
    
    def __init__(self, start, end, orientation):
        self.start = start
        self.end = end
        self.orientation = orientation

    def makeRelative(self, rStart):
        self.start = self.start - rStart
        self.end = self.end - rStart
        
    def __lt__(self, other):
        return self.start < other.start

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
