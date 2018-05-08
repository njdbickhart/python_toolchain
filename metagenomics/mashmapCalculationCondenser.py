"""
This is a script that takes mashmap output and estimates the percentage
coverage of query sequence

@author: dbickhart
"""

import argparse
from collections import defaultdict
from functools import reduce

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program that condenses mashmap alignments to estimate coverage"
            )
    parser.add_argument('-i', '--input', 
                        help="Mashmap input file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename",
                        required=True, type=str
                        )
    
    return parser.parse_args()

def subtract_segs(qlen, qlist):
    # List of unaligned query segs
    qunalign = [] # [start, end]
    # List of aligned reference segs
    ralign = []
    
    #qunalign[0] = [1, qlen]
    qunalign.append([1, qlen])
    for a in qlist:
        ralign.append('-'.join([a.rchr, str(a.rstart), str(a.rend)]))
        if a.qstart > qunalign[-1][0]:
            # We need to resection the qualign
            prevend = qunalign[-1][1]
            qunalign[-1][1] = a.qstart
            qunalign.append([a.qend, prevend])
        
        if a.qend == qlen:
            # Remove the last section
            qunalign.pop()
        elif a.qend > qunalign[-1][0]:
            qunalign[-1][0] = a.qend
            
    qunalign = [[0,0]] if not qunalign else qunalign
    unalnlen = 0
    for row in qunalign:
        unalnlen += row[1] - row[0]

    #unalnlen = reduce(lambda x,y: x + y, map(lambda x,y: y - x, qunalign)) if len(qunalign) > 1 else map(lambda x,y: y - x, qualign)
    lastlist = []
    for start, end in qunalign:
        lastlist.append('-'.join([str(start), str(end)]))
    return processedSeg(qlist[0].qchr, qlen, lastlist, unalnlen, ralign)
        
def main(args):
    # Table of chr lengths and arrays of alignments
    queryLens = {}
    aligns = defaultdict(list)
    
    with open(args.input, "r") as f:
        for l in f:
            l.rstrip("\n")
            segs = l.split()
            queryLens[segs[0]] = int(segs[1])
            aligns[segs[0]].append(align(segs[0], int(segs[2]), int(segs[3]), segs[4],
                                segs[5], int(segs[7]), int(segs[8])))
    
    segsfile = args.output + '.segs.tab'
    with open(segsfile, "w") as o:
        for query, alist in aligns.items():
            # Process each query chromosome in turn
            alist.sort()
            processed = subtract_segs(queryLens[query], alist)
            
            o.write("{}\n".format(processed.printOut()))
            
    print("Finished with file processing")

class align:
    def __init__(self, qchr = None, qstart = -1, qend = -1, orient = "-", 
                 rchr = None, rstart = -1, rend = -1):
        self.qchr = qchr
        self.qstart = qstart
        self.qend = qend
        self.orient = orient
        self.rchr = rchr
        self.rstart = rstart
        self.rend = rend
        
    def __lt__(self, other):
        return self.qstart < other.qstart
        
class processedSeg:
    def __init__(self, qchr = None, origlen = -1, unalgSegs = [], unalgLen = -1, algnSegs = []):
        self.qchr = qchr
        self.origlen = origlen
        self.unalgSegs = unalgSegs
        self.unalgLen = unalgLen
        self.algnSegs = algnSegs
        
    def __lt__(self, other):
        return self.qchr < other.qchr
    
    def printOut(self):
        return '\t'.join([self.qchr, str(self.origlen), str(self.unalgLen), 
                         ';'.join(self.unalgSegs),
                         ';'.join(str(x) for x in self.algnSegs)])
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
