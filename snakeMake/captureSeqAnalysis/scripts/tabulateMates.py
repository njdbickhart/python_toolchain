# This script is designed to count and consolidate mate maps from the vector sequence
# This is a one-shot script designed to work in a snakemake workflow
from collections import defaultdict
import re
import os
import sys

usage = f'python {sys.argv[0]} <vec genome (indexed)> <vector reads (sam)> <reference reads (bed)> <outputname>'
# 1 = vector sequence file (indexed)
# 2 = Vector links
# 3 = reference links
if len(sys.argv) != 5:
    print(usage)
    sys.exit(-1)

def main(args):
    # get vector length
    vecsize = 0
    with open(args[1] + '.fai', 'r') as input:
        for l in input:
            s = l.rstrip().split()
            vecsize = int(s[1])

    # Now run through worker class
    worker = evidence(vecsize)

    # sam file first
    with open(args[2], 'r') as input:
        for l in input:
            s = l.rstrip().split()
            worker.addEvidence(s, True)

    # Now the bed file
    with open(args[3], 'r') as input:
        for l in input:
            s = l.rstrip().split()
            worker.addEvidence(s, False)

    # Cluster and organize the reads
    unpaired = worker.clusterReads()
    print(f'Identified {unpaired} unpaired reads.')

    # Print out the table
    worker.dumpAllObs(args[4])

# Class to store node data and edge weights
class node_edge:
    def __init__(self, chr, pos, side, orient):
        self.chr = chr
        self.pos = pos
        self.side = side
        self.orient = orient
        self.weight = 0

    def comp_and_increment(self, side, orient):
        data = (side, orient)
        if data == (self.side, self.orient):
            self.weight += 1
            return True
        else:
            return False
        
    def getOutStr(self):
        return "\t".join([self.chr, str(self.pos), self.side, self.orient, str(self.weight)])

class evidence:
    def __init__(self, veclength):
        self.veclength = veclength
        self.pairs = dict() # readname -> read_pair class
        self.unpaired = 0
        self.clusters = defaultdict(list) # chromosome -> (start, end)
        self.chrpos = defaultdict(lambda : defaultdict(list)) # chromosome -> position -> node_edge class

    def addEvidence(self, segs, sam):
        r = read()
        if sam:
            r.loadSam(segs)
        else:
            r.loadBed(segs)
        
        if r.name not in self.pairs:
            self.pairs[r.name] = read_pair(r.name, self.veclength)
        self.pairs[r.name].loadread(r)

    def clusterReads(self):
        # First we need to generate position coordinates for merger because I did not plan the data structures
        # very well
        for k, v in self.pairs.items():
            v.determine_orientation()
            if not v.paired:
                self.unpaired += 1
                continue
            (chr, pos) = v.getRefPos()

            # update existing clusters if possible
            if chr in self.clusters:
                found = False
                for i, d in enumerate(self.clusters[chr]):
                    if pos < d[0] and pos >= d[0] - 300:
                        self.clusters[chr][i][0] = pos
                        found = True
                    elif pos > d[1] and pos <= d[1] + 300:
                        self.clusters[chr][i][1] = pos
                        found = True
                if not found:
                    # Not found, so we create a new cluster
                    self.clusters[chr].append((pos, pos))
            else:
                self.clusters[chr].append((pos, pos))

        # Creating node_edge lists from clusters
        for chr, v in self.clusters.items():
            for i, d in enumerate(self.clusters[chr]):
                self.chrpos[chr][i] = []

        # Next we loop through the clusters and generate the nodes and edge weights
        for k, v in self.pairs.items():
            if not v.paired:
                continue
            (chr, pos) = v.getRefPos()

            for i, d in enumerate(self.clusters[chr]):
                (start, end) = d
                nedgepos = end if v.side == "left" else start
                if len(self.chrpos[chr][i]) == 0:
                    # Create new node_edge
                    self.chrpos[chr][i].append(node_edge(chr + '', nedgepos + 0, v.side + '', v.orient + ''))
                else:
                    for ne in self.chrpos[chr][i]:
                        added = ne.comp_and_increment(v.side, v.orient)
                        if not added:
                            # Create new node_edge
                            self.chrpos[chr][i].append(node_edge(chr + '', nedgepos + 0, v.side + '', v.orient + ''))
        return self.unpaired


    def dumpAllObs(self, out):
        with open(out, 'w') as output:
            output.write("chr\tpos\trefside\tvec_orient\tweight\n")
            for k, v in sorted(self.chrpos.items()):
                for i in v:
                    output.write(i.getOutStr() + "\n")
            



class read_pair:
    def __init__(self, name, veclength):
        self.pairs = dict() # dict with two entries loaded by read number
        self.name = name
        self.veclength = veclength

        # generated statistics
        self.start = False # orientation on the construct. Is it at the beginning?
        self.end = False # orientation on the construct. Is it at the end? If NEITHER, it is in the middle, so partial insertion
        self.orient = "+"  # orientation on the chromosome. Is the construct in the for or rev orientation?
        self.paired = False # If this is a proper read pair or an anchor
        self.side = "None" # If this read pair is on the right or left side of the reference sequence

    def loadread(self, read):
        self.pairs[read.num] = read

    def getRefPos(self):
        for k, v in self.pairs.items():
            if v.issam:
                continue
            else:
                return (self.pairs[k].chr, self.pairs[k].pos)

    def determine_orientation(self):
        if len(self.pairs.keys()) < 2:
            return
        else:
            self.paired = True 

        samside = None
        bedside = None
        for k, v in self.pairs.items():
            if v.issam:
                samside = self.pairs[k]
            else:
                bedside = self.pairs[k]

        # If the read alignment is on the forward side of the vector (within the first 300 bp)
        if samside.pos <= 300:
            self.start = True
        # Also check if the read alignment is on the end of the construct
        elif samside.pos >= self.veclength - 300:
            self.end = True

        # Determining insertion orientation
        orientations = (bedside.orient, samside.orient)
        if orientations == ("-", "+") and self.start:
            # disconnected join, should not be possible
            self.orient = "?"
        elif orientations == ("+", "-") and self.end:
            # disconnected join, should not be possible
            self.orient = "?"
        elif orientations == ("+", "+") and self.end:
            # ref seq <-> reversed construct
            self.orient = "-"
            self.side = "left"
        elif orientations == ("-", "-") and self.start:
            # reversed construct <-> ref seq
            self.orient = "-"
            self.side = "right"
        elif orientations == ("+", "-") and self.start:
            # ref seq <-> construct
            self.orient = "+"
            self.side = "left"
        elif orientations == ("-", "+") and self.end:
            # construct <-> ref seq
            self.orient = "+"
            self.side = "right"
        else:
            self.orient = "?"

    def printPairStats(self):
        # Name\t(read 1 str)\t(read 2 str)\tstart?\tend?\torient\tpaired\tside
        data = list()
        data.append(self.name)
        for i in range(2):
            data.append(self.pairs[i + 1].getStr())
        data.extend([self.start, self.end, self.orient, self.side])



class read:
    def __init__(self):
        self.name = "None"
        self.num = 0
        self.chr = "None"
        self.pos = 0
        self.orient = "None"
        self.issam = False # indicates if this is from the alignment to the construct

    def loadSam(self, segs):
        self.num = 1 if int(segs[1]) & 0x40 else 2
        self.name = segs[0]
        self.chr = segs[2]
        self.pos = int(segs[3])
        self.orient = "-" if int(segs[1]) & 0x10 else "+"
        self.issam = True 

    def loadBed(self, segs):
        self.num = 1 if segs[3].endswith('/1') else 2
        self.name = segs[3].replace(r'\/\d{1}', '')
        self.chr = segs[0]
        self.pos = int(segs[1])
        self.orient = segs[5]

    def getStr(self):
        return "\t".join([self.num, self.chr, self.pos, self.orient, self.issam])

if __name__ == "__main__":
    main(sys.argv)