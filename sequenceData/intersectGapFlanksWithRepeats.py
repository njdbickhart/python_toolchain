#!/usr/bin/python3
# This is a script designed to intersect gap flanking sequence anlaysis with repeats
# The goal is to try to predict the formation mechanism of gaps

import re
import argparse
from collections import defaultdict
import subprocess as sp

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate associations of gap flanking regions with repetitive elements"
            )
    parser.add_argument('-g', '--gap',
                        help="The input gap association file.",
                        type=str, required=True
                        )
    parser.add_argument('-r', '--repeat',
                        help="The input repeat bed file.",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Output file base name.",
                        type=str, required=True
                        )

    return parser.parse_args()

def main(args):
    workhorse = regions(args.repeat)

    # Read file and load into containers
    with open(args.gap, 'r') as input:
        for l in input:
            l = l.rstrip()
            workhorse.loadCoords(l)

    # Process left, middle and right hand coords
    for i in [1,2,3]:
        workhorse.runComp(i)

    # print them out
    workhorse.printOut(args.output)

class regions:

    def __init__(self, repeat):
        self.repeat = repeat
        self.regions = defaultdict(list)
        # list entries: 0 -> type {close, trans}, 1-> left side, 2-> center, 3-> right side
        self.output = defaultdict(list)
        # Same order as above
        self.sre = re.compile('[:-]')
        self.clip = re.compile('Clip_')

    def loadCoords(self, line : str):
        segs = line.split()

        loc = '-'.join(segs[1:3])
        segs[0] = re.sub(self.clip, '', segs[0])
        self.regions[loc].append(segs[0])
        if segs[0] != "Unmapped":
            for i in [6,-1,7]:
                if i > 0:
                    tcord = re.split(self.sre, segs[i])
                    temp = coord(tcord[0], int(tcord[1]), int(tcord[2]))
                else:
                    if segs[0] != "Trans":
                        fcord = re.split(self.sre, segs[6])
                        rcord = re.split(self.sre, segs[7])
                        tlist = [int(x) for x in [fcord[1:2], rcord[1:2]]]
                        tlist = tlist.sort()
                        temp = coord(fcord[0], tlist[1], tlist[2])
                    else:
                        temp = None

                self.regions[loc].append(temp)

    def runComp(self, idx : int = 1):
        with open('temp.bed', 'w') as out:
            for k, v in self.regions.items():
                if v[0] != "Trans":
                    out.write(v[idx].getStr)
                    out.write(f'\t{k}\n')
                elif v[0] == "Trans" and idx != 2:
                    out.write(v[idx].getStr)
                    out.write(f'\t{k}\n')
                self.output[k] = [v[0], None, None, None]

        sp.run(f'module load bedtools; bedtools intersect -a {self.repeat} -b temp.bed -wa -wb > temp.out', shell=True)

        with open('temp.out', 'r') as input:
            for line in input:
                line = line.rstrip()
                segs = line.split()

                self.output[segs[6]][idx] = segs[3]

    def printOut(self, outfile : str):
        with open(outfile, 'w') as out:
            for k, v in self.output.items():
                out.write(k + "\t" + '\t'.join(v) + '\n')

class coord:

    def __init__(self, chr : str, start : int, end : int):
        self.chr = chr
        self.start = start
        self.end = end

    def getStr(self) -> str:
        return f'{self.chr}\t{self.start}\t{self.end}'

if __name__ == "__main__":
    args = parse_user_input()
    main(args)
