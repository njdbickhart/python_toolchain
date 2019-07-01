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
    workhorse.printOut(args.output + ".comp")
    
    # Calc statistics
    workhorse.calcStats(args.output + ".stats")

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

        loc = '-'.join(segs[1:4])
        segs[0] = re.sub(self.clip, '', segs[0])
        if segs[0] != "Unmapped":
            self.regions[loc].append(segs[0])
            self.output[loc] = [segs[0], None, None, None]
            for i in [6,-1,7]:
                if i > 0:
                    tcord = re.split(self.sre, segs[i])
                    temp = coord(tcord[0], int(tcord[1]), int(tcord[2]))
                else:
                    if segs[0] != "Trans":
                        fcord = re.split(self.sre, segs[6])
                        rcord = re.split(self.sre, segs[7])
                        tlist = []
                        for i in fcord[1:] + rcord[1:]:
                            tlist.append(int(i))
                        tlist.sort()
                        temp = coord(fcord[0], tlist[1], tlist[2])
                    else:
                        temp = None

                self.regions[loc].append(temp)

    def runComp(self, idx : int = 1):
        with open('temp.bed', 'w') as out:
            for k, v in self.regions.items():
                if len(v) < idx + 1:
                    print(k)
                if v[0] != "Trans":
                    out.write(v[idx].getStr())
                    out.write(f'\t{k}\n')
                elif v[0] == "Trans" and idx != 2:
                    out.write(v[idx].getStr())
                    out.write(f'\t{k}\n')

        sp.run(f'module load bedtools; bedtools intersect -a {self.repeat} -b temp.bed -wa -wb > temp.out', shell=True)

        with open('temp.out', 'r') as input:
            for line in input:
                line = line.rstrip()
                segs = line.split()
                if len(segs) < 7:
                    continue
                self.output[segs[-1]][idx] = segs[5]

    def printOut(self, outfile : str):
        with open(outfile, 'w') as out:
            for k, v in self.output.items():
                out.write(k + "\t" + '\t'.join([str(x) for x in v]) + '\n')
                
    def calcStats(self, outfile : str):
        cats = defaultdict(dict)
        catkeys = ["RepeatConsistent", "RepeatComplex", "None"]
        for k, v in self.output.items():
            isNone = True if v[1] == v[2] == v[3] == None else False
            isSame = True if v[1] == v[2] == v[3] else False
            noMid = True if v[1] == v[3] else False
            
            # Use booleans to count categories
            if isNone:
                cats[v[0]]['None'] = cats[v[0]].get('None', 0) + 1
            elif isSame:
                cats[v[0]]['RepeatConsistent'] = cats[v[0]].get('RepeatConsistent', 0) + 1
            elif noMid and v[0] == 'Trans':
                cats[v[0]]['RepeatConsistent'] = cats[v[0]].get('RepeatConsistent', 0) + 1
            else:
                cats[v[0]]['RepeatComplex'] = cats[v[0]].get('RepeatComplex', 0) + 1
        
        # Generate a total tally and print out to file
        with open(outfile, 'w') as out:
            out.write("Class\t" + '\t'.join(catkeys) + '\n')
            totals = dict()
            for k, v in cats.items():
                out.write(k)
                for j in catkeys:
                    totals[j] = totals.get(j, 0) + v[j]
                    out.write("\t" + str(v[j]))
                out.write("\n")
            
            # Now print out the totals
            out.write("Total")
            for k in catkeys:
                out.write("\t" + str(totals[k]))
                
            out.write("\n")
            
                

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
