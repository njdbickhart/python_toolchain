# -*- coding: utf-8 -*-
"""
This is a script that converts the chromosome.data coordinates of Paul's 
marker sim program into different chromosome coordinates

@author: dbickhart
"""

import sys
import argparse

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program to simulate genotype markers from pedigrees"
            )
    parser.add_argument('-c', '--ctable', 
                        help="Chromosome.data file",
                        required=True, type=str
                        )
    parser.add_argument('-f', '--faidx',
                        help="Fasta index entry",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename",
                        required=True, type=str
                        )
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_user_input()
    
    # Read in fasta index and set chromosome lengths
    chrlens = {}
    with open(args.faidx, "r") as f:
        for c, s, b, l, n in f.readlines():
            chrlens[c] = s
            
    # Process chromosome.data file
    defChrLen = 100000000
    with open(args.ctable, "r") as c, open(args.output, "w") as o:
        header = c.readline()
        o.write(header)
        for line in c:
            line.rstrip('\n')
            segs = line.split("\s+")
            
            if segs[1] in chrlens:
                newpos = int(chrlens[segs[1]] * (int(segs[4]) / defChrLen))
                o.write(segs[0])
                for i, val in zip([5,9,9,12,5,9,9,9], \
                                  [segs[1], segs[2], segs[3], segs[4], newpos, \
                                   segs[5], segs[6], segs[7], segs[8]]):
                    o.write('{:>base}'.format(val, base=i))
                o.write("\n")