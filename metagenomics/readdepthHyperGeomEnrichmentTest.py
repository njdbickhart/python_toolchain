# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 09:25:14 2018

@author: dbickhart
"""

import argparse
from scipy.stats import hypergeom
import re
from collections import defaultdict
from typing import List


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Calculate Hypergeometric two-tailed probability for read coverage of a single sample against all others using groups of specific categories"
            )
    parser.add_argument('-f', '--file', 
                        help="The input tab delimited file for testing",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--column',
                        help="The columns (zero-based) containing the additional samples, separated by commas",
                        type=str, required=True
                        )
    parser.add_argument('-s', '--sample',
                        help="The test sample column (zero-based) that will compared against the 'columns' samples",
                        type=int, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Text output file name",
                        type=str, required=True
                        )
    parser.add_argument('-g', '--group',
                        help="The grouping column (zero-based) that will be used to segregate the samples",
                        type=int, required=True
                        )
    
    return parser.parse_args()

def main(args):
    # Start with separating column information
    cols = [int(x) for x in args.column.split(",")]
    
    # open files and group entries
    data = defaultdict(list)
    with open(args.file, 'r') as fh:
        head = fh.readline()
        for l in fh:
            l = l.rstrip()
            segs = re.split("\t", l)
            data[segs[args.group]].append(container(float(segs[args.sample]), 
                 [float(segs[x]) for x in cols], segs[args.group]))
    
    # calculate per group stats
    with open(args.output, 'w') as out:
        out.write("Group\tTotCount\tNumSampisMax\tNumSampisMin\tpMax\tpMin\n")
        for g, clist in data.items():
            M = len(clist)
            n_lower = len(filter(lambda x: x.lesser, clist))
            n_greater = len(filter(lambda x: x.greater, clist))
            
            p_lower = hypergeom.sf(n_lower - 1, M, n_lower, M)
            p_greater = hypergeom.sf(n_greater - 1, M, n_greater, M)
            
            out.write(f'{g}\t{M}\t{n_greater}\t{n_lower}\t{p_greater}\t{p_lower}\n')
            
    exit
            

class container:
    
    def __init__(self, val : float, others : List, group : str) -> None:
        self.val = val
        self.data = others
        self.group = group
        
        first, last = sorted(others)[0],sorted(others)[-1]
        self.greater = True if val >= last else False
        self.lesser = True if val <= first else False
        
        
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)