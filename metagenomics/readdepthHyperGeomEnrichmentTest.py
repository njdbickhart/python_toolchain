# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 09:25:14 2018

@author: dbickhart
"""

import argparse
from scipy.stats import hypergeom, fisher_exact
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
    parser.add_argument('-m', '--meta',
                        help="[Optional] Metadata columns (zero-based) for the sample groups, separated by commas",
                        type=str, default="None"
                        )
    parser.add_argument('-o', '--output',
                        help="Text output file name",
                        type=str, required=True
                        )
    parser.add_argument('-g', '--group',
                        help="The grouping column (zero-based) that will be used to segregate the samples",
                        type=int, required=True
                        )
    parser.add_argument('-e', '--exact',
                        help="Use a Fisher's exact test instead of the hypergeometric test",
                        action='store_true'
                        )
    
    return parser.parse_args()

def main(args):
    # Start with separating column information
    cols = [int(x) for x in args.column.split(",")]
    hasmeta = False
    if args.meta != "None":
        meta = [int(x) for x in args.meta.split(",")]
        hasmeta = True
    
    # open files and group entries
    totcontigs = 0
    data = defaultdict(list)
    mstr = dict()
    with open(args.file, 'r') as fh:
        head = fh.readline()
        for l in fh:
            l = l.rstrip()
            segs = re.split("\t", l)
            totcontigs += 1
            data[segs[args.group]].append(container(float(segs[args.sample]), 
                 [float(segs[x]) for x in cols], segs[args.group]))
            if hasmeta:
                mstr[segs[args.group]] = [segs[x] for x in meta]

    # Normalize data by proportion of aligned reads
    valtot = 0
    othertot = [0 for x in cols]
    groupcount = 0
    for g, clist in data.items():
        groupcount += 1
        for c in clist:
            valtot += c.val
            for i in range(0, len(cols)):
                othertot[i] += c.data[i]

    # Now to recalculate rankings based on proportion of dataset aligned reads
    for g, clist in data.items():
        for c in clist:
            c.calc_prop(valtot, othertot)
    
    # calculate per group stats
    pvalues = defaultdict()
    ranks = defaultdict(dict)
    with open(args.output, 'w') as out:
        if hasmeta:
            for x in meta:
                out.write(f'metacol{x}\t')
        out.write("Group\tTotCount\tNumSampisMax\tNumSampisMin\tpMax\tpMin\tMaxSig\tMinSig\n")
        totgreat = 0
        totlower = 0
        for g, clist in data.items():
            M = len(clist)
            n_lower = len(list(filter(lambda x: x.lesser, clist)))
            n_greater = len(list(filter(lambda x: x.greater, clist)))

            p_lower = 1.0
            p_greater = 1.0
            
            pvalues[g] = [n_lower, n_greater, p_lower, p_greater]
            totgreat += n_greater
            totlower += n_lower

        for g, clist in data.items():
            N = len(clist)
            n_lower = pvalues[g][0]
            n_greater = pvalues[g][1]

            p_lower = hypergeom.sf(n_lower - 1, totcontigs, totlower, N)
            p_greater = hypergeom.sf(n_greater - 1, totcontigs, totgreat, N)
            pvalues[g][2] = p_lower
            pvalues[g][3] = p_greater
            ranks["lower"][p_lower] = 1
            ranks["upper"][p_greater] = 1

        # Calculate Benjamini-Hockberg correction
        cutoff = {}
        totentries = len(list(data.keys()))
        for i in ["lower", "upper"]:
            keylist = list(ranks[i].keys())
            keylist.sort()
            # This needs to be the value at the index of 0.05 of the total count
            index = int(totentries * 0.05)
            if index + 1 > len(keylist):
                cutoff[i] = 0.0
            else:
                cutoff[i] = keylist[int(totentries * 0.05)]
            print(f'{i} value cutoff is > {cutoff[i]}')

        # Now we can print
        for g, clist in data.items():
            M = len(clist)
            msig = pvalues[g][3] <= cutoff["upper"]
            lsig = pvalues[g][2] <= cutoff["lower"]
            
            if hasmeta:
                out.write("\t".join(mstr[g]))
                out.write("\t")
            out.write(f'{g}\t{M}\t{pvalues[g][1]}\t{pvalues[g][0]}\t{pvalues[g][3]}\t{pvalues[g][2]}\t{msig}\t{lsig}\n')
            
    exit
            

class container:
    
    def __init__(self, val : float, others : List, group : str) -> None:
        self.val = val
        self.data = others
        self.group = group
        
    def calc_prop(self, val : int, others : List) -> None:
        self.val /= val
        for i in range(0,len(others)):
            self.data[i] /= others[i]

        first, last = sorted(self.data)[0],sorted(self.data)[-1]
        self.greater = True if self.val >= last else False
        self.lesser = True if self.val <= first else False
        
        
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
