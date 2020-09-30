# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:35:25 2020

@author: derek.bickhart-adm
"""

import argparse
import sys
from collections import defaultdict
import os
from os.path import basename
import subprocess as sp
import upsetplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A wrapper script for my modified version of meryl"
            )
    parser.add_argument('-m', '--meryl', 
                        help="The full path to the meryl executable.",
                        default="meryl", type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Produces a text file (.raw.tab) and a plot (.pdf)",
                        required=True, type=str,
                        )
    parser.add_argument('-d', '--db',
                        help="Meryl database folder name. At least two are required",
                        action="append", default=[],
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    if len(args.db) < 2:
        print("Must provide two or more meryl databases!")
        parser.print_usage()
        sys.exit(-1)
        
    worker = merylWrapper(args.meryl, args.db)
    
    worker.run(args.output)
    
    print(f'Finished! Output is in {args.output}')

class merylWrapper:
    
    def __init__(self, meryl, dbs):
        self.meryl = meryl 
        self.dbs = dbs 
        self.numdbs = len(dbs)
        
        # Non final attributes
        self.counter = defaultdict(int)
        
    def run(self, output):
        dcount = 0
        dbstr = " ".join(self.dbs)
        
        if os.path.exists(output + ".raw.tab"):
            print("Starting from previous task")
            with open(output + ".raw.tab", 'r') as input:
                for l in input:
                    s = l.rstrip().split()
                    self.counter[s[0]] = int(s[1])
        else:
            with sp.Popen(f'{self.meryl} print venn {dbstr}', shell=True, stdout=sp.PIPE, bufsize=1, universal_newlines=True) as sf:
                for h in sf.stdout:
                    s = h.split()
                    self.counter[s[1]] += 1
                    dcount += 1
                    if dcount % 10000000 == 0:
                        print(f'Progress: {dcount}')
                    
            # print out raw data
            with open(output + ".raw.tab", 'w') as out:
                for w in sorted(self.counter, key=self.counter.get, reverse=True):
                    out.write(f'{w}\t{self.counter[w]}\n')
                    
            print("Created raw output file")
        
        # Prepare membership df
        array = list()
        data = list()
        for k, v in self.counter.items():
            tlist = list()
            for i, e in enumerate(self.dbs):
                if (int(k) & (1 << int(i))):
                    # The bit is set, add the file name
                    tlist.append(basename(e).split('.')[0])
            array.append(tlist)
            data.append(v)
            
        # Plot things out
        dataset = upsetplot.from_memberships(array, data=data)
        print(dataset)
        
        upset = upsetplot.UpSet(dataset, sort_by='cardinality', show_percentages=True)
        upset.plot()
        
        plt.savefig(output + ".pdf")       


if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
