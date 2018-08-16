# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:20:31 2018

@author: dbickhart
"""
import argparse
import sys
from typing import Dict, List
from collections import defaultdict


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Join tab files with similar header rows together"
            )
    parser.add_argument('-f', '--file', 
                        help="An append style input of files to join together (in order). Must be specified at least twice",
                        action='append', required=True
                        )
    parser.add_argument('-c', '--column',
                        help="The columns (zero-based) by which the data should be merged (assumed to be similar among all files). Can be specified more than once",
                        action='append', required=True
                        )
    parser.add_argument('-m', '--merge',
                        help="The columns (zero-based) that should be kept from each file for the merger. Can be specified more than once. No specification means all columns except the -c columns",
                        action='append', default=[]
                        )
    parser.add_argument('-o', '--output',
                        help="Text output file name",
                        type=str, required=True
                        )
    parser.add_argument('-t', '--header',
                        help="[Flag] Is there a header string?",
                        type=bool, default=False
                        )
    
    return parser.parse_args()

def main(args):
    if len(args.file) < 2:
        print(f'Error! This program requires at least two input files! Please see the help message for more details')
        sys.exit()
    
    # Load files into classes
    workers = []
    for f in args.file:
        workers.append(colFile(f, args.head, args.cols, args.merge))
        
    for w in workers:
        w.loadFile()
        
    print(f'Loaded files. Now joining...')
    
    with open(args.output, 'w') as out:
        # Take care of any headers
        if args.header:
            select = [workers[0].header[i] for i in list(map(int, args.cols))]
            remain = [i for i in range(len(workers[0].header) -1) if i not in set(args.cols)]
            for w in workers:
                select.extend([str(w.file + w.header[i]) for i in remain])
            out.write('\t'.select + "\n")
            
        # Get sorted list of keys
        keysort = sorted(workers[0].data.keys())
        for k in keysort:
            out.write(k)
            for w in workers:
                out.write("\t" + '\t'.join(w.getValue(k)))
            out.write('\n')
            
    print(f'Finished!')
    
class colFile:
    
    def __init__(self, file : str, head : bool, cols : List, merge : List):
        self.file = file
        self.head = head
        self.cols = list(map(int, cols))
        self.merge = list(map(int, merge)) if len(merge) > 0 else []
        self.data = defaultdict(list)
        self.header = []
        
        self.max_v = max(self.cols + self.merge)
        
    def loadFile(self)->None:
        warning = False
        with open(self.file, 'r') as fh:
            for l in fh:
                l = l.rstrip()
                segs = l.split()
                
                if len(self.merge) == 0:
                    self.merge = [i for i in range(len(segs) -1 ) if i not in set(self.cols)]
                
                if len(self.header) == 0:
                    if self.head:
                        self.header = [segs[i] for i in self.cols + self.merge]
                        continue
              
                
                if len(segs) < self.max_v:
                    warning = True
                    continue
                
                key = '\t'.join([segs[x] for x in self.cols])
                self.data[key] = [segs[i] for i in self.merge]
                
        if warning:
            print(f'Skipped lines in {self.file} that did not have consistent columns!')
        
    def getValue(self, key : str) -> List:
        if key in self.data:
            return self.data[key]
        else:
            return ['-' for i in self.merge]

if __name__ == "__main__":
    args = parse_user_input()
    main(args)