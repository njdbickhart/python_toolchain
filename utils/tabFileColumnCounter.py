# -*- coding: utf-8 -*-
"""
This is a port of my perl column counter script
Created on Mon Aug  6 12:30:43 2018

@author: dbickhart
"""

import argparse
import os.path
import sys
import contextlib
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A selective tab file grep analog script"
            )
    parser.add_argument('-f', '--file', 
                        help="Input tab delimited file or \"stdin\"",
                        required=True, type=str
                        )
    parser.add_argument('-c', '--column',
                        help="The tab file column to parse \[starts with zero\]",
                        required=True, type=int
                        )
    parser.add_argument('-o', '--output',
                        help="[Optional] Text output file name",
                        type=str, default="stdout"
                        )
    parser.add_argument('-i', '--ignore',
                        help="[Optional] Ignore lines that begin with this comment",
                        type=str, default = "None"
                        )
    parser.add_argument('-m', '--markdown',
                        help="Markdown flag. Formats output into table format",
                        action='store_true'
                        )
    
    return parser.parse_args()

def main(args):
    # generate main worker class
    worker = columnCounter(args.column, args.markdown, delim = args.delim, out = args.output, ignore = args.ignore)
    
    # Now, read the file
    worker.readFile(args.file)
    
    # And print out the results (either to file or stdout)
    worker.writeOut()

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'r'):
    if filename == 'stdin' or filename == 'stdout':
        if filename == 'stdin':
            fh = sys.stdin
        else:
            fh = sys.stdout
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename is not 'stdin' and filename is not 'stdout':
            fh.close()
            
class columnCounter:
    
    def __init__(self, colnum : int, mkdwn : bool = False, numeric : bool = False, delim : str = "\s{1}", out : str = "stdout", ignore : str = "None"):
        self.colnum = colnum
        self.mkdwn = mkdwn
        self.numeric = numeric
        self.delim = delim
        self.out = out
        self.ignore = ignore
        
        # Set up counters and other variables
        self.counter = defaultdict(int)
        
    def readFile(self, file : str):
        if not os.path.isfile(file):
            print(f'Error! Could not open {file} file for counting!')
            sys.exit(-1)
        with smartFile(file, 'r') as fh:
            for l in fh:
                l = l.rstrip()
                
                if self.ignore != "None":
                    if l.startswith(self.ignore):
                        next
                segs = l.split(sep=self.delim)
                
                if len(segs) < self.colnum + 1:
                    next # Ignore lines where the column length is less than one
                    
                self.count[segs[self.colnum]] += 1
                
    def writeOut(self):
        with smartFile(self.out, 'w') as o:
            if self.mkdwn:
                collen = 5
                conlen = 5
                for k, v in self.counter.items():
                    collen = len(k) if len(k) > collen else collen
                    conlen = len(v) if len(v) > conlen else conlen
                    
                collen += 1
                conlen += 1
                
                sepstr = '-' * (collen - 1)
                sepcon = '-' * (conlen - 1)
                
                o.write('|{0: <{collen}}|{1: >{conlen}}|'.format("Entry", "Value", collen=collen, conlen=conlen))
                o.write('|:{}|{}:|'.format(sepstr, sepcon))
                
                for w in sorted(self.counter, key=self.counter.get, reverse = True):
                    o.write('|{0: <{collen}}|{1: >{conlen}}|'.format(w, self.counter[w], collen=collen, conlen=conlen))
            else:
                o.write(f'Entry\tValue')
                for w in sorted(self.counter, key=self.counter.get, reverse = True):
                    o.write('{}\t{}'.format(w, self.counter[w]))

if __name__ == '__main__':
    # I hate how Python does this!
    args = parse_user_input()
    main(args)