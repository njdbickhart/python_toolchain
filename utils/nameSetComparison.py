# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:09:21 2024

@author: Derek.Bickhart
"""

import argparse
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A command line script to compare the contents of two files"
            )
    parser.add_argument('-f', '--file', 
                        help="Iput Comparison file",
                        action='append', default=[]
                        )
    parser.add_argument('-o', '--output',
                        help="Output type [stdout, upset]",
                        default='stdout', type=str,
                        )
    return parser.parse_args(), parser

def main(args, parser):
    # Make sure that the input conforms to expectations
    if len(args.file < 2):
        print('Must provide more than 1 file for comparison')
        print(f'Only {len(args.file)} files were given!')
        print(parser.print_usage())
    
    # Create workhorse object
    manager = entryManager(args.file)
    
    # Update entries from list of files
    for f in args.file:
        manager.update_entries(f)
        
    # Print to Stdout for now
    manager.printStdout()
    
class entryManager():
    
    def __init__(self, filelist):
        self.fileconverter = {v : k for k, v in enumerate(filelist)}
        self.filelist = filelist
        self.entrydict = dict()
        
        # Output specific lists
        self.categories = defaultdict(int)
        
    def update_entries(self, file):
        filenum = self.fileconverter[file] + 1
        with open(file, 'r') as input:
            for l in input:
                # Assuming just the first column is used for now
                s = l.rstrip().split()
                if s[0] not in self.entrydict:
                    self.entrydict(entry(s[0],filenum))
                else:
                    self.entrydict[s[0]].add(filenum)
    
        
    def _count_entries(self):
        for k, v in self.entrydict.items():
            self.categories[v.get_ownerstring()] += 1
        
    def printStdout(self):
        # First, print out filenames and assignment
        print('Files:')
        for i, f in enumerate(self.filelist):
            print(f'File {i + 1}: {f}')
        print()
        
        # Next, count entries
        self._count_entries()
        
        # print out ordered list of entries
        cats = [k for k, v in self.categories.items()]
        cats.sort()
        for c in cats:
            print(f'{c}\t{self.categories[c]}')
        
    
class entry():
    
    def __init__(self, name, filenum):
        self.name = name
        self.owners = [filenum]
        
    def add(self, filenum):
        self.owners.append(filenum)
        
    def get_ownerstring(self):
        self.owners.sort()
        return ";".join(self.owners)

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)