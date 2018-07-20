# -*- coding: utf-8 -*-
"""
This is a script designed to process tab delimited files and selectively grep out 
specific entries from a flatfile or direct input
Created on Fri Jul 20 12:37:06 2018

@author: dbickhart
"""

import argparse
import os.path
import sys
import contextlib

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
                        type=str
                        )
    parser.add_argument('-l', '--list',
                        help="Either a flat file or a comma delimited list of entries to screen",
                        required=True, type=str
                        )
    parser.add_argument('-s', '--short',
                        help="[Optional, flag] Return just the name of the item, not the whole line",
                        action='store_true')
    
    return parser.parse_args()

def main(args):
    # Determine if the list is a file or a comma separated list of entries
    isFile = os.path.isfile(args.list)
    
    search = generateParseDict(args.list, isFile)
    
    # Now to parse the file and grep out only the portions in our search set
    with smartFile(args.file) as fHandle, smartFile(args.output, 'w') as o:      
        try:
            for l in fHandle:
                l = l.rstrip('\n')
                if args.ignore is not None:
                    if l.startswith(args.ignore):
                        continue
                segs = l.split()
                if len(segs) - 1 < args.column:
                    raise parserException("len(segs) - 1 < args.column", "Fewer columns than expected!", "l")
                elif segs[args.column] in search:
                    if not args.short:
                        o.write(f"{l}\n")
                    else:
                        o.write(f"{segs[args.column]}\n")
        except parserException as e:
            print(e.formatMsg)
        

        
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
    
def generateParseDict(data : str, isFile : bool):
    search = set()
    if isFile:
        with open(data, 'r') as f:
            for l in f:
               l = l.rstrip('\n')
               search.add(l)
    else:
        for x in data.split(sep=','):
            search.add(x)
    return search
    
class parserException(Exception):
    """ Exception for parsing a tab delimited file
	
	Attributes:
		expression -- message from error
        line -- line that was parsed
		message -- explanation of message
        """    
    def __init__(self, expression, line, message):
        self.expression = expression
        self.line = line
        self.message = message

    def formatMsg():
        return f'Parser error! Terminated at {self.expression}\nWith the following message: {self.message}\nAt file line:\n{self.line}\n'
    

if __name__ == '__main__':
    args = parse_user_input()
    main(args)