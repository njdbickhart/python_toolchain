# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 12:21:02 2019

@author: dbickhart
"""

import argparse;

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Select or remove entries from a format38 file"
            )
    parser.add_argument('-f', '--file', 
                        help="The input format38 file",
                        type=str, required=True
                        )
    parser.add_argument('-l', '--list', 
                        help="A newline delimited list of entries to search for",
                        type=str, required=True
                        )
    parser.add_argument('-v', '--reverse',
                        help="[Optional, flag] Remove entries from the list, do not display them",
                        action='store_true')
    parser.add_argument('-o', '--output',
                        help="Output file",
                        type=str, required=True
                        )
    return parser.parse_args()

def main(args):
    # Create the set for lookup
    entries = set()
    with open(args.list, 'r') as input:
        for l in input:
            seventeen = l[3:20]
            entries.add(seventeen)
                    
    with open(args.file, 'r') as input, open(args.output, 'w') as out:
        for l in input:
            seventeen = l[3:20]
            contained = seventeen in entries
            if (contained and not args.reverse) or (not contained and args.reverse):
                out.write(l)
                
    print("Finished filtering. Output is in: " + args.output)
    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
