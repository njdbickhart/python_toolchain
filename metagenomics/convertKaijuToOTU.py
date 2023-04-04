# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 10:41:49 2023

@author: Derek.Bickhart
"""

import argparse
import sys
from collections import defaultdict 

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script designed to consolidate Kaiju read count data into an OTU table"
            )
    parser.add_argument('-f', '--file', 
                        help="Input Kaiju count file. Can be specified more than once. Must be accompanied by an '-s' tag",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--sample', 
                        help="Kaiju count file sample name. Can be specified more than once. Must be accompanied by an '-f' tag",
                        action="append", default=[]
                        )
    parser.add_argument('-o', '--output',
                        help="Output OTU file name",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser

def main(args, parser):
    print("hey")
    # Check the argument requirements
    proceed = False
    if len(args.file) == len(args.sample) and len(args.file) > 0:
        proceed = True
    
    if not proceed:
        print("Error! Please specify equal numbers of -f and -s in your input files!")
        print(parser.print_usage())
        sys.exit(-1)
        
    container = defaultdict(dict)
    taxset = set()
    for f, s in zip(args.file, args.sample):
        container[s] = create_tax_dict(f)
        taxset.update(list(container[s].keys()))
        
    with open(args.output, 'w') as output:
        output.write("Tax\t" + "\t".join(args.sample) + "\n")
        for t in taxset:
            output.write(f'"{t}"')
            for s in args.sample:
                count = container[s].get(t, 0)
                output.write(f'\t{count}')
            output.write('\n')

def create_tax_dict(file):
    temp = dict()
    with open(file, 'r') as input:
        h = input.readline()
        for l in input:
            s = l.rstrip().split("\t")
            temp[s[4]] = int(s[2])
    return temp

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)