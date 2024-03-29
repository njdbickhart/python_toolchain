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
                        help="Input Kaiju count file. Can be specified more than once. Must be accompanied by an '-s' tag. Preferred over '-r' option",
                        action="append", default=[]
                        )
    parser.add_argument('-r', '--read', 
                        help="Input Kaiju read file. Can be specified more than once. Must be accompanied by an '-s' tag",
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
    # Check the argument requirements
    proceed = False
    total_files = len(args.file) + len(args.read)
    if total_files == len(args.sample) and total_files > 0:
        proceed = True
    
    if not proceed:
        print("Error! Please specify equal numbers of -f and -s in your input files!")
        print(parser.print_usage())
        sys.exit(-1)
        
    container = defaultdict(dict)
    taxset = set()
    if len(args.file) > 1:
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
    else:
        for f, s in zip(args.read, args.sample):
            container[s] = count_reads(f)
            taxset.update(list(container[s].keys()))
            
        with open(args.output, 'w') as output:
            output.write("sk\tcl\tp\tc\to\tf\tg\tsg\ts\t" + "\t".join(args.sample) + "\n")
            for t in taxset:
                tsegs = taxset.split(";")
                output.write("\t".join(tsegs))
                for s in args.sample:
                    count = container[s].get(t, 0)
                    output.write(f'\t{count}')
                output.write('\n')
    print("Fini")

def create_tax_dict(file):
    temp = dict()
    with open(file, 'r') as input:
        h = input.readline()
        for l in input:
            s = l.rstrip().split("\t")
            temp[s[4]] = int(s[2])
    return temp

def count_reads(file):
    temp = defaultdict(int)
    with open(file, 'r') as input:
        for l in input:
            s = l.rstrip().split("\t")
            if s[0] == 'U':
                temp['Unclassified;;;;;;;;'] += 1
            else:
                bsegs = s[3].split('; ')
                tstr = ''
                for i in range(9):
                    if i + 2 > len(bsegs):
                        tstr += ';'
                    else:
                        tstr += bsegs[i+1] + ';'
                if tstr[-1] == ';':
                    tstr = tstr[:-2]
                temp[tstr] += 1
    return temp

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)