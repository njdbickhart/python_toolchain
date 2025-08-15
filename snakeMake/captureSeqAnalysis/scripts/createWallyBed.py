import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import argparse

ucsc = re.compile(r'(.+):(\d+)-(\d+)')

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to generate bed files for wally screenshots"
            )
    parser.add_argument('-f', '--file', 
                        help="Input refinement file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output bed file name",
                        required=True, type=str,
                        )
    parser.add_argument('-d', '--depth',
                        help="Calibrator depth files",
                        required=True, type=str
                        )
    parser.add_argument('-s', '--samples', 
                         help="Sample name", 
                         required=True, type=str
                         )
    parser.add_argument('-e', '--extension', 
                        help="number of bases to extend the BED file for the snapshot", 
                        required=True, type=int
                        )
    return parser.parse_args(), parser


def main(args, parser):
    # Depth value
    dpvalue = ""
    with open(args.depth, 'r') as input:
        dpvalue = input.readline().rstrip()

    samp = args.samples

    # Nodes
    clusters = list()
    with open(args.file, 'r') as input:
        head = input.readline()
        for l in input:
            s = l.rstrip().split()
            m = re.search(ucsc, s[1])
            (chr, start, end) = m.group(1,2,3)
            #chr = s[1].split(':')[0]
            #interval = s[1].split(':')[1]
            start = int(start) - args.extension
            end = int(end) + args.extension
            clusters.append([f'{chr}\t{start}\t{end}', int(s[3])])
       
    clusters = sorted(clusters, key=lambda x: x[1], reverse=True)

    counter = len(clusters)
    with open(args.output, 'w') as bed:
        for (cluster, depth) in clusters:
            csegs = cluster.split()
            # Output file is in a folder of the sample_name
            outname = f"snapshots/{samp}/{csegs[0]}_{csegs[1]}_{csegs[2]}_{depth}cov"
            bed.write(f'{csegs[0]}\t{csegs[1]}\t{csegs[2]}\t{outname}\n')
    
    print(f'Successfully created {args.output} with {counter} entries')        



if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)