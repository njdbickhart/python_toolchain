import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to generate bed files for wally screenshots"
            )
    parser.add_argument('-f', '--file', 
                        help="Input cluster file",
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
        lines = input.readlines()
        for i in range(0, len(lines), 2):
            p = lines[i].rstrip().split()
            c = lines[i+1].rstrip().split()
            end = c[1].split(':')
            start = p[0].split(':')
            chr = start[0]
            truestart = min(int(start[1]), int(end[1])) - args.extension
            trueend = max(int(start[1]), int(end[1])) + args.extension
            meandp = int(p[2]) + int(c[2]) / 2
            # Filter to remove integration sites with weak evidence
            if meandp < float(dpvalue) / 5:
                continue
            if int(start[1]) == int(end[1]):
                continue # skip areas that appear to be artifacts
            else:
                clusters.append([f'{chr}\t{truestart}\t{trueend}', meandp])
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