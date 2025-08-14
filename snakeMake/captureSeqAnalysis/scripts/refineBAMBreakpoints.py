import os
import sys
import pysam
import numpy as np
from collections import defaultdict

import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to identify candidate insertion points in a bam file from a list of bedpe locations"
            )
    parser.add_argument('-f', '--file', 
                        help="Input cluster file (bedpe)",
                        required=True, type=str
                        )
    parser.add_argument('-b', '--bam',
                        help="BAM file for pileups",
                        required=True, type=str,
                        )
    parser.add_argument('-o', '--output',
                        help="Output tab file",
                        required=True, type=str,
                        )
    parser.add_argument('-d', '--depth',
                        help="Calibrator depth estimate for the bam", 
                        required=True, type=str,
                        )

    return parser.parse_args(), parser


def main(args, parser):
    # First generate clusters
    clusters = generate_clusters(args.depth, args.file)

    # Next, loop through clusters to get insertion regions
    refinement = list()
    for c in clusters:
        chrom = c[0].split(':')[0]
        positions = get_softclipped_bps(args.bam, c[0], args.depth)
        
        # Finally, merge bases into coordinate range if there are more than one in a row
        # Note: zero regions will print an empty list
        if len(positions) > 1:
            numbers = [int(p[0].split(':')[1]) for p in positions]
            ratios = [float(p[1]) for p in positions]
            depths = [int(p[2]) for p in positions]
            coords = getMinMax(numbers)
            refinement.append((c[0], f'{chrom}:{coords[0]}-{coords[1]}', str(np.mean(ratios)), str(np.mean(depths))))
        elif len(positions) == 1:
            start = int(positions[0].split(':')[1])
            end =  start + 1
            refinement.append((c[0], f'{chrom}:{start}-{end}', str(positions[1]), str(positions[2])))
        
    with open(args.output, 'w') as output:
        output.write('OldCluster\tRefined\tRatioClipped\tTotDepth\n')
        for r in refinement:
            output.write("\t".join(r) + "\n")

    print("Fini")



def getMinMax(numbers):
    return (min(numbers), max(numbers))


def get_softclipped_bps(bamfile, region, depthfile, log = True):
    depth = 0
    with open(depthfile, 'r') as input:
        depth = float(input.readline().rstrip())
    chrsegs = region.split(':')

    positions = list()
    bamreader = pysam.AlignmentFile(bamfile, 'rb')
    for pileup in bamreader.pileup(region=region, stepper="nofilter"):
        pos = pileup.reference_pos
        tot = 0
        clip = 0
        for i in pileup.get_query_sequences(mark_matches=True, mark_ends=True):
            for bp in i:
                if '^' in bp or '$' in bp:
                    clip += 1
                tot += 1
        if tot == 0:
            if log:
                print(f'{pos}\t{clip}\t{tot}\tZero')
            continue
        if clip/tot > 0:
            if log:
                print(f'{pos}\t{clip}\t{tot}\t{clip/tot}')
            if clip/tot > 0.25 and tot > depth / 2:
                print(f'{pos}\t{clip}\t{tot}\t{clip/tot}\tSelected')
                positions.append((f'{chrsegs[0]}:{pos}', clip/tot, tot))
    return positions # list of tuples (singlebase coordinate, ratio of clipped bases, total bases)


def generate_clusters(depthfile, clusterfile):
     # Depth value
    dpvalue = ""
    with open(depthfile, 'r') as input:
        dpvalue = input.readline().rstrip()
    
    # Nodes
    clusters = list()
    with open(clusterfile, 'r') as input:
        head = input.readline()
        lines = input.readlines()
        for i in range(0, len(lines), 2):
            p = lines[i].rstrip().split()
            c = lines[i+1].rstrip().split()
            end = c[1].split(':')
            meandp = int(p[2]) + int(c[2]) / 2
            # Filter to remove integration sites with weak evidence
            if meandp < float(dpvalue) / 5:
                continue
            # Filter to remove cluster pairs that are too distant to be real
            #if abs(int(p[0].split(':')[1]) - int(end[1])) > 100000:
            #    continue
            #else:
            clusters.append([f'{p[0]}-{end[1]}', meandp])
    return sorted(clusters, key=lambda x: x[1], reverse=True)  # returns list of lists: ['chrstring', 'meandepth']

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)