import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to consolidate capture seq information for gene editing"
            )
    parser.add_argument('-f', '--file', 
                        help="Input refinement file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name. Tab delimited",
                        required=True, type=str,
                        )
    parser.add_argument('-d', '--depth',
                        help="Calibrator depth files",
                        required=True, type=str
                        )
    parser.add_argument('-c', '--contaminant',
                        help="Contaminant files",
                        required=True, type=str
                        )
    parser.add_argument('-s', '--samples', 
                         help="Sample names", 
                         required=True, action='append'
                         )
    parser.add_argument('-g', '--genes', 
                         help="Gene intersections", 
                         required=True, action='append'
                         )
    return parser.parse_args(), parser


def main(args, parser):
    data = defaultdict(list)

    dpfiles = list()
    cpfiles = list()
    contfiles = list()
    genefiles = list()
    with open(args.file, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            cpfiles.extend(s)

    with open(args.depth, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            dpfiles.extend(s)

    with open(args.contaminant, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            contfiles.extend(s)

    with open(args.genes, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            genefiles.extend(s)


    for samp, file, depth, contaminants, genes in zip(args.samples, cpfiles, 
                                                      dpfiles, contfiles, genefiles):
        # Depth value
        dpvalue = ""
        with open(depth, 'r') as input:
            dpvalue = input.readline().rstrip()
        
        # Nodes
        clusters = list()
        with open(file, 'r') as input:
            head = input.readline()
            for l in input:
                s = l.rstrip().split()
                clusters.append([s[1], float(s[3])])
            # lines = input.readlines()
            # for i in range(0, len(lines), 2):
            #     p = lines[i].rstrip().split()
            #     c = lines[i+1].rstrip().split()
            #     end = c[1].split(':')
            #     meandp = int(p[2]) + int(c[2]) / 2
            #     # Filter to remove integration sites with weak evidence
            #     if meandp < float(dpvalue) / 5:
            #         continue
            #     else:
            #         clusters.append([f'{p[0]}-{end[1]}', meandp])
        clusters = sorted(clusters, key=lambda x: x[1], reverse=True)

        cstr = list()
        dstr = list()
        chrs = set()
        for c in clusters:
            chrsegs = c[0].split(':')
            chrs.add(chrsegs[0])
            cstr.append(c[0])
            dstr.append(str(c[1]))        

        # Contaminants
        cont = list()
        with open(contaminants, 'r') as input:
            for l in input:
                cont.append(l.rstrip())

        # Gene intersections
        totcount = 0
        foundcount = 0
        overlaps = dict()
        sortedOverlaps = list()
        with open(genes, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                overlaps[s[0]] = s[1]
                totcount += 1
            for c in clusters:
                if c[0] in overlaps:
                    foundcount += 1
                sortedOverlaps.append(overlaps.get(c[0], ""))

        # Logging print statement to check dictionary lookup of cluster
        print(f'{samp} gene overlaps {totcount} and found: {foundcount}')

        data['SAMPLE'].append(samp)
        data['CHROMOSOME'].append(";".join(list(chrs)))
        data['SITES'].append(";".join(cstr))
        data['DEPTHS'].append(";".join(dstr))
        data['CAL_DEPTH'].append(dpvalue)
        data['CONTAMINANTS'].append(";".join(cont))
        data['GENE_OVERLAPS'].append(";".join(sortedOverlaps))

    # Creating the dataframe and printing the consolidated table
    df = pd.DataFrame(data)

    df.to_csv(args.output, index=False)
            

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)