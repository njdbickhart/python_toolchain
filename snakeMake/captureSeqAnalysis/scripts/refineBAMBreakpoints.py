import pysam
import matplotlib as mlp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.signal import find_peaks
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
                        help="Output base file name",
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
    plotDFs = list()
    for c in clusters:
        chrom = c[0].split(':')[0]
        (positions, plotDF) = get_softclipped_bps(args.bam, c[0], args.depth)
        plotDFs.append(plotDF)
        # Finally, merge bases into coordinate range if there are more than one in a row
        # Note: zero regions will print an empty list
        if len(positions) > 1:
            numbers = [int(p[0].split(':')[1]) for p in positions]
            ratios = [float(p[1]) for p in positions]
            depths = [int(p[2]) for p in positions]
            coords = getMinMax(numbers)
            refinement.append((c[0], f'{chrom}:{coords[0]}-{coords[1]}', str(np.mean(ratios)), str(np.mean(depths))))
        elif len(positions) == 1:
            start = int(positions[0][0].split(':')[1])
            end =  start + 1
            refinement.append((c[0], f'{chrom}:{start}-{end}', str(positions[0][1]), str(positions[0][2])))
    
    # THis should eliminate most of the duplicate records
    refinement = set(refinement)    
    with open(args.output + ".tab", 'w') as output:
        output.write('OldCluster\tRefined\tRatioClipped\tTotDepth\n')
        for i, r in enumerate(refinement):
            createDiagnosticPlot(plotDFs[i], r[1], args.output)
            output.write("\t".join(r) + "\n")

    print("Fini")



def createDiagnosticPlot(df, ucsc, outbase):
    # Failsafes to avoid errors printing empty databases
    if df is None:
        return
    if df.empty:
        return
    if len(df.index) == 0:        
        return
    (fig, axis) = plt.subplots(nrows=2, sharex=True, sharey=False)
    ucsc = ucsc.replace(':', '_').replace('-', '_')
    #plt.ticklabel_format(style = 'plain')
    plt.xticks(rotation=15)
    plt.suptitle(f'Insertion evidence for {ucsc}')

    caxis = axis[0]
    sns.scatterplot(data=df, x='Pos', y='Tot', hue='Selected', ax = caxis)
    caxis.get_xaxis().set_major_formatter(mlp.ticker.StrMethodFormatter('{x:,.0f}'))
    caxis.get_yaxis().set_major_formatter(mlp.ticker.StrMethodFormatter('{x:,.0f}'))
    caxis.set_ylabel('Read Depth')

    caxis = axis[1]
    sns.scatterplot(data=df, x='Pos', y='Ratio', hue='Selected', ax = caxis)
    caxis.get_xaxis().set_major_formatter(mlp.ticker.StrMethodFormatter('{x:,.0f}'))
    caxis.get_yaxis().set_major_formatter(mlp.ticker.StrMethodFormatter('{x:,.0f}'))
    caxis.set_ylabel('Read ends / Read Depth')
    ylabels = ['{:,.1f}'.format(x) for x in caxis.get_yticks()]
    caxis.set_yticklabels(ylabels)
    caxis.set_xlabel('Basepair Position')

    fig.tight_layout()
    plt.savefig(f'{outbase}.{ucsc}.png')


def getMinMax(numbers):
    return (min(numbers), max(numbers))

def addToDFList(data, pos, clip, tot, sel):
    data['Pos'].append(pos)
    if tot == 0:
        data['Ratio'].append(0.0)
    else:
        data['Ratio'].append(clip/tot)
    data['Tot'].append(tot)
    data['Selected'].append(sel)
    return data


def get_softclipped_bps(bamfile, region, depthfile, log = True):
    if log:
        print(f'Working on {region} of {bamfile}')
    depth = 0
    with open(depthfile, 'r') as input:
        depth = float(input.readline().rstrip())
    chrsegs = region.split(':')

    positions = list()
    plotData = defaultdict(list)
    depths = list()
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
        depths.append(tot)
        #positions.append((pos, clip, tot))
        if tot == 0:
            plotData = addToDFList(plotData, pos, clip, tot, 'Zero')
        else:
            plotData = addToDFList(plotData, pos, clip, tot, 'NotSelected')
    if len(depths) == 0:
        return ([], None)
    q99 = getQ99Depth(depths) # This is the threshold of depth to determine if a breakpoint is valid
    final = list()
    df = pd.DataFrame(plotData)
    peaks = getPeaks(df, 0.10)

    df.loc[(df['Tot'] >= q99) & (df.index.isin(peaks)), 'Selected'] = 'Selected'
    if log:
        print(df)

    for row in df[df['Selected'] == 'Selected'].iterrows():
        final.append((f'{chrsegs[0]}:{row['Pos']}', row['Ratio'], row['Tot'], row['Selected']))
    # for (pos, clip, tot) in positions:
    #     if tot == 0:
    #         if log:
    #             print(f'{pos}\t{clip}\t{tot}\tZero')
    #             plotData = addToDFList(plotData, pos, clip, tot, 'Zero')
    #         continue
    #     if clip/tot > 0:
    #         if log:
    #             print(f'{pos}\t{clip}\t{tot}\t{clip/tot}')
    #         if clip/tot > 0.20 and tot >= q99:
    #             # We only accept breakpoint coordinates that have > 20% clipped start/end reads and depth of coverage above 99%
    #             print(f'{pos}\t{clip}\t{tot}\t{clip/tot}\tSelected')
    #             final.append((f'{chrsegs[0]}:{pos}', clip/tot, tot))
    #             plotData = addToDFList(plotData, pos, clip, tot, 'Selected')
    #         else:
    #             plotData = addToDFList(plotData, pos, clip, tot, 'NotSelected')
    return (final, df) # first is list of tuples (singlebase coordinate, ratio of clipped bases, total bases), second is a diagnostic dataframe for plotting

def getPeaks(df, prominence):
    nrows = len(df.index)
    minBurn = int(nrows * 0.10)
    maxBurn = int(nrows * 0.90)

    peaks, _ = find_peaks(df['Ratio'], prominence=prominence)

    filtpeaks = set()
    for i in peaks:
        if i > minBurn and i < maxBurn:
            filtpeaks.add(i)
    return filtpeaks

def getQ99Depth(depths):
    nrows = len(depths)
    minBurn = int(nrows * 0.10)
    maxBurn = int(nrows * 0.90)
    return np.quantile(depths[minBurn:maxBurn], 0.99)

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
            chr = p[0].split(':')[0]
            start = int(p[0].split(':')[1])
            ncoords = getMinMax((start, int(end[1])))
            # Filter to remove cluster pairs that are too distant to be real
            # Or are too small to be real
            if abs(ncoords[0] - ncoords[1]) > 100000:
                continue
            elif abs(ncoords[0] - ncoords[1]) < 20:
                continue
            else:
                clusters.append([f'{chr}:{ncoords[0]}-{ncoords[1]}', meandp])
    return sorted(clusters, key=lambda x: x[1], reverse=True)  # returns list of lists: ['chrstring', 'meandepth']

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)