# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:14:28 2020
@author: derek.bickhart-adm
"""

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib import cm
from itertools import cycle
from collections import defaultdict
import argparse
import pandas
import numpy as np
import pysam
import seaborn as sns
import logging as log

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to plot average bed file depth of coverage within a series of non-sliding windows"
            )
    parser.add_argument('-f', '--fai',
                        help="Input reference fasta index file for the bin",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are {output}.wins and {output}.pdf",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bam',
                        help="Input read depth bam file (may be plot multiple times)",
                        action="append", default=[]
                        )
    parser.add_argument('-u', '--chrom',
                        help="List of chromosomes to plot together",
                        action="append", default=[]
                        )
    parser.add_argument('-i', '--binsize',
                        help="Bin size in bases [5000 bp]",
                        type = int, default=5000
                        )
    return parser.parse_args(), parser

def main(args, parser):
    # Get the list of chromosomes to plot
    log.basicConfig(level=log.INFO)

    plotAll = True
    chrToPlot = set()
    if len(args.chrom) > 0:
        plotAll = False
        for f in args.chrom:
            chrToPlot.add(f)
        log.info(f'Plotting {len(chrToPlot)} chromosomes in the same plot')
    else:
        log.info(f'Attempting to plot entire bam file')


    # Get the contig length list
    ctglens = dict()
    with open(args.fai, 'r') as fai:
        for l in fai:
            s = l.rstrip().split()
            ctglens[s[0]] = int(s[1])
        if plotAll:
            chrToPlot = set([x for x in ctglens.keys()])

    # Create windows
    winlist = defaultdict(list)
    # offset bp to add for stitching contigs together in one line
    ctgoffset = dict()
    breaks = list()
    lastbp = 0
    for c in chrToPlot:
        c = c.rstrip()
        ctgoffset[c] = lastbp + 100
        for i in range(0, ctglens[c], args.binsize):
            winlist[c].append(window(c, i, i + args.binsize))
        lastbp += ctglens[c]
        breaks.append(lastbp)

    # Iterate through each bam file to generate the count of reads in each window
    bamToIdx = dict()
    # Max count value to plot contig text later
    maxy = 0
    for bamidx, bamfile in enumerate(args.bam):
        log.info(f'Reading from bamfile: {bamfile}')
        bamToIdx[bamidx] = bamfile
        # read each sam region and count the reads
        with pysam.AlignmentFile(bamfile, 'rb') as bamhandle:
            for c, w in winlist.items():
                for i, win in enumerate(w):
                    count = 0
                    try:
                        for s in bamhandle.fetch(c, win.start, win.end):
                            if s.is_secondary:
                                continue
                            count += 1
                    except:
                        log.info(f'Could not find {c} in bam file: {bamfile}!')
                        continue
                    if count > maxy:
                        maxy = count
                    winlist = updateWin(winlist, c, i, count, bamidx)

    # OK, data is in! Let's try plotting
    raw = defaultdict(list)
    bars = list()
    for c, w in winlist.items():
        bars.append([ctgoffset[c], ctglens[c], c])
        for win in w:
            for i, h in enumerate(sorted(win.count)):
                raw["contig"].append(c)
                raw["realSeg"].append(f'{win.start}-{win.end}')
                raw["start"].append(win.start + ctgoffset[c])
                raw["end"].append(win.end + ctgoffset[c])
                raw["bam"].append(bamToIdx[i])
                raw["count"].append(h)

    df = pandas.DataFrame(raw)
    df.to_csv(args.output + '.wins', sep='\t', header=True)

    log.info(f'Dataframe complete, head: {df.head()}')

    fig, ax = plt.subplots()
    sns.lineplot(data=df, x='start', y='count', ax=ax, palette='dark')
    plt.xticks(rotation=45)

    # Now to add contig break points
    for b in bars:
        midpoint = b[0] - int(b[1] / 2)
        if midpoint < 0:
            midpoint = int(b[1] / 2)
        ax.axvline(b[0], ls='-', c='k', zorder=-1)
        ax.text(midpoint, maxy * 0.80, b[2], rotation = 'vertical', verticalalignment='top', size=4.0)

    log.info(f'Plot saved to: {args.output}.pdf')
    plt.savefig(args.output + '.pdf')

def updateWin(winlist, contig, winidx, count, bamidx):
    winlist[contig][winidx].addCount(bamidx, count)
    return winlist


class window:

    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start
        self.end = end
        self.count = list()

    def addCount(self, idx, count):
        if idx < len(self.count):
            self.count[idx] = count
        else:
            self.count.append(count)

    def getCount(self, idx):
        if idx in self.count:
            return self.count[idx]
        else:
            return 0


if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
