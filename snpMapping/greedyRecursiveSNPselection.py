# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 11:01:53 2022

@author: Derek.Bickhart
"""

import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from itertools import cycle
import seaborn as sns
import pandas as pd
import numpy as np

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A recursive SNP selection utility"
            )
    parser.add_argument('-f', '--fai',
                        help="Input reference fasta index file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file base name",
                        required=True, type=str,
                        )
    parser.add_argument('-n', '--number',
                        help="Number of SNPs to select across the entire genome",
                        required=True, type=int,
                        )
    parser.add_argument('-s', '--snp',
                        help="Input SNP data frame file",
                        required=True, type=str,
                        )
    parser.add_argument('-m', '--min',
                        help="Minimum length between SNPs [50 kb]",
                        default = 50000, type=int,
                        )
    parser.add_argument('-e', '--exclude',
                        help="Newline delimited list of markers to exclude",
                        default = "None", type=str
                        )
    return parser.parse_args(), parser

def main(args, parser):
    chrlens = dict()
    total = 0
    # Read in fasta index file and get chromosome lengths
    with open(args.fai, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            chrlens[s[0]] = int(s[1])
            total += int(s[1])

    # create dictionary of snps for selection
    snpObjs = defaultdict(list)
    # Markers to exclude
    exclude = set()
    if args.exclude != "None":
        with open(args.exclude, 'r') as input:
            for l in input:
                l = l.rstrip()
                exclude.add(l)

    with open(args.snp, 'r') as input:
        for i in range(2):
            # Remove two header lines
            next(input)
        for l in input:
            s = l.rstrip().split()
            unmask = 0 if s[0] in exclude else 1
            snpObjs[s[3]].append(SNP(s[3], s[0], int(s[4]), s[5], s[6], float(s[1]), float(s[2]), unmask))

    selected = defaultdict(list)
    bedstyle = list()
    # Now divide chromosomes and identify SNPs based on the proportion of the total length of the chromosome
    with open(args.output + ".bim", 'w') as out:
        for c, s in snpObjs.items():
            prop = chrlens[c] / total
            maxsnps = int(prop * args.number)
            print(f'Working on chromosome {c} of proportional total length {prop} and {maxsnps} out of {args.number}')
            selsnps = greedyRecursive(snpObjs[c], 1, chrlens[c], 1, maxsnps, args.min, set())

            selsnps = list(selsnps)
            selsnps.sort()
            selected[c] = selsnps
            for g in selsnps:
                out.write(g.bimFormat() + '\n')
                bedstyle.append([g.chrom, g.pos, g.pos + 50, g.name])
            print(f'Wrote {len(selsnps)} snps for chrom: {c} out of an expected {maxsnps}')

    # Plot output
    averageDistancePlots(snpObjs, selected, args.output)
    mafDistributionPlot(snpObjs, selected, args.output)

    ideoworker = chrIdeogramPlot(args.fai, 18)

    snptable = ideoworker.straight_segment_df(bedstyle)
    ideoworker.plot_it(snptable, args.output + '.ideogram.pdf')

    print("Fini!")

def SNPSelectionWorkflow(fai, snp, number, excludeFile = "None", minimum = 50000):
    chrlens = dict()
    total = 0
    # Read in fasta index file and get chromosome lengths
    with open(fai, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            chrlens[s[0]] = int(s[1])
            total += int(s[1])

    # create dictionary of snps for selection
    snpObjs = defaultdict(list)
    # exclude list
    exclude = set()
    if excludeFile != "None":
        with open(excludeFile, 'r') as input:
            for l in input:
                l = l.rstrip()
                exclude.add(l)

    with open(snp, 'r') as input:
        for i in range(2):
            # Remove two header lines
            next(input)
        for l in input:
            s = l.rstrip().split()
            unmask = 0 if s[0] in exclude else 1
            snpObjs[s[3]].append(SNP(s[3], s[0], int(s[4]), s[5], s[6], float(s[1]), float(s[2]), unmask))

    selected = defaultdict(list)
    bedstyle = list()

    for c, s in snpObjs.items():
        prop = chrlens[c] / total
        maxsnps = int(prop * number)
        print(f'Working on chromosome {c} of proportional total length {prop} and {maxsnps} out of {number}')
        selsnps = greedyRecursive(snpObjs[c], 1, chrlens[c], 1, maxsnps, minimum, set())

        selsnps = list(selsnps)
        selsnps.sort()
        selected[c] = selsnps
        for g in selsnps:
            bedstyle.append([g.chrom, g.pos, g.pos + 50, g.name])
        print(f'Selected {len(selsnps)} snps for chrom: {c} out of an expected {maxsnps}')

    return (selected, bedstyle)

def averageDistancePlots(snpObjs, selsnps, outbase):
    cols = ["Chr", "Set", "Dist"]
    dists = defaultdict(list)
    for n, j in zip(["Original", "Selected"], [snpObjs, selsnps]):
        positions = defaultdict(list)
        for c in j.keys():
            for i in j[c]:
                positions[c].append(i.pos)

        distances = list()
        for c in positions.keys():
            positions[c].sort()
            for i in range(1, len(positions[c])):
                distances.append(positions[c][i] - positions[c][i-1])
        for d in distances:
            dists["Chr"].append("chr1")
            dists["Set"].append(n)
            dists["Dist"].append(d)

    tempdf = pd.DataFrame(dists, columns=cols)
    tempdf['logNorm'] = np.log10(tempdf['Dist'])

    oqvals = np.quantile(tempdf.loc[tempdf["Set"] == "Original"]['Dist'], (0.0,0.25, 0.50, 0.75, 1.0))
    qvals = np.quantile(tempdf.loc[tempdf["Set"] == "Selected"]['Dist'], (0.0, 0.25, 0.50, 0.75, 1.0))
    print(qvals)

    #print(tempdf.describe)
    ax = plt.subplot()
    sns.violinplot(x="Set", y="Dist", data=tempdf,  ax=ax)

    #ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    #ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("${{{x}}} Mbp$"))
    jitter = 0.1
    for i in qvals:
        ax.text(1 + jitter, i, f'{i / 1000:.0f} kbp')
        if jitter >= 0.3:
            jitter -= 0.5
        else:
            jitter += 0.2

    jitter = 0.1
    for i in oqvals:
        ax.text(0 + jitter, i, f'{i / 1000:.0f} kbp')
        if jitter >= 0.3:
            jitter -= 0.5
        else:
            jitter += 0.2

    plt.ylabel('Distance between markers (log10)')
    plt.savefig(outbase + '.distances.pdf')
    ax.cla()
    plt.close()

def mafDistributionPlot(snpObjs, selsnps, outbase):
    cols = ["Set", "MAF"]
    mafs = defaultdict(list)
    for n, i in zip(["Original", "Selected"], [snpObjs, selsnps]):
        for c in i.keys():
            for s in i[c]:
                mafs["Set"].append(n)
                mafs["MAF"].append(s.maf)

    tempdf = pd.DataFrame(mafs)
    ax = plt.subplot()
    sns.violinplot(data=tempdf, x="Set", y="MAF", ax=ax)
    plt.savefig(outbase + '.mafs.pdf')
    ax.cla()
    plt.close()


def greedyRecursive(snppool, start, end, curselect, maxselect, minlength, selsnps):
    num = len(snppool)
    if num <= 0 or end - start <= minlength or curselect >= maxselect:
   		# termination conditions:
   		# 1. no more SNPs in the region
   		# 2. The end/start region is smaller than the minimum length between snps
   		# 3. We have selected more than enough SNPs for now
   		return selsnps

    for i in snppool:
        i.calcScore(start, end)

    # Run through and select the highest snp; add it to pool if it is beyond the minimum length from previous SNP entries
    loop = True
    selection = None
    while(loop):
        topSNP = None
        topScore = 0
        for s in snppool:
            if s.curscore > topScore:
                topSNP = s
                topScore = s.curscore

        if topScore == 0:
            # We ran out of SNPs and could not bridge the min distance
            return selsnps
        elif topSNP.pos - start < minlength or end - topSNP.pos < minlength:
            topSNP.curscore = 0
            topScore = 0
        else:
            selection = topSNP
            break

    selsnps.add(selection)

    # Create two pools of SNPs, before and after the selected SNP and recurse
    spool = []
    for s in snppool:
        if s.pos < selection.pos:
            spool.append(s)

    epool = []
    for s in snppool:
        if s.pos > selection.pos:
            epool.append(s)

    #TODO: decrease penalty for selecting markers closer to the ends of the chromosomes.
    #TODO: find out if the termini of the chromosomes should be defined at 1 mb or 2 mb or other value.
    selsnps = greedyRecursive(spool, start, selection.pos, curselect * 2, maxselect, minlength, selsnps)
    selsnps = greedyRecursive(epool, selection.pos, end , curselect * 2, maxselect, minlength, selsnps)

    return selsnps

class chrIdeogramPlot:

    def __init__(self, fai, threshold):
        self.fai = fai
        self.threshold = threshold

        self.colcycle = cycle([ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
                  '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
                  '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ])

        # Decide which chromosomes to use
        self.list_chromosomes, list_length = self.get_chromosomes_names(self.fai)
        self.chr_assoc = {k : v for k, v in zip(self.list_chromosomes, list_length)}
        self.list_chromosomes = [x for x, j in sorted(self.chr_assoc.items(), key=lambda item: item[1], reverse=True)]
        self.chromosome_list = self.list_chromosomes[:self.threshold]
        list_length = [v for k, v in sorted(self.chr_assoc.items(), key=lambda item: item[1], reverse=True)]
        self.chromosome_size = list_length[:self.threshold]


    def plot_it(self, track, output):
        # Height of each ideogram
        chrom_height = 0.1

        # Spacing between consecutive ideograms
        chrom_spacing = 0.6

        # Height of the gene track. Should be smaller than `chrom_spacing` in order to
        # fit correctly
        track_height = 0.4

        # Padding between the top of a gene track and its corresponding ideogram
        track_padding = 0.1

        # Width, height (in inches)
        figsize = (6, 8)


        # Keep track of the y positions for ideograms and genes for each chromosome,
        # and the center of each ideogram (which is where we'll put the ytick labels)
        ybase = 0
        chrom_ybase = {}
        track_ybase = {}
        chrom_centers = {}

        # Iterate in reverse so that items in the beginning of `chromosome_list` will
        # appear at the top of the plot
        for chrom in self.chromosome_list[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + chrom_height / 2.
            track_ybase[chrom] = ybase - track_height - track_padding
            ybase += chrom_height + chrom_spacing

        # Read in ideogram.txt, downloaded from UCSC Table Browser

        ideo = pd.DataFrame({'chrom' : self.chromosome_list,
                                'start' : [0 for x in range(len(self.chromosome_list))],
                                'width' : self.chromosome_size,
                                'colors': ['#bd2309' for x in range(len(self.chromosome_list))]})

        # Note, I am plotting exact intervals here instead of windows
        # If the plots are too washed out, I may have to return to window-based clustering


        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        # Now all we have to do is call our function for the ideogram data...
        print("adding ideograms...")
        for collection in self.chromosome_collections(ideo, chrom_ybase, chrom_height):
            ax.add_collection(collection)

        # ...and the gene data
        print("adding tracks...")
        for collection in self.chromosome_collections(
            track, track_ybase, track_height, alpha=0.5, linewidths=0
        ):
            ax.add_collection(collection)

        # Axes tweaking
        ax.set_yticks([chrom_centers[i] for i in self.chromosome_list])
        ax.set_yticklabels(self.chromosome_list)
        ax.axis('tight')
        plt.savefig(output)
        plt.close()

    def heatmap_segment_df(self, bedsegs, winsize = 500000, colors=None):
        gtable = defaultdict(list)
        window = defaultdict(list)
        for c, v in zip(self.chromosome_list, self.chromosome_size):
            for i in range(int(v / winsize) + 1):
                window[c].append(0)

        for s in bedsegs:
            if s[0] not in self.chromosome_list:
                continue
            index = int(int(s[2]) / winsize)
            window[s[0]][index] += 1

        # calculate quantiles for color mapping
        density = []
        for c, l in window.items():
            density.extend(l)

        qvals = np.quantile(density, (0.0, 0.25, 0.75, 0.95, 1.0))
        print(f'Quartiles: {qvals}')
        qcols = ['#EAF2F8','#EAF2F8', '#E2E73C', '#E74C3C']

        ## TODO rewrite!
        with open(args.bed, 'w') as bed:
            for c, l in window.items():
                for j, v in enumerate(l):
                    gtable['chrom'].append(c)
                    gtable['start'].append(j * winsize)
                    gtable['end'].append((j * winsize) + winsize)
                    color = '#EAF2F8'
                    for i, q in enumerate(qvals):
                        if i == 0:
                            continue
                        if v <= int(qvals[i]) and v >= int(qvals[i-1]):
                            color = qcols[i-1]
                        if i == len(qvals) - 1:
                            # The last quartile range
                            start = j * winsize
                            end = start + winsize
                            bed.write(f'{c}\t{start}\t{end}\t{v}\n')
                    gtable['colors'].append(color)
        return gtable

    def straight_segment_df(self, bedsegs, colors = None):
        gtable = defaultdict(list)
        if colors == None:
            colors = self.colcycle
        else:
            colors = cycle(colors)
        names = set()

        for s in bedsegs:
            if s[0] not in self.chromosome_list:
                continue
            gtable['chrom'].append(s[0])
            gtable['start'].append(int(s[1]))
            gtable['end'].append(int(s[2]))
            gtable['name'].append(s[3])
            names.add(s[3])

        ctranslate = dict()
        for i in names:
            ctranslate[i] = next(colors)

        genes = pd.DataFrame(gtable)
        genes['colors'] = genes['name'].map(ctranslate)
        return genes

    def chromosome_collections(self, df, y_positions, height,  **kwargs):
        """
        Yields BrokenBarHCollection of features that can be added to an Axes
        object.
        Parameters
        ----------
        df : pandas.DataFrame
            Must at least have columns ['chrom', 'start', 'end', 'color']. If no
            column 'width', it will be calculated from start/end.
        y_positions : dict
            Keys are chromosomes, values are y-value at which to anchor the
            BrokenBarHCollection
        height : float
            Height of each BrokenBarHCollection
        Additional kwargs are passed to BrokenBarHCollection
        """
        del_width = False
        if 'width' not in df.columns:
            del_width = True
            df['width'] = df['end'] - df['start']
        for chrom, group in df.groupby('chrom'):
            print(chrom)
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield BrokenBarHCollection(
                xranges, yrange, facecolors=group['colors'], **kwargs)
        if del_width:
            del df['width']

    def get_chromosomes_names(self, input):
        list_chromosomes = []
        list_length = []
        with open(input, 'r') as fai:
            for l in fai:
                segs = l.rstrip().split()
                list_chromosomes.append(segs[0])
                list_length.append(int(segs[1]))
        return list_chromosomes, list_length


class SNP:

    def __init__(self, chrom, name, pos, a1, a2, maf, miss, unmask):
        self.chrom = chrom
        self.name = name
        self.pos = pos
        self.a1 = a1
        self.a2 = a2
        self.maf = maf
        self.call = 1 - miss
        self.unmask = unmask  # Value is either 1 (use in selection) or 0 (do not use)

        self.curscore = 0

    def bimFormat(self):
        return "\t".join([self.chrom, self.name, "0", str(self.pos), str(self.a1), str(self.a2)])

    def calcScore(self, start, end):
        length = end - start
        # TODO: set hard limits on MAF to exclude values that are above and below 0.5 and 0.1
        # TODO: Set a local weighted MAF selection algorithm
        self.curscore = self.unmask * self.maf * self.call * (length - abs(2 * self.pos - (end + start)))

    def __repr__(self):
        return f'SNP({self.chrom}, {self.name}, {self.pos})'

    def __eq__(self, other):
        if isinstance(other, SNP):
            return ((self.chrom == other.chrom) and (self.name == other.name) and (self.pos == other.pos))
        else:
            return False

    def __ne__(self, other):
        return (not self.__eq__(other))

    def __hash__(self):
        return hash(self.__repr__())

    def __lt__(self, other):
        return self.pos < other.pos

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
