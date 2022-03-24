# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:08:35 2022

@author: derek.bickhart-adm
"""

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from itertools import cycle
from collections import defaultdict
import argparse
import pandas
import numpy as np

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A wrapper script for the FlashFry CRISPR design tool"
            )
    parser.add_argument('-f', '--fai', 
                        help="Input FAI file with the names and lengths of chromosomes",
                        required=True, type=str
                        )
    parser.add_argument('-b', '--bed', 
                        help="Input bed file with names to color code",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name (.pdf or .png)",
                        required=True, type=str,
                        )
    parser.add_argument('-t', '--threshold',
                        help="How many of the largest chromosomes to plot? [30]",
                        type=int, default=30
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Height of each ideogram
    chrom_height = 1

    # Spacing between consecutive ideograms
    chrom_spacing = 1

    # Height of the gene track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    gene_height = 0.4

    # Padding between the top of a gene track and its corresponding ideogram
    gene_padding = 0.1

    # Width, height (in inches)
    figsize = (6, 8)

    # Decide which chromosomes to use
    list_chromosomes, list_length = get_chromosomes_names(args.fai)
    chr_assoc = {k : v for k, v in zip(list_chromosomes, list_length)}
    list_chromosomes = [x for x, j in sorted(chr_assoc.items(), key=lambda item: item[1], reverse=True)]
    chromosome_list = list_chromosomes[:args.threshold]
    list_length = [v for k, v in sorted(chr_assoc.items(), key=lambda item: item[1], reverse=True)]
    chromosome_size = list_length[:args.threshold]
    
    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}
    
    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # Read in ideogram.txt, downloaded from UCSC Table Browser
    colcycle = cycle([ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
              '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
              '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ])
    ideo = pandas.DataFrame({'chrom' : chromosome_list,
                            'start' : [0 for x in range(len(chromosome_list))],
                            'width' : chromosome_size,
                            'colors': [next(colcycle) for x in range(len(chromosome_list))]})
    
    # Note, I am plotting exact intervals here instead of windows
    # If the plots are too washed out, I may have to return to window-based clustering
    genes = straight_plot_df(chromosome_list, chromosome_size, args.bed, ['#EAF2F8', '#E2E73C', '#005295', '#2B0C1C'])
    
    print(genes.head())
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    # ...and the gene data
    print("adding genes...")
    for collection in chromosome_collections(
        genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
    ):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')
    plt.savefig(args.output)
    

def straight_plot_df(chromosome_list, chromosome_size, bed, colors):
    gtable = defaultdict(list)
    cycler = cycle(colors)
    names = set()
    with open(bed, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            gtable['chrom'].append(s[0])
            gtable['start'].append(int(s[1]))
            gtable['end'].append(int(s[2]))
            gtable['name'].append(s[3])
            names.add(s[3])
            
    ctranslate = dict()
    for i in names:
        ctranslate[i] = next(cycler)
        
    genes = pandas.DataFrame(gtable)
    genes['color'] = genes['name'].map(ctranslate)
    return genes
    
def chromosome_collections(df, y_positions, height,  **kwargs):
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

def get_chromosomes_names(input):
    list_chromosomes = []
    list_length = []
    with open(input, 'r') as fai:
        for l in fai:
            segs = l.rstrip().split()
            list_chromosomes.append(segs[0])
            list_length.append(int(segs[1]))
    return list_chromosomes, list_length

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
