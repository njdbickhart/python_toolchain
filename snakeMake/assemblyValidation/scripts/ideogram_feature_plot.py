import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from itertools import cycle
from collections import defaultdict
import argparse
import pandas
import numpy as np
import pysam

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Plot features on a chromosome ideogram"
            )
    parser.add_argument('-f', '--file',
                        help="The input, frc Features.txt file",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Output png file for plot",
                        type=str, required=True
                        )
    parser.add_argument('-e', '--bed',
                        help="Output bed file of upper quantile regions",
                        type=str, required=True
                        )
    parser.add_argument('-b', '--bam',
                        help="Input indexed and sorted bam file",
                        type=str, required=True
                        )
    parser.add_argument('-t', '--threshold',
                        help="How many of the largest chromosomes to plot? [30]",
                        type=int, default=30
                        )

    return parser.parse_args()

# Code refactored from Ryan Dale's example: https://gist.github.com/daler/c98fc410282d7570efc3
# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
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

    # opening the bam file with pysam
    bamfile = pysam.AlignmentFile(input, 'rb')
    # query all the names of  the chromosomes in a list
    list_chromosomes = bamfile.references
    list_length = bamfile.lengths
    bamfile.close()
    return list_chromosomes, list_length

if __name__ == "__main__":
    args = parse_user_input()

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
    list_chromosomes, list_length = get_chromosomes_names(args.bam)
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

    # ideo = pandas.read_table(
    #     'ideogram.txt',
    #     skiprows=1,
    #     names=['chrom', 'start', 'end', 'name', 'gieStain']
    # )
    #
    # # Filter out chromosomes not in our list
    # ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]
    #
    # # Add a new column for width
    # ideo['width'] = ideo.end - ideo.start
    #
    # # Colors for different chromosome stains
    # color_lookup = {
    #     'gneg': (1., 1., 1.),
    #     'gpos25': (.6, .6, .6),
    #     'gpos50': (.4, .4, .4),
    #     'gpos75': (.2, .2, .2),
    #     'gpos100': (0., 0., 0.),
    #     'acen': (.8, .4, .4),
    #     'gvar': (.8, .8, .8),
    #     'stalk': (.9, .9, .9),
    # }
    #
    # # Add a new column for colors
    # ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


    # Same thing for genes
    gtable = defaultdict(list)
    window = defaultdict(list)
    for c, v in zip(chromosome_list, chromosome_size):
        for i in range(int(v / 500000) + 1):
            window[c].append(0)
    with open(args.file, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            if s[0] not in chromosome_list:
                continue
            index = int(int(s[2]) / 500000)
            window[s[0]][index] += 1

    # calculate quantiles for color mapping
    density = []
    for c, l in window.items():
        density.extend(l)

    qvals = np.quantile(density, (0.0, 0.25, 0.75, 0.95, 1.0))
    print(f'Quartiles: {qvals}')
    qcols = ['#EAF2F8','#EAF2F8', '#E2E73C', '#E74C3C']
    with open(args.bed, 'w') as bed:
        for c, l in window.items():
            for j, v in enumerate(l):
                gtable['chrom'].append(c)
                gtable['start'].append(j * 500000)
                gtable['end'].append((j * 500000) + 500000)
                color = '#EAF2F8'
                for i, q in enumerate(qvals):
                    if i == 0:
                        continue
                    if v <= int(qvals[i]) and v >= int(qvals[i-1]):
                        color = qcols[i-1]
                    if i == len(qvals) - 1:
                        # The last quartile range
                        start = j * 500000
                        end = start + 500000
                        bed.write(f'{c}\t{start}\t{end}\t{v}\n')
                gtable['colors'].append(color)

    genes = pandas.DataFrame(gtable)
    print(genes.head())
    # genes = pandas.read_table(
    #     'ucsc_genes.txt',
    #     names=['chrom', 'start', 'end', 'name'],
    #     usecols=range(4))
    # genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
    # genes['width'] = genes.end - genes.start
    # genes['colors'] = '#2243a8'


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
