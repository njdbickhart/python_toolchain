# These functions are meant only to be called in an ipython notebook in conjunction with the
# "greedyRecursiveSNPselection" script. These are helper functions for data organziation and
# visualization
from upsetplot import UpSet, from_contents
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict

def markerUpsetPlot(selsnplist, selsnpcats):
    gtypesets = dict()

    for cat, l in zip(selsnpcats, selsnplist):
        gtypesets[cat] = list()
        for chrom, snpl in l.items():
            for s in snpl:
                gtypesets[cat].append(s.name)

    UpSet(from_contents(gtypesets), subset_size='count', sort_by='cardinality', show_percentages=True).plot()

def flattenSNPlist(selsnplist, target):
    flist = list()
    for chrom, snpl in selsnplist[target].items():
        for s in snpl:
            flist.append(s)
    return flist;

def mergeMarkers(selsnplist):
    seen = set()
    flist = list()
    overlaps = 0
    for l in selsnplist:
        for chrom, snpl in l.items():
            for s in snpl:
                if s.name not in seen:
                    seen.add(s.name)
                    flist.append(s)
                else:
                    overlaps += 1;
    return (flist, overlaps)

def NBmafDistributionPlot(selsnps, categories):
    cols = ["SNP", "Set", "MAF"]
    mafs = defaultdict(list)
    for n, i in zip(categories, selsnps):
        for c in i.keys():
            for s in i[c]:
                mafs["SNP"].append(s.name)
                mafs["Set"].append(n)
                mafs["MAF"].append(s.maf)

    tempdf = pd.DataFrame(mafs)
    pdf = tempdf.pivot(index= 'SNP', columns='Set', values='MAF')

    values = list()
    minimum = list()
    maximum = list()
    for i in pdf.columns:
        values.append(np.nanmedian(pdf[[i]]))
        maximum.append(float(np.max(pdf[[i]], axis=0)))
        minimum.append(float(np.min(pdf[[i]], axis=0)))

    ax = plt.subplot()
    sns.violinplot(data=tempdf, x="Set", y="MAF", ax=ax)
    for i, v in enumerate(values):
        plt.gca().text(i + 0.2, v, f'{v:.1f}')
        plt.gca().hlines(y=v, xmin=i, xmax=i + 0.2, linestyles = 'dashed', linewidth =1, color='k')
        plt.gca().text(i + 0.2, maximum[i], f'{maximum[i]:.1f}')
        plt.gca().hlines(y=maximum[i], xmin=i, xmax=i + 0.2, linestyles = 'dashed', linewidth =1, color='k')
        plt.gca().text(i + 0.2, minimum[i], f'{minimum[i]:.1f}')
        plt.gca().hlines(y=minimum[i], xmin=i, xmax=i + 0.2, linestyles = 'dashed', linewidth =1, color='k')

def NBaverageDistancePlots(selsnps, categories):
    cols = ["Chr", "Set", "Dist"]
    dists = defaultdict(list)
    for n, i in zip(categories, selsnps):
        positions = defaultdict(list)
        for c in i.keys():
            for x in i[c]:
                positions[c].append(x.pos)

        distances = list()
        for c in positions.keys():
            positions[c].sort()
            for x in range(1, len(positions[c])):
                distances.append(positions[c][x] - positions[c][x-1])
        for d in distances:
            dists["Chr"].append("chr1")
            dists["Set"].append(n)
            dists["Dist"].append(d)

    tempdf = pd.DataFrame(dists, columns=cols)
    tempdf['logNorm'] = np.log10(tempdf['Dist'])

    #oqvals = np.quantile(tempdf.loc[tempdf["Set"] == "Original"]['Dist'], (0.0,0.25, 0.50, 0.75, 1.0))
    qvals = dict()
    for c in categories:
        qvals[c] = np.quantile(tempdf.loc[tempdf["Set"] == c]['Dist'], (0.0, 0.25, 0.50, 0.75, 1.0))
        print(f'{c}: ', qvals)

    #print(tempdf.describe)
    ax = plt.subplot()
    sns.violinplot(x="Set", y="Dist", data=tempdf,  ax=ax)

    #ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    #ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("${{{x}}} Mbp$"))
    jitter = 0.5
    for x in range(len(categories)):
        for i in range(0, 5, 2):
            v = qvals[categories[x]][i]
            ax.text(x + jitter, v, f'{v / 1000:.0f} kbp')
            ax.hlines(y=v, xmin=x, xmax=x + jitter, linestyles = 'dashed', linewidth =1, color='k')
            if jitter == 0.5:
                jitter = -0.5
            else:
                jitter = 0.5


    plt.ylabel('Distance between markers (log10)')

def add_value_labels(ax, spacing=5):
    """Add labels to the end of each bar in a bar chart.

    Arguments:
        ax (matplotlib.axes.Axes): The matplotlib object containing the axes
            of the plot to annotate.
        spacing (int): The distance between the labels and the bars.
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        label = "{:.0f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            rotation = 45.0,
            va=va)

def NBplotSelCounts(selected):
    data = defaultdict(list)
    for c, s in selected.items():
        data[c] = len(s)

    plt.bar(list(data.keys()), data.values(), 1.0, color=(0.1, 0.1, 0.1, 0.1),  edgecolor='blue')
    add_value_labels(plt.gca(), 3)

def NBwriteOutList(selected, outfile):
    with open(outfile, 'w') as output:
        for s in selected:
            output.write(s.bimFormat() + "\n")

def NBwriteOutput(selected, outfile):
    with open(outfile, 'w') as output:
        for c, s in selected.items():
            for g in s:
                output.write(g.bimFormat() + '\n')

def writeOutTempBIM(tdf, colname, filt, mapfile, outfilename):
    conv = {k : v for k, v in zip(tdf['SNP'], tdf['CHR'])}
    bpconv = dict()
    a1conv = dict()
    a2conv = dict()
    with open(mapfile, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            bpconv[s[1]] = int(s[3])
            a1conv[s[1]] = 1
            a2conv[s[1]] = 2

    tDf = tdf[tdf[colname] == filt][['SNP', 'MAF', 'MISS', 'CHR']]
    tDf['Pos'] = tDf['SNP'].apply(lambda x : bpconv.get(x))
    tDf['A1'] = tDf['SNP'].apply(lambda x : a1conv.get(x))
    tDf['A2'] = tDf['SNP'].apply(lambda x : a2conv.get(x))
    tDf.sort_values(['CHR', 'SNP', 'Pos']).to_csv(outfilename, sep="\t", index=False)
    print(f'Wrote to: {outfilename}')

def selectedCountsChrs(selsnplist, selsnpcats):
    counters = defaultdict(list)
    for cat, l in zip(selsnpcats, selsnplist):
        for chrom, snpl in l.items():
            counters["Chr"].append(chrom)
            counters["Type"].append(cat)
            counters["Value"].append(len(snpl))

    tdf = pd.DataFrame(counters)
    sns.barplot(data=tdf, x="Chr", y="Value", hue="Type",
                palette=sns.color_palette("Set2"), order=[str(x) for x in range(1, 34)])
    plt.ylabel('Count of markers selected per chromosome')
    plt.xlabel("Chromosome")
    plt.xticks(rotation=45, ha='right')
