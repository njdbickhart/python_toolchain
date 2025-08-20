import os
import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import pandas as pd
import numpy as np

usage = f'python3 {sys.argv[0]} <input tab summary> <samples per page> <output pdf filename>'

if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)

def number_formatter(number, pos=None):
    """Convert a number into a human readable format."""
    magnitude = 0
    while abs(number) >= 1000:
        magnitude += 1
        number /= 1000.0
    return '%.0f%s' % (number, ['', 'K', 'M', 'B', 'T', 'Q'][magnitude])

def space_labels(ax, nsample):
    displayidx = int(nsample/20)
    for i, l in enumerate(ax.xaxis.get_ticklabels()):
        if i % displayidx != 0:
            l.set_visible(False)

def cov_thresh(ax, df):
    ax.set_title("Count of bases under coverage threshold", loc="left")

    ax.set_ylabel("Bases")
    ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')

    #df[["SNAME", "Zbp"]].plot(kind='bar', x="SNAME", y="Zbp", ax=ax)
    sns.barplot(data=df, x='SNAME', y='Zbp', ax=ax)
    ax.set_xlabel("")

def iq_values(ax, df):
    ax.set_title("Interquartile coverage values", loc="left")

    ax.set_ylabel("Coverage")
    boxes = [{'label' : r['SNAME'], 'whislo' : 3.0, 'q1' : r['Q25'], 'med' : r['Median'],
    'q3' : r['Q75'], 'whishi' : r['Q75'] * 1.5, 'fliers' : []} for i, r in df.iterrows()]
    #df[["SNAME", "Q25", "Median", "Q75"]].plot(x= "SNAME", ax=ax)

    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels())
    ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))
    ax.bxp(boxes, positions=[x for x in range(len(ax.get_xticks()))], showfliers=False)
    ax.set_xlabel("")

def cov_plot(ax, df):
    ax.set_title("Coverage Mean", loc="left")

    ax.set_ylabel("Mean X Coverage")

    avg = df['Mean'].mean()
    maxv = df['Mean'].max()

    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
    ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))
    sns.barplot(data=df, x='SNAME', y='Mean', ax=ax)
    ax.set_xlabel("")
    ax.axhline(avg, linestyle='--')
    ax.annotate(f'Avg X Coverage: {avg:.2f}', (int(len(df.index) / 2), avg), xytext=(0.35, 0.55), textcoords='axes fraction', arrowprops=dict(facecolor='blue', arrowstyle='wedge'), fontsize=14, color='b')

def heatmap_plot(ax, df, cbar_ax):
    ax.set_xlabel("Samples")
    tdf = df[["SNAME", 'Zbp', 'sub15', 'sub30', 'gt30']].set_index('SNAME').copy()
    tdf.index = [''.join(col) for col in tdf.index.values]
    tdf = tdf.T
    sns.heatmap(tdf, ax=ax, cbar_ax=cbar_ax, yticklabels=['Zbp', '<15X', '15-30X', '>30X'], cbar_kws={'format' : FuncFormatter(number_formatter)})

cols = ["SNAME", "ThreshBp", "Zbp", "Mean", "Q25", "Median", "Q75", "Max", "Stdev",
'sub15', 'sub30', 'gt30']

data = defaultdict(list)
linecount = 0
with open(sys.argv[1], 'r') as input:
    head = input.readline().rstrip().split()
    for l in input:
        s = l.rstrip().split()
        linecount += 1
        for i, k in enumerate(s):
            data[head[i]].append(k)

df = pd.DataFrame(data)
df["Zbp"] = df["Zbp"].astype('int32')
for i in ("Q25", "Mean", "Median", "Q75", "Stdev", "Max"):
    df[i] = df[i].astype('float')

for i in ('sub15', 'sub30', 'gt30'):
    df[i] = df[i].astype('float')

print(df)
splits = int(sys.argv[2])

figures = list()
# New feature: split up the list of samples if they exceed a theshold and print their plots on separate PDF pages
for i in range(0, linecount, splits):
    end = i + splits
    if end >= linecount:
        end = linecount -1
    tdf = df.iloc[i:end, :].copy()
    (fig, axis) = plt.subplots(ncols=2, nrows=4, sharex='col', constrained_layout=True, gridspec_kw={'width_ratios' : [100, 1]})

    axis[0,1].remove()
    axis[1,1].remove()
    axis[2,1].remove()

    ax = axis[0,0]
    cov_thresh(ax, tdf)

    ax = axis[1,0]
    iq_values(ax, tdf)

    ax = axis[2,0]
    cov_plot(ax, tdf)

    ax = axis[3,0]
    heatmap_plot(ax, tdf, axis[3,1])

    fig.set_size_inches(21,12)
    fig.tight_layout()
    figures.append(fig)

with PdfPages(sys.argv[3]) as pdf:
    for f in figures:
        pdf.savefig(sys.argv[3], bbox_inches='tight')

print('Fini')
