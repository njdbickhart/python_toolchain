#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip
from collections import defaultdict, Counter
import seaborn as sns
import pandas as pd
import numpy as np
import glob
import argparse

version = "0.1"

plt.style.use("seaborn-colorblind")

matplotlib.rcParams['axes.titlesize']=13
matplotlib.rcParams['axes.labelsize']=12
matplotlib.rcParams['xtick.labelsize']=11
matplotlib.rcParams['ytick.labelsize']=11
matplotlib.rcParams['figure.titlesize']=20

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Compare two ONT fastqs to plot differences. Version:" + version
            )
    parser.add_argument('-f', '--first', 
                        help="A single fastq file (can be gzipped). The first fastq for comparison.",
                        type=str, required=True
                        )
    parser.add_argument('-s', '--second',
                        help="A single fastq file (can be gzipped). The second fastq for comparison.",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Output file Basename",
                        type=str, required=True
                        )
    
    return parser.parse_args()

def main(args):
    # Get a list of fastq files in the directory
    files = [args.first, args.second]

    data = fastq_comp()
    
    for fnum in range(len(files)):	
        data = get_fastq_info(files[fnum], data, fnum + 1)
    
    df = data.get_pddf()
    print_plot(df, args.output)
    print(df.describe())

def smartFile(filename : str, mode : str = 'r'):
    fh = None
    if filename.endswith('.gz'):
        fh = gzip.open(filename, mode='rt')
    else:
        fh = open(filename, mode)
    return fh
    
def print_plot(df, outbase, outsuffix = "total"):
    with open(outbase + "_" + outsuffix + "_stats.tab", 'w') as out:
        out.write(df.describe() + "\n")
    #df.to_csv(outbase + "_" + outsuffix + "_stats.tab")
    plt = plot_fastq_info(df, outbase)
    plt.savefig(outbase + "_" + outsuffix + "_plots.png")
    
def getRString(fqname):
    fsegs = fqname.split()
    fdict = {}
    for f in fsegs[1:]:
        t = f.split('=')
        fdict[t[0]] = t[1]
    
    if not(set(["read", "ch"]) - fdict.keys()):
        return "{}_{}".format(fdict["read"], fdict["ch"])
    else:
        return "NULL"
            
def get_fastq_info(filename, fqinfo, number):
    infile = smartFile(filename, 'r')

    for head, seq, qual in fastq_reader_fh(infile):
        rstring = getRString(head)
        if rstring == "NULL":
            # The fastq file is either improperly formatted, or we ran into a read error
            continue
        
        flen = len(seq)
        mqual = round(np.mean(bytearray(qual, "ascii")) - 33, 2)
        
        ntc = Counter()
        ntc.update(seq)
        
        ats = ntc['A'] + ntc['T'] + ntc['U']
        gcs = ntc['G'] + ntc['C']
        
        fqinfo.add_data(rstring, number, flen, mqual, ats, gcs)
        
    infile.close()

    return fqinfo

def fastq_reader_fh(infile):
  name = infile.readline().rstrip()
  while True:
    seq = ""
    for s in infile:
      if s[0] == '+':
        break
      else:
        seq += s.rstrip()
    qual = ""
    for q in infile:
      if len(qual) > 0 and  q[0] == '@':
        yield name, seq, qual
        name = q.rstrip()
        break
      else:
        qual += q.rstrip()
    else:
      yield name, seq, qual
      return
  
def xlabel_pos_right(ax):
    label = ax.xaxis.get_label()
    x_lab_pos, y_lab_pos = label.get_position()
    label.set_position([1.0, y_lab_pos])
    label.set_horizontalalignment('right')
    ax.xaxis.set_label(label)

def compare_quals(ax, df, cols):
    ax.set_title("Difference in Histogram quality scores", loc="left")
    data = df.loc[df['Category'] == "meanQuality"].loc[:,cols]
    ax.hist(data, 160, range=(-40, 40))
    ax.set_xlim(left=-40, right=40)
    ax.set_xlabel("Quality difference")
    
def compare_gcs(ax, df, cols):
    ax.set_title("Comparative heatmap for GC percentages", loc = "left")
    columns = [str(x) for x in range(11)]
    data = df.loc[df['Category'] == "GCs"].loc[:, cols]
    data = data.apply(lambda x : x * 10)
    data = data.astype(int)
    mat = np.zeroes(shape=(11,11))
    for x, y in zip(data.loc[:, 0], data.loc[:,1]):
        mat[x, y] += 1
    ax.imshow(mat)
    ax.set_xticklabels(columns)
    ax.set_yticklabels(columns)
    ax.set_xlabel("File 1 GC")
    ax.set_ylabel("File 2 GC")
    
def compare_lengths(ax, df, cols):
    ax.set_title("Comparative scatterplot for seq Lengths", loc="left")
    data = df.loc[df['Category'] == "seqLength"].loc[:, cols]
    data.assign(stdev = data.std(axis=1))
    data.assign(colors = "green" if data.stdev >= 500 else "blue")
    ax.scatter(data.iloc[:,0] / 1000 ** 3, data.iloc[:,1] / 1000 * 3, c=data.colors)
    ax.set_xlabel("File 1 Gbp")
    ax.set_ylabel("File 2 Gbp")

def plot_qualities(ax,df, col):
    ax.set_title("Histogram of Mean Qualities for {}".format(col), loc="left")
    data = df.loc[df['Category'] == "meanQuality"].loc[:,col]
    ax.hist(data, 160, range=(0,40)) #, s=point_size) #, range=(0,20000))
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Quality")
    #xlabel_pos_right(ax)

def plot_nucleotides_per_length(ax, df, col):
    blens = {}
    data = df.loc[df['Category'] == "seqLength"].loc[:,col]
    for l in data.loc[:,col]:
        bin_l = int(l / 100) * 100
        if bin_l in blens:
            blens[bin_l] += l
        else:
            blens[bin_l] = l
            
    len_sum = data.loc[:,col].sum()
    mean_sum = len_sum / 2.0
    n50 = len_sum
    l = 0
    for k in sorted(data.loc[:,col], reverse=True):
        l += k
        if l > mean_sum:
            break
        n50 = k

    print("N50\t" + str(n50))

    ax.set_title("Sequence Length in 100nt bins for {} (N50:".format(col) + str(n50) + ")", loc="left")
    ax.set_xlabel("Length")
    ax.set_ylabel("GBases")
    lens = np.asfarray(list(blens.values()))
    ax.bar(list(blens.keys()), lens / float(1000 ** 3), width=100) #, s=point_size) #, range=(0,20000))
    ax.axvline(x=n50, color="red")
    xlabel_pos_right(ax)
    
def plot_gc_bins(ax, df, col, num_bins=20):
    gcs = df.loc[df['Category'] == "GCs"].loc[:,col]
    ats = df.loc[df['Category'] == "ATs"].loc[:,col]
    gcpercs = (gcs / (gcs + ats)) * 100
    ax.set_title("Histogram of GC percentage bins for file: {}".format(col), loc="left")
    ax.set_xlabel("GC Bins")
    ax.set_ylabel("Count")
    ax.hist(gcpercs, num_bins)
  
def plot_fastq_info(df, outbase):
    (fig, axis) = plt.subplots(ncols=3, nrows=3)
    axs = axis.flatten().tolist()
    axs.reverse()
    
    cols = ["1", "2"]
    
    # right side
    compare_gcs(axs.pop(), df, cols)
    
    compare_quals(axs.pop(), df, cols)
    
    compare_lengths(axs.pop(), df, cols)
    
    # Middle and left
    for c in reversed(cols):
        plot_gc_bins(axs.pop(), df, c)
        
        plot_qualities(axs.pop(), df, c)
        
        plot_nucleotides_per_length(axs.pop(), df, c)
    
    
    plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)
    #plt.suptitle(outbase + " (" + str(len(df)) + " reads, " + str(round(df.seq_length.sum() / (1000.0**3), 2)) + " Gbp)")
    fig.set_size_inches(21,12)
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    return plt
    
    
class fastq_comp:

    def __init__(self):
        self.seq_lengths = defaultdict(list)
        self.mean_qualities = defaultdict(list)
        self.nts_AT = defaultdict(list)
        self.nts_GC = defaultdict(list)
        self.cols = 1
        self.check = False
        self.cats = {"seqLength" : self.seq_lengths, "meanQuality" : self.mean_qualities,
                     "ATs" : self.nts_AT, "GCs" : self.nts_GC}

    def add_data(self, fstr, col, length, qual, ATs, GCs):
        if self.cols < col:
            self.cols += 1
            self.check = True
        
        if self.check:
            if fstr not in self.seq_lengths:
                # Fill in missing values
                for x in range(self.cols - 1):
                    self.seq_lengths[fstr].append(0)
                    self.mean_qualities[fstr].append(0.0)
                    self.nts_AT[fstr].append(0)
                    self.nts_GC[fstr].append(0)
            elif len(self.seq_lengths[fstr]) != self.cols - 1:
                for x in range(len(self.seq_lengths[fstr]), self.cols - 1):
                    self.seq_lengths[fstr].append(0)
                    self.mean_qualities[fstr].append(0.0)
                    self.nts_AT[fstr].append(0)
                    self.nts_GC[fstr].append(0)
        self.seq_lengths[fstr].append(length)
        self.mean_qualities[fstr].append(qual)
        self.nts_AT[fstr].append(ATs)
        self.nts_GC[fstr].append(GCs)
            

    def get_pddf(self):
        df = pd.concat([self._melt(k, v) for k, v in self.cats.items()])
        return df.sort_values(by=['Category'])
    
    def _melt(self, cat, data):
        melted = []
        cols = ["Category"]
        for x in range(self.cols + 1):
            melted.append([])
            if x > 0:
                cols.append("{}".format(x))
            
        for k, v in data.items():
            melted[0].append(cat)
            for x in range(1, self.cols + 1):
                melted[x].append(v[x-1])
        
        tdf = pd.DataFrame(melted, columns = cols)
        return tdf
        
if __name__ == "__main__":
    args = parse_user_input()
    main(args)