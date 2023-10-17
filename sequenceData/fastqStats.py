#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip
import collections
import seaborn as sns
import pandas as pd
import numpy as np
import glob
import argparse

version = "0.2"

plt.style.use("seaborn-colorblind")

matplotlib.rcParams['axes.titlesize']=13
matplotlib.rcParams['axes.labelsize']=12
matplotlib.rcParams['xtick.labelsize']=11
matplotlib.rcParams['ytick.labelsize']=11
matplotlib.rcParams['figure.titlesize']=20

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate summary information and PDF plots of Fastq Stats. Version:" + version
            )
    parser.add_argument('-f', '--file',
                        help="A single fastq file (can be gzipped). If specified at the same time as 'b', this is preferentially run.",
                        type=str, default=None
                        )
    parser.add_argument('-b', '--base',
                        help="A folder containing fastq files. Not processed if 'f' is also specified",
                        type=str, default=None
                        )
    parser.add_argument('-s', '--separate',
                        help="A flag to print plots for each fastq file separately [default: False]",
                        type=bool, default=False
                        )
    parser.add_argument('-o', '--output',
                        help="Output file Basename",
                        type=str, required=True
                        )

    return parser.parse_args()

def main(args):
    # Get a list of fastq files in the directory
    files = []
    if args.file != None or args.base != None:
        if args.file != None:
            # Single file checking
            files.append(args.file)
        else:
            # Folder checking
            for exts in ('*.fq', '*.fastq', '*.fq.gz', '*.fastq.gz'):
                files.extend(glob.glob('{}/{}'.format(args.base, exts)))

    data = fastq_info()

    for filename in files:
        data = get_fastq_info(filename, data)
        if args.separate:
            prefix = filename.split('.')
            df = data.get_pddf()
            print_plot(df, args.output, prefix[0])
            data = fastq_info()

    if not args.separate:
        df = data.get_pddf()
        print_plot(df, args.output)
        print(df.describe())

def smartFile(filename : str, mode : str = 'r'):
    fh = None
    if filename.endswith('.gz') and mode == 'r':
        fh = gzip.open(filename, mode='rt')
    elif filename.endswith('.gz') and mode == 'w':
        fh = gzip.open(filename, mode='wt')
    else:
        fh = open(filename, mode)
    return fh

def print_plot(df, outbase, outsuffix = "total"):
    df.to_csv(outbase + "_" + outsuffix + "_raw.tab")
    df.describe().to_csv(outbase + "_" + outsuffix + "_stats.tab")
    plt = plot_fastq_info(df, outbase)
    plt.savefig(outbase + "_" + outsuffix + "_plots.png")

def get_fastq_info(filename, fqinfo):
    infile = smartFile(filename, 'r')

    for head, seq, qual in fastq_reader_fh(infile):
        fqinfo.seq_lengths.append(len(seq))

        fqinfo.mean_qualities.append(round(np.mean(bytearray(qual, "ascii")) - 33, 2))

        fqinfo.kmers_start.append(seq[0:4])
        fqinfo.kmers_end.append(seq[-4:])

        ntc = collections.Counter()
        ntc.update(seq)
        fqinfo.nts_A.append(ntc['A'])
        fqinfo.nts_G.append(ntc['G'])
        fqinfo.nts_T.append(ntc['T'])
        fqinfo.nts_C.append(ntc['C'])
        fqinfo.nts_U.append(ntc['U'])

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

def plot_nt_content(ax, df):
    ax.set_title("Nucleotide Content", loc="left")
    ax.set_xlabel("NT")
    ax.set_ylabel("Percent")
    nt = {}
    nt['A'] = np.sum(df.nt_A)
    nt['T'] = np.sum(df.nt_T)
    nt['G'] = np.sum(df.nt_G)
    nt['C'] = np.sum(df.nt_C)
    nt['U'] = np.sum(df.nt_U)
    s = np.sum((nt['A'], nt['T'], nt['G'], nt['C'], nt['U']))
    for k,v in nt.items(): nt[k] = (v / s) * 100
    if nt['U'] > 0:
        sns.barplot(x=['G', 'C', 'A', 'U'], y=[nt['G'], nt['C'], nt['A'], nt['U']], ax=ax)
    else:
        sns.barplot(x=['G', 'C', 'A', 'T'], y=[nt['G'], nt['C'], nt['A'], nt['T']], ax=ax)

def plot_length_binned(ax, df, num_bins = 1000):
    ax.set_title("Histogram of sequence lengths 1000 bins", loc="left")
    ax.set_xlabel("Lengths")
    ax.set_ylabel("Count (Log)")
    ax.hist(df.seq_length, num_bins)
    ax.set_yscale('log')
    xlabel_pos_right(ax)

def plot_length_percentile(ax, df, percentile=80):
    cutoff = int(round(np.percentile(df.seq_length, percentile)))
    num_bins = int(round(cutoff))
    ax.set_title("Histogram of Sequence Lengths from 0 to " + str(cutoff) + " (" + str(percentile) + "%ile) in " + str(num_bins) + " bins", loc="left")
    ax.set_xlabel("Length")
    ax.set_ylabel("Count (Log)")
    ax.set_yscale('log')
    xlabel_pos_right(ax)
    ax.hist(df.seq_length, num_bins, range=(0,cutoff))

def plot_gc_bins(ax, df, num_bins=20):
    gcpercs = ((df['nt_G'] + df['nt_C']) / (df['seq_length'])) * 100
    ax.set_title("Histogram of GC percentage bins", loc="left")
    ax.set_xlabel("GC Bins")
    ax.set_ylabel("Count")
    ax.hist(gcpercs, num_bins)

def plot_kmer_start(ax, df):
    df.kmers_start.value_counts().head(n=40).plot.bar(ax=ax)
    ax.set_title("Start of Sequence k-mer Content (k=4) top 40", loc="left")
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    xlabel_pos_right(ax)

def plot_kmer_end(ax, df):
    df.kmers_end.value_counts().head(n=40).plot.bar(ax=ax)
    ax.set_title("End of Sequence k-mer Content (k=4) top 40", loc="left")
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    xlabel_pos_right(ax)

def plot_qualities(ax,df):
    ax.set_title("Histogram of Mean Qualities", loc="left")
    ax.hist(df.mean_quality, 160, range=(0,40)) #, s=point_size) #, range=(0,20000))
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Quality")
    #xlabel_pos_right(ax)

def plot_nucleotides_per_length(ax, df):
    blens = {}
    for l in df.seq_length:
        bin_l = int(l / 100) * 100
        if bin_l in blens:
            blens[bin_l] += l
        else:
            blens[bin_l] = l

    len_sum = df["seq_length"].sum()
    mean_sum = len_sum / 2.0
    n50 = len_sum
    l = 0
    for k in sorted(df["seq_length"], reverse=True):
        l += k
        if l > mean_sum:
            break
        n50 = k

    print("N50\t" + str(n50))

    ax.set_title("Nucleotides per Sequence Length in 100nt bins (N50:" + str(n50) + ")", loc="left")
    ax.set_xlabel("Length")
    ax.set_ylabel("GBases")
    lens = np.asfarray(list(blens.values()))
    ax.bar(list(blens.keys()), lens / float(1000 ** 3), width=100) #, s=point_size) #, range=(0,20000))
    ax.axvline(x=n50, color="red")
    xlabel_pos_right(ax)

def plot_fastq_info(df, outbase):
    (fig, axis) = plt.subplots(ncols=2, nrows=4)
    axs = axis.flatten().tolist()
    axs.reverse()

    # right side
    plot_length_binned(axs.pop(), df, 1000)

    plot_length_percentile(axs.pop(), df, 80)

    plot_qualities(axs.pop(), df)

    plot_nucleotides_per_length(axs.pop(), df)

    # left side
    plot_kmer_start(axs.pop(), df)

    plot_kmer_end(axs.pop(), df)

    plot_nt_content(axs.pop(), df)

    plot_gc_bins(axs.pop(), df, 20)


    plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)
    plt.suptitle(outbase + " (" + str(len(df)) + " reads, " + str(round(df.seq_length.sum() / (1000.0**3), 2)) + " Gbp)")
    fig.set_size_inches(21,12)
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    return plt


class fastq_info:

    def __init__(self):
        self.seq_lengths = []
        self.mean_qualities = []
        self.kmers_start = []
        self.kmers_end = []
        self.nts_A = []
        self.nts_G = []
        self.nts_T = []
        self.nts_C = []
        self.nts_U = []

    def get_pddf(self):
        df = pd.DataFrame({ "seq_length": self.seq_lengths,
                            "mean_quality": self.mean_qualities,
                            "kmers_start": self.kmers_start,
                            "kmers_end": self.kmers_end,
                            "nt_A": self.nts_A,
                            "nt_G": self.nts_G,
                            "nt_T": self.nts_T,
                            "nt_C": self.nts_C,
                            "nt_U": self.nts_U})
        return df.sort_values(by=['seq_length'])


if __name__ == "__main__":
    args = parse_user_input()
    main(args)
