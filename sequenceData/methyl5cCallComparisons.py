# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:10:30 2022

@author: derek.bickhart-adm
"""

import argparse
import networkx as nx
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict
import contextlib
import gzip
import math
from scipy.stats.kde import gaussian_kde
matplotlib.use('Agg')

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A wrapper script for the FlashFry CRISPR design tool"
            )
    parser.add_argument('-b', '--bisulfite', 
                        help="Input bisulfite sequence data file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename for plots and tables",
                        required=True, type=str,
                        )
    parser.add_argument('-p', '--pacbio', 
                        help="Input PacBio bed file.",
                        required=True, type=str
                        )
    parser.add_argument('-n', '--nanopore',
                        help="Input Nanopore bed file",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    
    worker = sparseArrayHelper()
    
    # Load Bisulfite data
    worker.loadBisulf(args.bisulfite, limit="haplotype1-0000001")
    print("Loaded Bisulfite")
    
    # Load PacBio data
    worker.loadPacB(args.pacbio, limit="haplotype1-0000001")
    print("Loaded PacB")
    
    # Load ONT data
    worker.loadONT(args.nanopore, limit="haplotype1-0000001")
    print("Loaded ONT")
    
    # Generate plots
    (fig, axis) = plt.subplots(ncols=2, nrows=4)
    axs = axis.flatten().tolist()
    axs.reverse()

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'r'):
    if filename.endswith('gz'):
        fh = gzip.open(filename, 'rt', encoding='utf-8')
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        fh.close()

def plotDistribution(workerdf, plotName):
    sns.violinplot(data=workerdf)
    plt.savefig(plotName)
 
def plotCorrelation(keys, workerdf, plotName):
    corrMatrix = workerdf.corr()
    sns.heatmap(corrMatrix, annot=True)
    plt.savefig(plotName)

def plotHeatMap(term1, term2, worker, plotName):
    x = list(worker.getCol(term1))
    y = list(worker.getCol(term2))
    
    Z, xedges, yedges = np.histogram2d(x, y)


    plt.pcolormesh(xedges, yedges, Z.T)
    plt.colorbar()
    plt.savefig(plotName)

def func(s):
    tokens = s.split('_')
    for i in range(2):
        token = s[i]
        try:
            v = int(token)
            return ('.'.join(tokens[0:i]), v)
        except ValueError:
            pass
    return (s, 0)
    
class sparseArrayHelper:
    
    def __init__(self):
        self.data = defaultdict(dict)
        self.ucscToIdx = dict()
        self.sampleToIdx = dict()
        self.sparse = None
        self.srtArray = []
        
    def getCol(self, second):
        if len(self.srtArray) == 0:
            self.srtArray = sorted(list(self.data.keys()), key=func)
        vec = []
        for i, c in enumerate(self.srtArray):
            self.ucscToIdx[c] = i
            if math.isnan(self.data[c].get(second, 0.0)):
                print(f'Error NaN! {c} {second}')
            vec.append(self.data[c].get(second, 0.0))
                
        return vec
        
    def toSparse(self):
        if len(self.srtArray) == 0:
            self.srtArray = sorted(list(self.data.keys()), key=func)
        numSamples = len(self.sampleToIdx.keys())
        
        
        for i, c in enumerate(self.srtArray):
            self.ucscToIdx[c] = i
            for j, k in enumerate(list(self.sampleToIdx.keys())):
                self.sparse[i, j] = self.data[c].get(k, 0.0)
                
        self.data = dict()
        
    def loadPacB(self, file, limit = None, ftype="PacB"):
        num = len(self.sampleToIdx.keys())
        self.sampleToIdx[ftype] = num
        with smartFile(file, 'r') as input:
            for l in input: 
                s = l.rstrip().split()
                if limit != None:
                    if s[0] != limit:
                        continue
                self.data[f'{s[0]}_{s[2]}'][ftype] = float(s[8])
                
    def loadBisulf(self, file, limit = None, ftype="BiSulf"):
        num = len(self.sampleToIdx.keys())
        self.sampleToIdx[ftype] = num
        
        with smartFile(file, mode='rt') as input:
            for l in input:
                s = l.rstrip().split()
                s[0] = s[0].replace('_', '-')
                if limit != None:
                    if s[0] != limit:
                        continue             
                
                self.data[f'{s[0]}_{s[2]}'][ftype] = float(s[5]) * 100.0
                
    def loadONT(self, file, limit=None, ftype="ONT"):
        num = len(self.sampleToIdx.keys())
        self.sampleToIdx[ftype] = num
        
        with smartFile(file, mode='rt') as input:
            for l in input:
                s = l.rstrip().split()
                if limit != None:
                    if s[0] != limit:
                        continue
                if s[10] == "nan":
                    continue
                self.data[f'{s[0]}_{s[2]}'][ftype] = float(s[10]) 
                
    def toDataFrame(self, cols):
        return pd.DataFrame({k : self.getCol(k) for k in cols}, columns=cols)

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
