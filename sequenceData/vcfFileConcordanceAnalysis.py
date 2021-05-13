# -*- coding: utf-8 -*-
"""
Created on Thu May 13 08:56:41 2021

@author: derek.bickhart-adm
"""

import argparse
import gzip
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Carrier concordance analysis from gzipped vcf files"
            )
    parser.add_argument('-v', '--vcf', 
                        help="Input gzipped vcf file. Required",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are (basename).tab and (basename).pdf",
                        required=True, type=str,
                        )
    parser.add_argument('-c', '--chromosome',
                        help="Chromosome for filtering coordinates",
                        required=True, type=str,
                        )
    parser.add_argument('-s', '--start',
                        help="Start base for filtering coordinates",
                        required=True, type=int,
                        )
    parser.add_argument('-e', '--end',
                        help="End base for filtering coordinates",
                        required=True, type=int,
                        )
    parser.add_argument('-l', '--list',
                        help="Newline delimited list of carrier animals. Carriers are assumed to be heterozygotes and all others homozygous reference",
                        required=True, type=str,
                        )
    
    parser.add_argument('-m', '--min',
                        help="Minimum y value to plot (float)",
                        default=0.90, type=float,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Read in carrier list
    clist = list()
    with open(args.list, 'r') as carriers:
        for l in carriers:
            clist.append(l.rstrip())
            
    # Create workhorse
    worker = variantList(clist, args.chromosome, args.start, args.end, args.output)
    
    # Process vcf file
    worker.processVCF(args.vcf)
    
    # Print out tab file
    worker.printOutTable()
    
    # Print out plot
    worker.plotData(args.min)
    
    print("Data has been printed and plotted!")
    
    
class variantList:
    
    def __init__(self, carriers, chrom, start, end, output):
        self.carriers = carriers
        self.chrom = chrom
        self.start = start
        self.end = end
        self.output = output
        
        self.carIdx = set()
        self.sampTot = 0
        self.varPos = list()
        self.varCon = list()
        self.varPerc = list()
        self.varAnn = list()
        self.hasAnn = True
        
        
    def _getCarrierIdx(self, cline):
        segs = cline.split()
        for x in range(len(segs)):
            self.sampTot += 1
            if segs[x] in self.carriers:
                self.carIdx.add(x)
        self.sampTot -= 9
        
        print(f'Total Samples: {self.sampTot}')
        if len(self.carIdx) == 0:
            print("Found zero carriers in the header! Exiting...")
            sys.exit(-1)
                
    def _estimateConcordance(self, segs):
        self.varPos.append(int(segs[1]))
        varCon = 0
        
        if self.hasAnn:
            dsegs = segs[7].split(';')
            astr = ''
            for d in dsegs:
                if d.startswith('ANN='):
                    isegs = d.split(',')[0].split('|')
                    astr = isegs[2]
            if astr == '':
                self.hasAnn = False
            else:
                self.varAnn.append(astr)
        
        for i in range(9, len(segs)):
            gtypes = segs[i][0:3:2]
            naltcount = 0
            for g in gtypes:
                if g != '0' and g != '.':
                    naltcount += 1
            if i in self.carIdx and naltcount == 1:
                varCon += 1
            elif not i in self.carIdx and naltcount == 0:
                varCon += 1
        self.varCon.append(varCon)        
        self.varPerc.append(float(varCon / self.sampTot))
            
                
    def processVCF(self, vcf):
        with gzip.open(vcf, mode='rt') as data:
            header = False
            for l in data:
                if l.startswith('##'):
                    continue
                elif l.startswith('#CHROM'):
                    header = True
                    self._getCarrierIdx(l)
                    continue
                elif l.startswith('#'):
                    continue
                    
                if not header:
                    print("Error! Did not encounter a parseable header line in this vcf! Exiting...")
                    sys.exit(-1)
                    
                segs = l.rstrip().split()
                pos = int(segs[1])
                if segs[0] == self.chrom and pos <= self.end and pos >= self.start:
                    self._estimateConcordance(segs)
        print(f'Gathered concordance stats for {len(self.varPos)} variant positions')
        
    def printOutTable(self):
        with open(self.output + ".tab", 'w') as out:
            out.write("Chr\tPos\tConcord\tPercent")
            if self.hasAnn:
                out.write("\tAnn\n")
            else:
                out.write("\n")
            for i in range(len(self.varPos)):
                out.write(f'{self.chrom}\t{self.varPos[i]}\t{self.varCon[i]}\t{self.varPerc[i]}')
                if self.hasAnn:
                    out.write(f'\t{self.varAnn[i]}\n')
                else:
                    out.write("\n")
                
    def plotData(self, minimum):
        df = pd.DataFrame({
            "Pos" : self.varPos,
            "Perc" : self.varPerc
            })
        df = df.sort_values(by=['Pos'])
        #df['Quantile'] = pd.qcut(df['Perc'], [0.90, 0.95, 1.0], labels=False)
        
        if self.hasAnn:
            df['Ann'] = self.varAnn
        
        ymin = minimum
        
        quantiles = df['Perc'].apply(lambda x : x if x > 0.95 else 0.95).rename("Magnitude")
        
        g = sns.scatterplot(x=df['Pos'], y=df['Perc'], hue=quantiles)
        g.set(ylim=(ymin, 1.01))
        
        if self.hasAnn:
            for row in range(0, df.shape[0]):
                if df.Ann[row] == "HIGH":
                    plt.text(df.Pos[row], df.Perc[row] + 0.01, f'{df.Ann[row]};{df.Pos[row]}')
                    plt.plot([df.Pos[row], df.Pos[row]], [df.Perc[row], df.Perc[row] + 0.01], 'k-')
                
        
        plt.savefig(self.output + '.pdf')       
                    
                

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
