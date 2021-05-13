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
    worker.plotData()
    
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
        
        
    def _getCarrierIdx(self, cline):
        segs = cline.split()
        for x in range(len(segs)):
            self.sampTot += 1
            if segs[x] in self.carriers:
                self.carIdx.add(x)
        self.sampTot -= 9
        
        if len(self.carIdx) == 0:
            print("Found zero carriers in the header! Exiting...")
            sys.exit(-1)
                
    def _estimateConcordance(self, segs):
        self.varPos.append(int(segs[1]))
        varCon = 0
        for i in range(9, len(segs)):
            gtypes = segs[i][0:3:2]
            naltcount = 0
            for g in gtypes:
                if g != '0':
                    naltcount += 1
            if i in self.carIdx and naltcount == 1:
                varCon += 1
            elif not i in self.carIdx and naltcount == 0:
                varCon += 1
        self.varCon.append(varCon)        
        self.varPerc.append(float(self.varCon / self.sampTot))
            
                
    def processVCF(self, vcf):
        with gzip.open(vcf, mode='rt') as data:
            header = False
            for l in data:
                if l.startswith('##'):
                    continue
                elif l.startswith('#CHROM'):
                    header = True
                    self._getCarrierIdx(l)
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
            out.write("Chr\tPos\tConcord\tPercent\n")
            for i in range(len(self.varPos)):
                out.write(f'{self.chrom}\t{self.varPos[i]}\t{self.varCon[i]}\t{self.varPerc[i]}\n')
                
    def plotData(self):
        df = pd.DataFrame({
            "Pos" : self.varPos,
            "Perc" : self.varPerc
            })
        df = df.sort_values(by=['Pos'])
        df['Quantile'] = pd.qcut(df['Perc'], [0.90, 0.95, 1.0], labels=False)
        
        ymin = float(df['Perc'].min())
        
        g = sns.scatterplot(x=df['Pos'], y=df['Perc'], hue=df['Quantile'])
        g.set(ylim=(ymin, None))
        plt.savefig(self.output + '.pdf')       
                    
                

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
