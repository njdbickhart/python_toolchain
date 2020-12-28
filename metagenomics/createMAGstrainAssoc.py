# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 08:04:08 2020

@author: derek.bickhart-adm
"""

import argparse
import sys
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to combine mash bin associations with strain level estimations"
            )
    parser.add_argument('-r', '--reference', 
                        help="Reference strain count file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are (basename).table and (basename).melt",
                        required=True, type=str,
                        )
    parser.add_argument('-n', '--name',
                        help="Assembly strain name file; May be entered more than once",
                        action="append", default=[]
                        )
    parser.add_argument('-a', '--association',
                        help="Mash distance association file; May be entered more than once",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--strain',
                        help="Association strain count file; May be entered more than once",
                        action="append", default=[]
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    if len(args.name) != len(args.association) or len(args.name) != len(args.strain):
        print("Argument input must match for the strain, associations and names!")
        parser.print_help()
        sys.exit()
    elif len(args.name) == 0 or len(args.association) == 0 or len(args.strain) == 0:
        parser.print_help()
        sys.exit()
        
    worker = StrainAssoc(args.reference, args.output)
    
    for n, a, s in zip(args.name, args.association, args.strain):
        print(f'{n}\t{a}\t{s}')
        worker.addassoc(n, a)
        
        worker.addstrain(n, s)
        
    worker.createOut()
        

class StrainAssoc:

    def __init__(self, refstrain, output):
        self.output = output
        
        self.straincount = dict()
        with open(refstrain, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                s[0] = s[0].replace('.', '_')
                self.straincount[s[0]] = int(s[1])
        
        # name -> refbin -> assocbin
        self.assoc = defaultdict(dict)
        # name -> assocbin -> distance
        self.distance = defaultdict(dict)
        # name -> assocbin -> assocstraincount
        self.astrain = defaultdict(dict)
        
    def addassoc(self, name, assocfile):
        with open(assocfile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                self.assoc[name][s[0]] = s[1]
                self.distance[name][s[1]] = float(s[2])
                
    def addstrain(self, name, strainfile):
        with open(strainfile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                s[0] = s[0].replace('.', '_')
                self.astrain[name][s[0]] = int(s[1])
                
    def createOut(self):
        names = set()
        refbins = set()
        with open(self.output + '.melt', 'w') as out:
            # a negative strain count means more in the reference
            out.write("Refbin\tAssocBin\tBinDist\tStrainDelta\n")
            for n in self.assoc:
                names.add(n)
                for r in self.straincount:
                    refbins.add(r)
                    if r in self.assoc[n]:
                        abin, dist, strain = self.getDetails(n, r)
                        sd = strain - self.straincount[r]
                        out.write('\t'.join([r, abin, str(dist), '{:.2f}'.format(sd)]) + "\n")
                        
        with open(self.output + '.table', 'w') as out:
            out.write("Refbin\tRefStrain\t" + "\t".join([f'{x}_bin\t{x}_dist\t{x}_strain' for x in names]) + "\n")
            for r in refbins:
                out.write(f'{r}\t{self.straincount[r]}\t')
                for n in names:
                    abin, dist, strain = self.getDetails(n, r)
                    out.write(f'{abin}\t{dist}\t{strain}\t')
                out.write("\n")
                        
    def getDetails(self, name, refbin):
        abin = self.assoc[name][refbin]
        return abin, self.distance[name][abin], self.astrain[name][abin]
                    

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
