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
        self.complete = dict()
        self.contam =dict()
        with open(refstrain, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                s[0] = s[0].replace('.', '_')
                self.straincount[s[0]] = int(s[1])
                self.complete[s[0]] = float(s[2])
                self.contam[s[0]] = float(s[3])
        
        # name -> refbin -> [assocbin]
        self.assoc = defaultdict(dict)
        # name -> assocbin -> distance
        self.distance = defaultdict(dict)
        # name -> assocbin -> assocstraincount
        self.astrain = defaultdict(dict)
        # name -> assocbin -> completeness
        self.acomp = defaultdict(dict)
        # name -> assocbin -> contamination
        self.acont = defaultdict(dict)
        
    def addassoc(self, name, assocfile):
        with open(assocfile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                if not s[0] in self.assoc[name]:
                    self.assoc[name][s[0]] = []
                self.assoc[name][s[0]].append(s[1])
                self.distance[name][s[1]] = float(s[2])
                
    def addstrain(self, name, strainfile):
        with open(strainfile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                s[0] = s[0].replace('.', '_')
                self.acomp[name][s[0]] = float(s[2])
                self.acont[name][s[0]] = float(s[3])
                self.astrain[name][s[0]] = int(s[1])
                
    def createOut(self):
        names = set()
        refbins = set()
        with open(self.output + '.melt', 'w') as out:
            # a negative strain count means more in the reference
            out.write("Refbin\tRefComp\tRefCont\tDataset\tAssocBin\tAssocComp\tAssocCont\tBinDist\tStrainDelta\n")
            for n in self.assoc:
                names.add(n)
                for r in self.straincount:
                    refbins.add(r)
                    if r in self.assoc[n]:
                        for abin in self.assoc[n][r]:
                            dist, comp, cont, strain = self.getDetails(n, abin)
                            if dist == -1:
                                print(f'Missing bin: {abin} from {n} in ref {r}')
                            sd = strain - self.straincount[r]
                            out.write('\t'.join([r, str(self.complete[r]), str(self.contam[r]), n, abin, str(comp), str(cont), str(dist), '{:.2f}'.format(sd)]) + "\n")
                        
        with open(self.output + '.table', 'w') as out:
            out.write("Refbin\tRefComp\tRefCont\tRefStrain\t" + "\t".join([f'{x}_bin\t{x}_comp\t{x}_cont\t{x}_dist\t{x}_strain' for x in names]) + "\n")
            for r in refbins:
                out.write(f'{r}\t{self.complete[r]}\t{self.contam[r]}\t{self.straincount[r]}\t')
                for n in names:
                    bstr = []; cstr = []; ostr = []; dstr = []; sstr = [];
                    if r in self.assoc[n]:
                        for abin in self.assoc[n][r]:
                            dist, comp, cont, strain = self.getDetails(n, abin)
                            bstr.append(abin)
                            cstr.append(f'{comp:.2f}')
                            ostr.append(f'{cont:.2f}')
                            dstr.append(f'{dist:.2f}')
                            sstr.append(strain)
                    else:
                        bstr.append('NA')
                        cstr.append('NA')
                        ostr.append('NA')
                        dstr.append('NA')
                        sstr.append('NA')
                    out.write("{}\t{}\t{}\t{}\t{}\t".format(";".join(bstr), ";".join(cstr), ";".join(ostr), ";".join(dstr), ";".join([str(x) for x in sstr])))
                out.write("\n")
                        
    def getDetails(self, name, abin):
        if not abin in self.distance[name] or not abin in self.astrain[name]:
            return -1, -1, -1, -1
        return self.distance[name][abin], self.acomp[name][abin], self.acont[name][abin], self.astrain[name][abin]
                    

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
