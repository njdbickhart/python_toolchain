# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:03:11 2018

@author: dbickhart
"""
import argparse
from typing import Dict, Tuple, List
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A script to filter out only the homozygous haplotypes within a list of samples"
            )
    parser.add_argument('-v', '--vcf', 
                        help="Input vcf file with variant calls",
                        required=True, type=str)
    parser.add_argument('-l', '--list', 
                        help="Input list of animals and haplotype associations",
                        required=True, type=str)
    parser.add_argument('-o', '--output', 
                        help="Base output directory",
                        required=True, type=str)
    parser.add_argument('-s', '--segments',
                        help="Segment coordinates for haplotypes",
                        required=True, type=str)
   
    #parser.add_argument('-m', '--meta',
                        #help="[Optional] Add multiple targets for metagenome assembly",
                        #action='store_true') 
    return parser.parse_args()

def generateAnIndexList(vcfHeader : str, anList : Dict[str, int]) -> {}:
    vcfHeader = vcfHeader.rstrip('\n')
    segs = vcfHeader.split()    
    index = {}
    
    for x in range(9, len(segs)):
        if segs[x] in anList:
            index[segs[x]] = x
            
    return index

class haploregions:
    def __init__(self, limit : int = 2):
        self.limit = limit
        # Lookup of chr segments per haplotype
        self.hapregions = defaultdict
        # Reverse lookup of chr segments by haplotype count
        self.chrregions = defaultdict
        
    def loadHapregion(self, seg : int, chrom : int, start : int, end : int) -> None:
        if seg not in self.hapregions:
            self.hapregions[seg] = defaultdict()
        if chrom not in self.hapregions[seg]:
            self.hapregions[seg][chrom] = defaultdict()
        self.hapregions[seg][chrom][start] = end
        
        if chrom not in self.hapregions:
            self.hapregions[chrom] = defaultdict()
        self.hapregions[chrom][start] = seg

def main(args):
    print("hey")

if __name__ == "__main__":
    args = parse_user_input()
    main(args)