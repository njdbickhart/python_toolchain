# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:21:48 2022

@author: derek.bickhart-adm
"""

import argparse

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script designed to count how many chromosomes are telomere capped"
            )
    parser.add_argument('-f', '--fai', 
                        help="Input reference fasta index file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name.",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bed',
                        help="Input repeatmasker bed file",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    chrs = dict() 
    with open(args.fai, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            chrs[s[0]] = chromosome(s[0], int(s[1]))
    
    with open(args.bed, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            chrs[s[0]].processRMLine(s)
            
    with open(args.output, 'w') as output:
        for k, v in chrs.items():
            output.write(v.getOutLine())


    
class chromosome:
    
    def __init__(self, name, length):
        self.name = name
        self.len = length
        self.beg = False
        self.end = False
        
    def checkTelomere(self, rep):
        rep = rep.replace('(', '').replace(')', '').replace('n', '').rstrip()
        isTel = False
        for i in range(len(rep)):
            if rep == 'AACCCT':
                isTel = True
                break 
            trep = rep.translate('ACGT', 'TGCA').reverse()
            if trep == 'AACCCT':
                isTel = True 
                break 
            rep = rep[-1] + rep[:-2]
        return isTel 
    
    def processRMLine(self, segs):
        beginning = int(segs[1])
        if beginning <= 10000:
            if self.checkTelomere(segs[4]):
                self.beg = True
        elif beginning >= self.len - 10000:
            if self.checkTelomere(segs[4]):
                self.end = True
    
    def getOutLine(self):
        return f'{self.name}\t{self.len}\t{self.beg}\t{self.end}\n'
    
if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
