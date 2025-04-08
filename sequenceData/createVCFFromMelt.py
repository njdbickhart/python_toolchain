# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 12:25:57 2025

@author: Derek.Bickhart
"""

import sys
import re
from collections import defaultdict

usage = f'python {sys.argv[0]} <Input melted file> <template vcf> <output file>'
if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
coords = dict()
with open(sys.argv[2], 'r') as input:
    for l in input:
        if l.startswith('#CHROM'):
            continue
        s = l.rstrip().split()
        # ref allele, alt allele, qual, info
        coords[f'{s[0]}:{s[1]}'] = (s[3], s[4], s[5], s[7])
        
vcfdata = defaultdict(lambda : defaultdict(list))
animals = set()
with open(sys.argv[1], 'r') as input:
    for l in input:
        s = l.rstrip().split()
        vcfdata[s[0]][s[2]] = s[1]
        animals.add(s[2])

vcfheader = sorted(list(animals))

with open(sys.argv[3], 'w') as output:
    vcfheadstr = "\t".join(vcfheader)
    output.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{vcfheadstr}\n')
    for p, v in vcfdata.items():
        cvals = p.split(":")
        cdata = coords[p]
        outstr = f'{cvals[0]}\t{cvals[1]}\t.\t{cdata[0]}\t{cdata[1]}\t{cdata[2]}\t.\t{cdata[3]}\tGT'
        for a in vcfheader:
            gt = v[a].split(':')
            acgts = list()
            for g in gt:
                if g.startswith('?'):
                    acgts.extend(['.', '.'])
                    break
                if g == cdata[0]:
                    acgts.append('0')
                else:
                    acgts.append('1')
            acgts = sorted(acgts)
            n = "/".join(acgts)
            outstr += f'\t{n}'
        outstr += '\n'
        output.write(outstr)
        
print('Fini!')