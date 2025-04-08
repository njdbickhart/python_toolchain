# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 11:48:58 2025

@author: Derek.Bickhart
"""

import sys
import gzip
import re


usage = f'python {sys.argv[0]} <Input gzipped vcf> <output file>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

animallist = dict()  
with gzip.open(sys.argv[1], mode='rt') as data, open(sys.argv[2], 'w') as output:
    for l in data:
        if l.startswith('##'):
            continue
        s = l.rstrip().split()
        if l.startswith('#CHROM'):
            for x in range(9, len(s)):
                animallist[x] = s[x]
            continue
                
        coord = f'{s[0]}:{s[1]}'
        for x in range(9, len(s)):
            gtlist = s[x].split(':')
            gts = re.split(r'[|/]', gtlist[0])
            gtstr = list()
            for g in gts:
                if g == '.':
                    gtstr.append('?')
                    break
                gtstr.append(s[int(g) + 3])
            fgt = ":".join(gtstr)
            output.write(f'{coord}\t{fgt}\t{animallist[x]}\n')