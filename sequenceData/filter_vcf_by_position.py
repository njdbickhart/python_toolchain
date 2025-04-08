# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 11:35:03 2025

@author: Derek.Bickhart
"""

import sys
import gzip


usage = f'python {sys.argv[0]} <Input gzipped vcf> <comma separated coords (chr:pos)> <output file>'
if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
coords = list()
for s in sys.argv[2].split(','):
    (c, d) = s.split(':')
    coords.append((c, d))
    
with gzip.open(sys.argv[1], mode='rt') as data, open(sys.argv[3], 'w') as output:
    for l in data:
        if l.startswith('##'):
            continue
        if l.startswith('#CHROM'):
            output.write(l)
        s = l.rstrip().split()
        for (c, d) in coords:
            if s[0] == c and s[1] == d:
                output.write(l)
            