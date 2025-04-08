import os
import sys

usage = f'python {sys.argv[0]} <input mash distances> <contaminant name file> <output binary detection>'

if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)

names = set()

with open(sys.argv[2], 'r') as input:
    for l in input:
        l = l.rstrip()
        names.add(l)

