import os
import sys

usage = f'python {sys.argv[0]} <input mash distances> <contaminant name file> <filtered output> <output binary detection>'

if len(sys.argv) != 5:
    print(usage)
    sys.exit(-1)

names = set()

with open(sys.argv[2], 'r') as input:
    for l in input:
        l = l.rstrip()
        names.add(l)

with open(sys.argv[1], 'r') as input, open(sys.argv[3], 'w') as filtered, open(sys.argv[4], 'w') as binary:
    #start = False
    for l in input:
        #if l.startswith('Writing output'):
        #    start = True
        #    continue

        #if start:
        s = l.rstrip().split()
        if s[4] in names:
            filtered.write(l)

            kvalues = s[1].split('/')
            kratio = int(kvalues[0]) / int(kvalues[1])

            # Create a binary state if the identification is very strong
            if kratio > 0.8 or int(s[2]) > 20:
                binary.write(f'{s[4]}\n')