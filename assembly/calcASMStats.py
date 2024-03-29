# This script is designed to run on a samtools indexed fai file
import numpy as np
import sys
import os
#from collections import defaultdict

if len(sys.argv) < 2:
    print(f'python {sys.argv[0]} <assembly.fai file> <...> <assembly.fai file>')
    sys.exit(-1)


print("File\tCount\tOver1Mbp\tSum\tMean\tMedian\tStdev\tContigN50\tMax\tMin")
files = sys.argv[1:]
for f in files:
    sample = os.path.basename(f)
    sizes = list()

    over1mp = 0

    # Now read file and store fasta length lists
    with open(f, 'r') as input:
        for l in input:
            l = l.rstrip()
            lsegs = l.split()
            sizes.append(int(lsegs[1]))
            if int(lsegs[1]) > 1000000:
                over1mp += 1

    sizes.sort(reverse=True)

    # Calculate stats and print
    maximum = np.max(sizes)
    minimum = np.min(sizes)
    count = len(sizes)
    sum = np.sum(sizes) if count > 0 else 0
    mean = np.mean(sizes) if count > 0 else 0
    median = np.median(sizes) if count > 0 else 0
    stdev = np.std(sizes) if count > 0 else 0

    l50v = int(sum / 2)
    csum = np.cumsum(sizes)
    l50idx = min(csum[csum >= l50v])
    c50idx = np.where(csum == l50idx)

    n50 = sizes[int(c50idx[0])]

    print(f'{sample}\t{count}\t{over1mp}\t{sum}\t{mean:.2f}\t{median:.2f}\t{stdev:.2f}\t{n50}\t{maximum}\t{minimum}')
