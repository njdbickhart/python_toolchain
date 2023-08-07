# This script is designed to run on a samtools indexed fai file
import numpy as np
import sys
#from collections import defaultdict

if len(sys.argv) != 2:
    print(f'python {sys.argv[0]} <assembly.fai file>')
    sys.exit(-1)


print("File\tCount\tSum\tMean\tMedian\tStdev\tContigN50\n")

sample = sys.argv[1]
sizes = list()

# Now read file and store fasta length lists
with open(sys.argv[1], 'r') as input:
    for l in input:
        l = l.rstrip()
        lsegs = l.split()
        sizes.append(int(lsegs[1]))

sizes.sort(reverse=True)

# Calculate stats and print
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

print(f'{sample}\t{count}\t{sum}\t{mean}\t{median}\t{stdev}\t{n50}\n')
