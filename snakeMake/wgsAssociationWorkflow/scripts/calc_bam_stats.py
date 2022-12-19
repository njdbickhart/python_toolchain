import os
import sys
import pysam
import time
import numpy as np

usage = f'python3 {sys.argv[0]} <input bam> <min coverage thresh> <sample name> <output stats>'

def get_chromosomes_names(input):

    # opening the bam file with pysam
    bamfile = pysam.AlignmentFile(input, 'rb')
    # query all the names of  the chromosomes in a list
    list_chromosomes = bamfile.references
    list_length = bamfile.lengths
    bamfile.close()
    return list_chromosomes, list_length


def count_depth(chr_name, size, threshold, input):
    """
    Count the depth of the read. For each genomic coordinate return the
    number of reads
    -----
    Parameters :
        chr : (str) name of the chromosome
        threshold : (int) minimum value to count pileup
    -----
    Returns :
        int : count of pileups above threshold
    """
    bp = 0
    nbp = 0
    bases = list()
    bamfile = pysam.AlignmentFile(input, 'rb')
    for pileupcolumn in bamfile.pileup(chr_name):
        depth = pileupcolumn.nsegments
        if depth >= threshold:
            bp += 1
        else:
            nbp += 1
        bases.append(depth)
    bamfile.close()
    count = len(bases)
    mean = np.mean(bases) if count > 0 else 0
    median = np.median(bases) if count > 0 else 0
    stdev = np.std(bases) if count > 0 else 0
    return (bp, nbp, mean, median, stdev)

if len(sys.argv) != 5:
    print(usage)
    sys.exit(-1)

bam = sys.argv[1]
threshold = int(sys.argv[2])
sname = sys.argv[3]
output = sys.argv[4]
print(f'Starting depth estimate for bam: {bam} at threshold {threshold}')

class container:

    def __init__(self):
        self.bp = list()
        self.nbp = list()
        self.mean = list()
        self.median = list()
        self.std = list()

    def add(self, bp, nbp, mean, median, std):
        self.bp.append(bp)
        self.nbp.append(nbp)
        self.mean.append(mean)
        self.median.append(median)
        self.std.append(std)

    def values(self):
        return (np.sum(self.bp), np.sum(self.nbp), np.mean(self.mean), np.mean(self.median), np.mean(self.std))

list_chrs, list_sizes = get_chromosomes_names(bam)

worker = container()

print("Found {} chromosomes to count".format(len(list_chrs)))

for chr, size in zip(list_chrs, list_sizes):
    (bp, nbp, mean, median, std) = count_depth(chr, size, threshold, bam)
    worker.add(bp, nbp, mean, median, std)
    print("Finished with {} {} {} {} {} {}".format(chr, bp, nbp, mean, median, std))

with open(output, 'w') as final:
    final.write(sname + "\t" + "\t".join(worker.values()) + "\n")
