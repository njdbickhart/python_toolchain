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

    n25 = np.quantile(bases, 0.25)
    n75 = np.quantile(bases, 0.75)

    min = min(bases)
    max = max(bases)
    mean = np.mean(bases) if count > 0 else 0
    median = np.median(bases) if count > 0 else 0
    values, vcounts = np.unique(np.digitize(bases, [1, 15, 30, max]), return_counts = True)
    vdict = dict(zip(values, vcounts))
    nZbp = count - vdict[0]
    hmCells = [vdict[x] for x in range(1,4)]
    stdev = np.std(bases) if count > 0 else 0
    return (bp, nZbp, mean, n25, median, n75, max, stdev, hmCells[0], hmCells[1], hmCells[2])

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
        self.nZbp = list()
        self.mean = list()
        self.n25 = list()
        self.median = list()
        self.n75 = list()
        self.max = list()
        self.std = list()
        self.hcell1 = list()
        self.hcell2 = list()
        self.hcell3 = list()

    def add(self, bp, nZbp, mean, n25, median, n75, max, std, hcell1, hcell2, hcell3):
        self.bp.append(bp)
        self.nZbp.append(nZbp)
        self.mean.append(mean)
        self.n25.append(n25)
        self.median.append(median)
        self.n75.append(n75)
        self.max.append(max)
        self.std.append(std)
        self.hcell1.append(hcell1)
        self.hcell2.append(hcell2)
        self.hcell3.append(hcell3)

    def values(self):
        return (np.sum(self.bp), np.sum(self.nZbp), np.mean(self.mean),
        np.mean(self.n25), np.mean(self.median), np.mean(self.n75), np.mean(self.max), np.mean(self.std),
        np.sum(self.hcell1), np.sum(self.hcell2), np.sum(self.hcell3))

list_chrs, list_sizes = get_chromosomes_names(bam)

worker = container()

print("Found {} chromosomes to count".format(len(list_chrs)))

for chr, size in zip(list_chrs, list_sizes):
    (bp, nZbp, mean, n25, median, n75, max, stdev, hcell1, hcell2, hcell3) = count_depth(chr, size, threshold, bam)
    worker.add(bp, nZbp, mean, n25, median, n75, max, stdev, hcell1, hcell2, hcell3)
    print("Finished with {} {} {} {} {} {} {} {} {}".format(chr, bp, nZbp, mean, median, stdev, hcell1, hcell2, hcell3))

with open(output, 'w') as final:
    final.write(sname + "\t" + "\t".join(worker.values()) + "\n")
