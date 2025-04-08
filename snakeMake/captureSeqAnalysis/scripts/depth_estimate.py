import os
import sys
import pysam
import numpy as np
from collections import defaultdict

usage = f'python {sys.argv[0]} <input bamfile> <coordinates for calibrator region> <output depth bedfile>'

if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)

def make_chromosome_histogram(input):

    # opening the bam file with pysam
    bamfile = pysam.AlignmentFile(input, 'rb')
    # query all the names of  the chromosomes in a list
    list_chromosomes = bamfile.references
    list_length = bamfile.lengths
    bamfile.close()

    # Create the histogram
    histo = defaultdict(list) # chr -> [list = length of 10kb segments]
    for chr, len in zip(list_chromosomes, list_length):
        for i in range(0, len, 10000):
            histo[chr].append(-1)
    return histo


def count_depth(region, input):
    """
    Count the depth of the read. For each genomic coordinate return the
    number of reads
    -----
    Parameters :
        region: the region of the calibrator sequence
        input: the bamfile to run on the calculation
    -----
    Returns :
        int : count of pileups above threshold
    """
    bp = list()
    bamfile = pysam.AlignmentFile(input, 'rb')
    for pileupcolumn in bamfile.pileup(region=region):
        depth = pileupcolumn.nsegments
        bp.append(depth)
    bamfile.close()
    return np.median(bp)


bam = sys.argv[1]

median = count_depth(sys.argv[2], bam)

print(f'Total bases: {median}')
with open(sys.argv[3], 'w') as final:
    final.write(f'{median}\n')
