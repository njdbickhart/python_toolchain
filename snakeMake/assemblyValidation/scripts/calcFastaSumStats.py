# This script is designed to run on a snakemake pipeline.
import numpy as np
from snakemake.shell import shell
#from collections import defaultdict

#data = defaultdict(list)
f = snakemake.input[0]
fai = f + '.fai'
shell("samtools faidx {f}")

with open(snakemake.output[0], 'w') as out:
    out.write("Num\tSumBp\tMeanBp\tMedianBp\tStdevBp\tContigN50\n")

    sizes = list()

    # Now read file and store fasta length lists
    with open(fai, "r") as input:
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

    out.write(f'{count}\t{sum}\t{mean}\t{median}\t{stdev}\t{n50}\n')
