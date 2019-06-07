# This script is designed to run on a snakemake pipeline.
import numpy as np
#from collections import defaultdict

#data = defaultdict(list)

with open(snakemake.output[0], 'w') as out:
    out.write("Sample\tCount\tSum\tMean\tMedian\tStdev\tContigN50\n")
    for f in snakemake.input:
        # Get samplename
        segs = f.split("/")
        bsegs = segs[0].split(".")
        sample = bsegs[0]
        sizes = list()

        # Now read file and store fasta length lists
        with open(f, "r") as input:
            for l in input:
                l = l.rstrip()
                lsegs = l.split()
                sizes.append(int(lsegs[1]))

        # Calculate stats and print
        count = length(sizes)
        sum = np.sum(sizes) if count > 0 else 0
        mean = np.mean(sizes) if count > 0 else 0
        median = np.median(sizes) if count > 0 else 0
        stdev = np.std(sizes) if count > 0 else 0

        l50v = int(sum / 2)
        csum = np.cumsum(sizes)
        l50idx = min(csum[csum >= l50v])
        c50idx = np.where(csum == l50idx)

        n50 = sizes[int(c50idx[0])]

        out.write(f'{sample}\t{count}\t{sum}\t{mean}\t{median}\t{stdev}\t{n50}\n')
