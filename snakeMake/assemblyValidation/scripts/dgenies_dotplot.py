import os
import sys
from dgenies.lib.paf import Paf
from dgenies.bin.index import Index, index_file
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
import numpy as np

plotdots = snakemake.params["pixels"]

paths=[]
idxs=[]
for i in (snakemake.input["query"], snakemake.input["reference"]):
    dir = os.path.dirname(i)
    file = os.path.basename(i)
    print(f'Indexing {dir}/{file}...')
    (success, numctgs, err) = index_file(dir, file, file + ".idx")
    if not success:
        print(err)
        sys.exit(-1)
    else:
        paths.append(dir)
        idxs.append(file)

# Load data structures

paf_file = snakemake.input["paf"]
idx1 = os.path.join(paths[0], idxs[0])
idx2 = os.path.join(paths[1], idxs[1])
paf = Paf(paf_file, idx1, idx2, False)
paf.sort()

# Calculate values for matrix
asize1 = map(sum, Index.load(idx1)[2])
asize2 = map(sum, Index.load(idx2)[2])

awinsize1 = asize1 / plotdots
awinsize2 = asize2 / plotdots
