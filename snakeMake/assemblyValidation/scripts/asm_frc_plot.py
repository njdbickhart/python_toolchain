import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from itertools import cycle
from collections import defaultdict
import pandas
import numpy as np


asms = snakemake.wildcards["asm"]
print(asms)

alist = snakemake.input
print(alist)

colors = [ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
          '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
          '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ]

# Start loading assemblies into dataframe
data = defaultdict(list)
for i, v in enumerate(alist):
    atext = asms[i]
    with open(v, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            data['asm'].append(atext)
            data['x'].append(float(s[0]))
            data['y'].append(float(s[1]))

df = pandas.DataFrame(data)

# Plot the FRC data
fig, ax = plt.subplots()
for k, g in df.groupby(['asm']):
    ax = g.plot(ax=ax, kind='line', x='x', y='y', c=k, label=k)

plt.legend(loc='best')
plt.savefig(snakemake.output["plot"])
