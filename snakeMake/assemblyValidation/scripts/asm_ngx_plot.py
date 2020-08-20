import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from itertools import cycle
from collections import defaultdict
import pandas
import numpy as np

asms = snakemake.params["asms"]
print(asms)

alist = snakemake.input
print(alist)

data = defaultdict(list)

for i, v in enumerate(alist):
    atext = asms[i]
    with open(v + ".fai") as input:
        for l in input:
            s = l.rstrip().split()
            data['asm'].append(atext)
            data['len'].append(int(s[1]))

df = pandas.DataFrame(data)
df = df.sort_values(by=['asm', 'len'], ascending=(True, False))
# Lazy implementation: assuming largest assembly is correct here!
# TODO: allow user input for max assembly size instead
largestctg = df['len'].max()
asmsize = 0
NGX = list()
for k, g in df.groupby(['asm']):
    if g['len'].sum() > asmsize:
        asmsize = g['len'].sum()
    temp = [0.0]
    for i in range(1,len(g)):
        temp.append(NGX[i-1] + (g.at[i, 'len'] / asmsize * 100))
    NGX.extend(temp)

df.assign(NGX = NGX)

# Plot the lines
fig, ax = plt.subplots()
for k, g in df.groupby(['asm']):
    ax = g.plot(ax=ax, kind='line', x='NGX', y='len', c=k, label=k)

ax.vlines(x=50.0, ymin=0, ymax=largestctg, linestyles='dashed')
plt.legend(loc='best')
plt.savefig(snakemake.output["plot"])
