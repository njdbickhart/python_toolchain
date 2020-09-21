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

alist = snakemake.input["asms"]
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
        temp.append(temp[i-1] + (g.iat[i, 1] / asmsize * 100))
    NGX.extend(temp)

df = df.assign(NGX = NGX)

colors = [ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
          '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
          '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ]

print(df.head())
# Plot the lines
fig, ax = plt.subplots()
i = 0
for k, g in df.groupby(['asm']):
    #ax = g.plot(ax=ax, kind='line', x='NGX', y='len', c=colors[i], label=k)
    ax = g.plot(ax=ax, marker='', x='NGX', y='len', c=colors[i], linewidth=1, label=k)
    i += 1

ax.vlines(x=50.0, ymin=0, ymax=largestctg, linestyles='dashed')
plt.legend(loc='best')
plt.savefig(snakemake.output["plot"])
