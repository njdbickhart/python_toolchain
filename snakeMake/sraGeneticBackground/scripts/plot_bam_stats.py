import os
import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

usage = f'python3 {sys.argv[0]} <input tab summary> <output pdf filename>'

if len(sys.argv) != 2:
    print(usage)
    sys.exit(-1)

cols = ["SNAME", "ThreshBp", "Zbp", "Mean", "Q25", "Median", "Q75", "Stdev"]

data = defaultdict(list)
with open(sys.argv[1], 'r') as input:
    head = l.readline().rstrip().split()
    for l in input:
        s = l.rstrip().split()
        for i, k in enumerate(s):
            data[head[i]].append(k)

df = pd.DataFrame(data)


(fig, axis) = plt.subplots(ncols=1, nrows=4)
axs = axis.flatten().tolist()
axs.reverse()

ax = axs.pop()
ax.set_title("Count of bases under coverage threshold", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Bases")
df[["SNAME", "Zbp"]].plot(x="SNAME", y="Zbp", ax=ax)

ax = axs.pop()
ax.set_title("Interquartile coverage values", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Coverage")
df[["SNAME", "Q25", "Median", "Q75"]].plot(x= "SNAME", ax=ax)

ax = axs.pop()
ax.set_title("Coverage Standard Deviation", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Stdev")
df[["SNAME", "Stdev"]].plot(x="SNAME", y="Stdev", ax=ax)

plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)
fig.set_size_inches(21,12)
fig.tight_layout()
plt.subplots_adjust(top=0.9)

plt.savefig(sys.argv[2])

print('Fini')
