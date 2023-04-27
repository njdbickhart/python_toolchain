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

if len(sys.argv) != 3:
    print(usage)
    sys.exit(-1)

cols = ["SNAME", "ThreshBp", "Zbp", "Mean", "Q25", "Median", "Q75", "Max", "Stdev",
'sub15', 'sub30', 'gt30']

data = defaultdict(list)
with open(sys.argv[1], 'r') as input:
    head = input.readline().rstrip().split()
    for l in input:
        s = l.rstrip().split()
        for i, k in enumerate(s):
            data[head[i]].append(k)

df = pd.DataFrame(data)
df["Zbp"] = df["Zbp"].astype('int32')
for i in ("Q25", "Median", "Q75", "Stdev"):
    df[i] = df[i].astype('float')

print(df)

(fig, axis) = plt.subplots(ncols=1, nrows=4)
axs = axis.flatten().tolist()
axs.reverse()

ax = axs.pop()
ax.set_title("Count of bases under coverage threshold", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Bases")
df[["SNAME", "Zbp"]].plot(kind='bar', x="SNAME", y="Zbp", ax=ax)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')

ax = axs.pop()
ax.set_title("Interquartile coverage values", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Coverage")
boxes = [{'label' : r['SNAME'], 'whislo' : 0, 'q1' : r['Q25'], 'med' : r['Median'],
'q3' : r['Q75'], 'whishi' : r['Max'], 'fliers' : []} for i, r in df.iterrows()]
#df[["SNAME", "Q25", "Median", "Q75"]].plot(x= "SNAME", ax=ax)
ax.bxp(boxes, showfliers=False)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')

ax = axs.pop()
ax.set_title("Coverage Standard Deviation", loc="left")
ax.set_xlabel("Samples")
ax.set_ylabel("Stdev")
df[["SNAME", "Stdev"]].plot(kind='bar', x="SNAME", y="Stdev", ax=ax)
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')

ax = axs.pop()
ax.set_xlabel("Samples")
tdf = df[["SNAME", 'Zbp', 'sub15', 'sub30', 'gt30']].set_index('SNAME').copy().T
sns.heatmap(tdf, ax=ax)

plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)
fig.set_size_inches(21,12)
fig.tight_layout()
plt.subplots_adjust(top=0.9)

plt.savefig(sys.argv[2])

print('Fini')
