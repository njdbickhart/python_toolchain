import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from collections import defaultdict
import pandas
import numpy as np
import re

asms = snakemake.params["asms"]
print(asms)

buscos = snakemake.input
print(list(buscos))

data = defaultdict(list)

for i, b in enumerate(list(buscos)):
    with open(b, 'r') as input:
        for l in b:
            l = l.strip()
            if l.startswith('#'):
                continue
            elif l.startswith('C'):
                m = re.match(r'C:.+%\[S:(.+)%,D:(.+)%\],F:(.+)%,M:(.+)%,n:.+', l)
                data["Assembly"].append(asms[i])
                data["CompleteSC"].append(float(m.group(1)))
                data["CompleteDup"].append(float(m.group(2)))
                data["Fragmented"].append(float(m.group(3)))
                data["Missing"].append(float(m.group(4)))
                break

df = pandas.DataFrame(data)
#df.set_index("Assembly")

fig, ax = plt.subplots()

ax = df[["CompleteSC", "CompleteDup", "Fragmented", "Missing"]].plot.hbar(stacked=True, edgecolor='none')
ax.legend(bbox_to_anchor=(1.03, 1.0))

plt.xticks(np.arange(len(asms)), df["Assembly"])
plt.savefig(snakemake.output)
