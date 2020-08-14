import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from collections import defaultdict
import pandas
import numpy as np

maxlines = snakemake.params["maxlines"]

inputfile = snakemake.input["associations"]
outputplot = snakemake.output["lineplot"]

# Organize the data into a pandas DataFrame
data = defaultdict(list)
with open(inputfile, 'r') as input:
    head = input.readline()
    for l in input:
        segs = l.rstrip().split()
        data["name"].append(segs[0])
        # Make sure that each line has a value up to the max lines
        lastcount = 0
        for i, val in enumerate(segs[2:]):
            bsegs = segs[i + 2].split(';')
            if i + 1 > maxlines:
                data['{}'.format(maxlines)][-1] += int(bsegs[1])
            else:
                data['{}'.format(str(i+1))].append(int(bsegs[1]))
                lastcount += 1
        if lastcount < maxlines - 1:
            # fill in empty values
            for i in range(maxlines - lastcount):
                data['{}'.format(str(i+1))].append(0)

# Generate  a list of column headers for plotting
colors = [ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000',
          '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
          '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ]
plotcols = list()
for i in range(maxlines):
    plotcols.append('{}'.format(str(i + 1)))

for k, v in data.items():
    print(f'{k}\t{len(v)}')

df = pandas.DataFrame(data)

# Sort the dataframe by reverse order
df = df.sort_values(by=plotcols, ascending=[False for x in plotcols])

print(df.head())

# Now plot that linegraph!
fig = plt.figure(figsize=(6,8))
ax = fig.add_subplot(111)

for i, v in enumerate(plotcols):
    ax.plot('name', v, data=df, marker='', color=colors[i], linewidth=2)

ax.set_xticklabels([])
ax.axis('tight')
ax.legend(loc='upper right')
plt.savefig(outputplot)
