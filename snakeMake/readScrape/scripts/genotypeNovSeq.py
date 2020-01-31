# This script is designed to run with a snakemake pipeline
# It reads in cd-hit clusters and chromosome insertion coordinate bed files
# It attempts to genotype each sample and then prepares a summary tab file
# and individual genotype files per sample
from collections import defaultdict
import re
import numpy as np

csub = re.compile(r'[\+\/\-\%]')
crem = re.compile(r'at\s')
csamp = re.compile(r'(.+)\..+\.')
cpat = re.compile(r'(\d+)\s+(\d+)nt, >(.+)\.{3} (.+)$')

log = open(snakemake.log[0], 'w')

class Cluster:
    def __init__(self, cname):
        self.cname = cname
        self.paragon = "None"
        self.sizes = list()
        self.scaffolds = list()
        self.percs = list()
        self.maxSize = 0
        self.samples = None

    def processCDHITLine(self, line):
        line = csub.sub('', line)
        m = cpat.match(line)
        if not m or len(m.groups()) < 4:
            log.write(f'Did not resolve cd-hit line: {line}\n')
            return
        if m.group(4) == "*":
            self.paragon = m.group(3)
            self.percs.append(100.0)
        else:
            temp = crem.sub('', m.group(4))
            self.percs.append(float(temp))
        self.sizes.append(int(m.group(2)))
        self.scaffolds.append(m.group(3))
        if int(m.group(2)) > self.maxSize:
            self.maxSize = int(m.group(2))

    def getScafNums(self):
        return len(self.scaffolds)

    def returnOutTab(self):
        # Cluster name \t paragon \t sizeMax \t size average \t size stdev \t num scaffolds \t avg percID \t min percID \t samples
        aboveOne = True if len(self.sizes) > 1 else False
        sizeAvg = np.average(self.sizes) if aboveOne else self.sizes[0]
        sizeStd = np.std(self.sizes) if aboveOne else 0.0
        percAvg = np.average(self.percs) if aboveOne else self.percs[0]
        percMin = np.amin(self.percs) if aboveOne else self.percs[0]
        self.samples = list()
        for x in self.scaffolds:
            m = csamp.match(x)
            self.samples.append(m.group(1))
        self.samplookup = {self.samples[i]:i for i in range(len(self.samples))}
        return [self.cname, self.paragon, str(self.maxSize), str(sizeAvg), str(sizeStd), str(len(self.scaffolds)), str(percAvg), str(percMin), ';'.join(self.samples)]

    def getSampInfo(self, samp):
        index = self.samplookup[samp]
        return [self.cname, self.scaffolds[index], self.sizes[index], self.percs[index]]

def getSampleName(file):
    s = file.rstrip().split(sep="/")
    d = s[-1].split(sep=".")
    return (d[0], file)

def getScaffIndex(entries, value):
    for i in range(len(entries)):
        if entries[i][1] == value:
            return i
    return -99

cluster = snakemake.input["c_file"]
flanks = {k: v for k, v in map(getSampleName, snakemake.input["flanks"])}
summary = snakemake.output["summary"]

log.write("Starting cluster condensing")

# Generate clusters. Cross ref samples
clist = list()
counter = 0
sampdata = defaultdict(list)
with open(cluster, 'r') as input, open(summary, 'w') as out:
    out.write("cluster\tparagon\tmaxsize\tavgsize\tsizestdev\tnumSamps\tavgPercID\tminPercID\tsampNames\n")
    for l in input:
        l = l.rstrip()
        if l.startswith('>'):
            l = re.sub(r'\s+', '', l)
            cluster = l[1:]
            clist.append(Cluster(cluster))
        else:
            clist[-1].processCDHITLine(l)

    for c in clist:
        counter += 1
        out.write('\t'.join(c.returnOutTab()) + '\n')
        for s in c.samples:
            sampdata[s].append(c.getSampInfo(s))

# Read flanks files and create individual genotype files
association = defaultdict(lambda: defaultdict(int)) # {cluster}->{chr} = count
for k, v in flanks.items():
    with open(v, 'r') as input, open(f'genotypes/{k}.genotype.tab', 'w') as out:
        out.write("sample\tcluster\tscaffname\tlength\tpercID\tchr\tstart\tend\tsupport\n")
        for l in input:
            s = l.rstrip().split()
            working = sampdata[k]
            index = getScaffIndex(working, s[3])
            if index == -99:
                log.write(f'Error finding scaffold: {s[3]} in {k}')
                continue
            temp = [str(x) for x in working[index]]
            association[temp[0]][s[0]] += 1
            out.write(k + "\t" + "\t".join(temp) + f'\t{s[0]}\t{s[1]}\t{s[2]}\t{s[4]}\n')

log.write("Done processing genotypes")
with open(snakemake.output["associations"], 'w') as out:
    out.write("Cluster\tCount\tAssociations\n")
    for k, d in association.items():
        out.write("{}\t{}".format(k, len(d.keys())))
        for j, n in sorted(d.items(), key=lambda item: item[1], reverse=True):
            out.write(f'\t{j};{n}')
        out.write("\n")

log.close()
