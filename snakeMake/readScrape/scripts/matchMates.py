# This script is designed to run with a snakemake pipeline
# It takes unmapped mate alignments, makes sure they map to contig ends
# and gives a bedpe file with the likely insertion site of the novel sequence
# It also prints only the contigs that fit the filtering criteria
from collections import defaultdict
import re

log = open(snakemake.log[0], 'w')
class evidence:
    def __init__(self, scaff):
        self.scaffold = scaff
        self.chrs = defaultdict(int)
        self.chrpos = defaultdict(lambda : defaultdict(int))

    def addEvidence(self, chr, pos):
        self.chrs[chr] += 1
        self.chrpos[chr][pos] += 1

    def dumpAllObs(self, out):
        for k, v in sorted(self.chrs.items(), key=lambda x: x[1], reverse=True):
            out.write(f'{self.scaffold}\t{k}\t{v}')
            for j, l in self.chrpos[k].items():
                out.write(f'\t{j}')
            out.write("\n")

    def getLikelyLoc(self, edge):
        # Select winner
        winner = "None"
        max = 0
        for k, v in self.chrs.items():
            if v > max:
                winner = k
                max = v

        # select all of the positions from the winner
        start = 0
        end = 0
        reads = 0
        # Sort by largest value to start with the most frequently observed
        worker = {k: v for k, v in sorted(self.chrpos[winner].items(), key=lambda x: x[1], reverse = True)}
        for k, v in worker.items():
            if start == 0:
                start = int(k)
                reads += v
                continue

            comp = int(k)
            if comp > start and comp <= start + edge:
                if comp > end:
                    end = comp
                    reads += v
            elif comp < start and comp >= start - edge:
                if start > end:
                    end = start
                    start = comp
                else:
                    start = comp
                reads += v

        # In the rare case of only one observance
        if end == 0:
            end = start
        # return bed with chr, start, end, novelcontig, support
        return (winner, str(start), str(end), self.scaffold, str(reads))


class read:
    def __init__(self):
        self.name = "None"
        self.num = 0
        self.chr = "None"
        self.pos = 0

    def loadSam(self, segs):
        self.num = 1 if int(segs[1]) & 0x40 else 2
        self.name = segs[0]
        self.chr = segs[2]
        self.pos = int(segs[3])

    def loadBed(self, segs):
        self.num = 1 if segs[3].endswith('/1') else 2
        self.name = segs[3]
        self.chr = segs[0]
        self.pos = int(segs[1])

unlinks = snakemake.input["unmaplinks"]
links = snakemake.input["links"]
fai = snakemake.input["rfai"]
contigs = snakemake.input["fasm"]
edge = int(snakemake.params["edge"])
sampname = snakemake.params["samp"]
raw = snakemake.output["raw"]

# Read in scaffold fai file for lengths
ctglens = dict()
with open(fai, 'r') as input:
    for l in input:
        s = l.rstrip().split()
        ctglens[s[0]] = int(s[1])

log.write("Finished reading fai file\n")

# Read unmapped read names and create data structure
unlist = dict()
filt = 0
total = 0
with open(unlinks, 'r') as input, open(snakemake.output["filtreads"], 'w') as out:
    for l in input:
        s = l.rstrip().split()
        # Test if read maps to ends of contig
        clen = ctglens[s[0]]
        pos = int(s[1])
        total += 1
        if pos > edge and pos < (clen - edge):
            filt += 1
        else:
            out.write(l)
            r = read()
            r.loadBed(s)
            unlist[s[3]] = r

log.write(f'UnmappedReads: Filtered {filt} reads out of {total}\n')

# Read links sam and try to associate with mates
evid = dict()   #Evidence class dictionary # {scaffold} -> evidence
keeps = dict()  #Scaffolds to keep
with open(links, 'r') as input:
    for l in input:
        if l.startswith('@'):
            continue
        s = l.rstrip().split()

        if s[0] in unlist:
            r = read()
            r.loadSam(s)
            comp = unlist[s[0]]
            if comp.num != r.num:
                if comp.chr not in evid:
                    evid[comp.chr] = evidence(comp.chr)
                evid[comp.chr].addEvidence(r.chr, r.pos)
                keeps[comp.chr] = 1

log.write("Finished collecting evidence of mate mapping\n")

# produce bed file with mate map associations
count = 0
with open(raw, 'w') as temp, open(snakemake.output["flanks"], 'w') as out:
    for k, v in evid.items():
        v.dumpAllObs(temp)
        temp = v.getLikelyLoc(edge)
        count += 1
        out.write('\t'.join(temp) + '\n')

log.write(f'Wrote out {count} entries with evidence\n')

# Filter out scaffolds that had evidence
filtered = 0
with open(contigs, 'r') as input, open(snakemake.output["mcontigs"], 'w') as output:
    for l in input:
        if l.startswith(">"):
            lsegs = l.split()
            lm = re.match(r">(.+)", lsegs[0]).group(1)

            if lm in keeps:
                write = True
            else:
                write = False
                filtered += 1
        if write: output.write(l)

log.write(f'Filtered {filtered} novel contigs at the end\n')
log.close()
