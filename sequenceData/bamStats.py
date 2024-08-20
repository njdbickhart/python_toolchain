# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 12:54:06 2020

@author: dbickhart
"""

import argparse
import sys
import pysam
import contextlib
import numpy as np
from os.path import basename
from collections import defaultdict


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Calculate and tabulate BAM file index stats"
            )
    parser.add_argument('-f', '--file',
                        help="The input, indexed bam file. May be specified more than once!",
                        action="append", default=[]
                        )
    parser.add_argument('-l', '--list',
                        help="New line delimited list of BAM files [preferred input]",
                        type=str, default="None"
                        )
    parser.add_argument('-o', '--output',
                        help="Output tab delimited file. [stdout]",
                        type=str, default="stdout"
                        )
    parser.add_argument('-r', '--reads',
                        help="How many reads to sample for determining average read length [1000]",
                        type=int, default=1000
                        )

    return parser.parse_args()

def main(args):
    fileHead = ["File", "GenomeSize", "AvgReadLen", "TotalReads", "MappedReads", "MapPerc", "RawXCov", "MapXCov", "AvgMapChrXCov", "AvgChrMapPerc"]

    # Validate input
    if len(args.file) < 1 and args.list == "None":
        print("Error! Must enter at least one file!")
        #args.print_help()
        sys.exit()

    # Process files
    classes = list()
    if args.list == "None":
        for f in args.file:
            worker = samStats(f, args.reads)
            try:
                worker.calculate()
            except BamFileException:
                continue # Continue if the file wasn't a proper bam file!
            classes.append(worker)
    else:
        with open(args.list, 'r') as input:
            for l in input:
                l = l.rstrip().split()[0]
                worker = samStats(l, args.reads)
                try:
                    worker.calculate()
                except BamFileException:
                    continue # Continue if the file wasn't a proper bam file!
                classes.append(worker)

    # Create the list of stats
    data = list()
    for w in classes:
        data.append(w.retList())

    # If the genome size and read length is the same, remove redundant info
    gsizes = set([x[1] for x in data])
    rlen = set([x[2] for x in data])
    cropg = True if len(gsizes) == 1 else False
    cropr = True if len(rlen) == 1 else False

    removes = list()
    if cropg: removes.append(1)
    if cropr: removes.append(2)

    if cropg or cropr:
        fileHead = [x for i, x in enumerate(fileHead) if i not in removes]
        for i, d in enumerate(data):
            d = [x for j, x in enumerate(d) if j not in removes]
            data[i] = d

    # Print output
    with smartFile(args.output, 'w') as out:
        if cropg or cropr:
            print(f'#Genome size: {gsizes} and avg read len: {rlen}. Cropping similar columns')
        out.write('\t'.join(fileHead) + "\n")
        for d in data:
            out.write('\t'.join(d) + "\n")


class samStats:

    def __init__(self, file : str, reads: int):
        self.file = file
        self.reads = reads

        # other attributes
        self.avgRLen = 0
        self.gsize = 0
        self.totReads = 0   # Total reads
        self.mapReads = 0   # Total mapped reads
        self.unmapReads = 0 # Total unmapped reads
        self.compUnmap = 0  # Unmapped, no paired reads
        self.chrLens = dict()
        self.mapChrR = defaultdict(int) # per chr map reads
        self.unmapChrR = defaultdict(int)   # per chr unmap reads

        # Final stats
        self.mapPerc = ""
        self.rawX = ""
        self.mapX = ""
        self.avgMapX = ""
        self.avgMapPerc = ""

    def calculate(self):
        # first determine if file is valid
        bam = pysam.AlignmentFile(self.file, "rb")

        if not bam.check_index():
            raise BamFileException("Index Check", "Bam file did not have index!")

        # Calculate average read length
        self._calcAvgRLen(bam)

        # Fill attribute containers
        self._pullStats()

        # Calculate final summary stats
        self._finalStats()

    def retList(self)-> list:
        return [basename(self.file), "{:,}".format(self.gsize), str(self.avgRLen),
                "{:,}".format(self.totReads), "{:,}".format(self.mapReads),
                self.mapPerc, self.rawX, self.mapX, self.avgMapX, self.avgMapPerc]


    def _finalStats(self):
        self.mapPerc = "{0:.6f}".format(self.mapReads / self.totReads)
        self.rawX = "{0:.6f}".format((self.totReads * self.avgRLen) / self.gsize)
        self.mapX = "{0:.6f}".format((self.mapReads * self.avgRLen) / self.gsize)

        clens = list(self.chrLens.values())
        maps = list(self.mapChrR.values())
        unmaps = list(self.unmapChrR.values())

        avgMX = [((x * self.avgRLen) / y) if y > 0 else 0 for x, y in zip(maps, clens)]
        avgPC = [x / (x + y) if y > 0 or x > 0 else 0 for x, y in zip(maps, unmaps)]
        if len(avgMX) < 1 or len(avgPC) < 1:
            print(f'Error with chr table. avgMx: {avgMX} avgPC: {avgPC}')
            self.avgMapX = "0"
            self.avgMapPerc = "0"
            return
        self.avgMapX = "{0:.6f}".format(np.mean(avgMX))
        self.avgMapPerc = "{0:.6f}".format(np.mean(avgPC))

    def _pullStats(self):
        text = pysam.idxstats(self.file)
        lines = text.split(sep="\n")
        for i in lines:
            segs = i.split()
            if len(segs) < 4:
                continue
            clen = int(segs[1])
            mapped = int(segs[2])
            unmapped = int(segs[3])
            self.totReads += mapped + unmapped
            self.mapReads += mapped
            self.unmapReads += unmapped
            self.gsize += clen

            if segs[0].startswith("*"):
                self.compUnmap = unmapped
            else:
                self.chrLens[segs[0]] = clen
                self.mapChrR[segs[0]] = mapped
                self.unmapChrR[segs[0]] = unmapped


    def _calcAvgRLen(self, bam):
        rlens = []
        for i in bam.head(self.reads):
            l = i.infer_read_length()
            if isinstance(l, int):
                rlens.append(l)

        self.avgRLen = int(np.mean(rlens))


class BamFileException(Exception):

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'r'):
    if filename == 'stdin' or filename == 'stdout':
        if filename == 'stdin':
            fh = sys.stdin
        else:
            fh = sys.stdout
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename != 'stdin' and filename != 'stdout':
            fh.close()

if __name__ == "__main__":
    args = parse_user_input()
    main(args)
