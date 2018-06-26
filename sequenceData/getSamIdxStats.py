# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:03:40 2018

@author: dbickhart
"""

import argparse
import os.path
import sys
import pysam
import contextlib

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program designed to tabulate information from singular or multiple bam files"
            )
    parser.add_argument('-b', '--bams', nargs='+',
                        help="BAM files separated by single spaces",
                        required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Output tab file with stats",
                        required=False, type=str, default="-"
                        )
    
    return parser.parse_args()

def main(args):
    # Fail-fast check for file existence
    for b in args.bams:
        if not os.path.isfile(b):
            print("Could not find file: {}!".format(str(b)))
            sys.exit()
    
    data = {}
    with smartOut(args.output) as out:
        for b in args.bams:
            data[b] = bamStats()
            for line in pysam.idxstats(b).split('\n'):
                segs = line.split('\t')
                if len(segs) < 4:
                    continue
                data[b].addChr(segs[0], int(segs[1]), int(segs[2]), int(segs[3]))
            
        out.write("BamName\tTotalReads\tMappedReads\tUnmappedReads\tMapProportion\tRawXcov\tMapXcov\tavgRawChrCov\tavgMapChrCov\n")
        for b in data:
            v = data[b]
            v.calculateStats()
            out.write("{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n"\
                      .format(b, v.totalReads, v.mappedReads, v.unmappedReads, v.mapProportion, v.rawXcov, v.mapXcov, v.avgRawChrCov, v.avgMapChrCov))
            
@contextlib.contextmanager            
def smartOut(filename : str = None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
        
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

class bamStats:
    def __init__(self):
        self.chrMapped = {}
        self.chrUnmapped = {}
        self.chrLen = {}
        self.totalReads = 0
        self.mappedReads = 0
        self.unmappedReads = 0
        self.mapProportion = 0.0
        self.rawXcov = 0
        self.mapXcov = 0
        self.avgRawChrCov = 0
        self.avgMapChrCov = 0
    
    def addChr(self, chrom : str, length : int, mapped : int, unmapped : int):
        self.chrLen[chrom] = length
        self.chrMapped[chrom] = mapped
        self.chrUnmapped[chrom] = unmapped
        
    def calculateStats(self):
        genomeLength = 0
        for k, v in self.chrLen.items():
            genomeLength += v
            
        mapChrXcov = []
        rawChrXcov = []
        numchrs = 0
        for k in self.chrMapped:
            v = self.chrMapped[k]
            self.totalReads += v
            numchrs += 1
            self.mappedReads += v
            if self.chrLen[k] <= 0:
                continue
            mapChrXcov.append(v / self.chrLen[k])
            
        numchrs -= 1
        for k in self.chrUnmapped:
            v = self.chrUnmapped[k]
            self.totalReads += v
            self.unmappedReads += v
            if self.chrLen[k] <= 0:
                continue
            rawChrXcov.append((v + self.chrMapped[k]) / self.chrLen[k])
            
        self.rawXcov = self.totalReads / genomeLength
        self.mapXcov = self.mappedReads / genomeLength
        self.mapProportion = self.mappedReads / self.totalReads
        self.avgRawChrCov = self.average(rawChrXcov)
        self.avgMapChrCov = self.average(mapChrXcov)
        
    def average(self, values):
        if not values:
            return 0.0
        s = 0
        for v in values:
            s += v
        return s / len(values)

if __name__ == '__main__':
    args = parse_user_input()
    main(args)
