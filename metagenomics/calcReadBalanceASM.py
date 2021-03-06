# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 08:41:55 2020

@author: derek.bickhart-adm
"""

import argparse
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pysam
import numpy as np
from collections import defaultdict
from sequenceData.reformatFasta import fastaReader

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to compare read depth profiles for long and short reads on metagenomes"
            )

    parser.add_argument('-o', '--output',
                        help="Output tab delimited table file",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bam',
                        help="Short-read, indexed bam files",
                        action="append", default=[]
                        )
    parser.add_argument('-p', '--paf',
                        help="long-read, aligned paf files",
                        action="append", default=[]
                        )
    parser.add_argument('-l', '--longreads',
                        help="long-read fastas used in alignments",
                        action="append", default=[]
                        )
    return parser.parse_args(), parser
    
class shortRead:
    
    def __init__(self, blist):
        #self.blist = [pysam.AlignmentFile(b, "rb") for b in blist]
        self.blist = blist
        self.rnum = 0
        self.rbases = 0
        self.ctglens = dict()
        self.algnr = defaultdict(int)
        self.algnbases = defaultdict(int)
        
    def _calcAvgRLen(self, bam):
        rlens = []
        b = pysam.AlignmentFile(bam, "rb")
        for i in b.head(100000):
            l = i.infer_read_length()
            if isinstance(l, int):
                rlens.append(l)
            
        return int(np.mean(rlens))
    
    def _pullStats(self, bam, rl):
        text = pysam.idxstats(bam)
        lines = text.split(sep="\n")
        for i in lines:
            segs = i.rstrip().split()
            if len(segs) < 4:
                continue 
            elif segs[0] == "*":
                if segs[0] == "*":
                    self.rnum += int(segs[3])
                    self.rbases += int(int(segs[3]) * rl)
                continue
            clen = int(segs[1])
            mapped = int(segs[2])
            unmapped = int(segs[3])
            
            self.ctglens[segs[0]] = clen
            self.rnum += mapped + unmapped
            self.rbases += int((mapped * rl) + (unmapped * rl))
            self.algnr[segs[0]] += 1
            self.algnbases[segs[0]] += int((mapped * rl))
        
    def calcBases(self):
        # First, collect total number of reads and calculate total readset bases for each bam
        bamRlens = [self._calcAvgRLen(b) for b in self.blist]
        
        # Next, calculate total reads for all bams, total bases and bases mapped
        for i, b in enumerate(self.blist):
            print(f'Estimating bases for {b} with avg read {bamRlens[i]}')
            self._pullStats(b, bamRlens[i])
            
        # Return total reads, total bases, ctglen dict, aligned bases dict and aligned read dict
        return self.rnum, self.rbases, self.ctglens, self.algnbases, self.algnr

class longRead:

    def __init__(self, plist, flist):
        self.plist = plist
        self.flist = flist
        self.rnum = 0
        self.rbases = 0
        self.algnr = defaultdict(int)
        self.algnbases = defaultdict(int)
        
    def calcBases(self):
        # First, collect total number of reads and bases from fastas
        for i in self.flist:
            if not os.path.exists(i + '.fai'):
                print(f'Fasta file: {i} is not samtools indexed! Aborting')
                sys.exit(-1)
            with open(i + '.fai', 'r') as fasta:
                for l in fasta:
                    seg = l.rstrip().split()
                    if len(seg) < 2:
                        continue
                    self.rbases += int(seg[1])
                    self.rnum += 1
        
        print(f'Estimated {self.rnum} reads and {self.rbases} bases from long-read fastas')
        # Next calculate total reads for all paf alignments and bases mapped
        for p in self.plist:
            print(f'Working on paf: {p}')
            with open(p, 'r') as paf:
                for l in paf:
                    s = l.rstrip().split()
                    # TODO: confirm that this column gives the expected values!
                    self.algnbases[s[5]] += int(s[10])
                    self.algnr[s[5]] += 1
                    
        # return total reads, total bases, aligned bases dict and aligned read dict
        return self.rnum, self.rbases, self.algnbases, self.algnr

def main(args, parser):
    # Sanity check
    if len(args.bam) < 1 or len(args.paf) < 1 or (len(args.paf) != len(args.longreads)):
        print("Must enter at least one bam and one paf file, along with an equal number of longread fastas!")
        parser.print_help()
        sys.exit(-1)
        
    sworker = shortRead(args.bam)
    (srnum, srbases, ctglens, salgns, salgnr) = sworker.calcBases()
    print(f'shortreads:\t{srnum}\t{srbases}')
    
    lworker = longRead(args.paf, args.longreads)
    (lrnum, lrbases, lalgns, lalgnr) = lworker.calcBases()
    
    
    print(f'longreads:\t{lrnum}\t{lrbases}')
    
    with open(args.output, 'w') as out:
        out.write("contig\tsreads\tsbases\tsrprop\tsbprop\tlreads\tlbases\tlrprop\tlbprop\n")
        for ctg, clen in ctglens.items():
            sreads = salgnr.get(ctg, 0)
            sbases = salgns.get(ctg, 0)
            srprop = "{0:.3f}".format(sreads / srnum)
            sbprop = "{0:.3f}".format(sbases / srbases)
            
            lreads = lalgnr.get(ctg, 0)
            lbases = lalgns.get(ctg, 0)
            lrprop = "{0:.3f}".format(lreads / lrnum)
            lbprop = "{0:.3f}".format(lbases / lrbases)
            
            out.write(f'{ctg}\t{sreads}\t{sbases}\t{srprop}\t{sbprop}\t{lreads}\t{lbases}\t{lrprop}\t{lbprop}\n')
            
    print("Done!")

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
