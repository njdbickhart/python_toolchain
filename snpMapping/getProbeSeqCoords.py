#!/usr/bin/python3
"""
This script is designed to process Infinium probe alignments to a reference genome
Inputs are a sam file alignment
Output is the tab delimited list of coordinate positions

@author: dbickhart
"""

import sys, getopt
import re

usage = "python3 getProbeSeqCoords.py -s <sam file> -o <output tab file> -c <1 bp correction flag>"
oddSubstitutions = {"[T/A]", "[A/T]", "[G/C]", "[C/G]"}
baseSub = re.compile("(\[.\/.\])")
cigarRe = re.compile("(\d+)([MIDNSHP=X])")

def main(argv):
    samFile = outTab = ''
    correction = False
    try:
        opts, leftover = getopt.getopt(argv, "hcs:o:")
    except getopt.GetoptError:
        print("Error processing command line arguments!\n{}".format(usage))
        sys.exit(2)
        
    for opt, arg in opts:
        if(opt == "-h"):
            print(usage)
            sys.exit(0)
        elif(opt == "-s"):
            samFile = arg
        elif(opt == "-o"):
            outTab = arg
        elif(opt == "-c"):
            correction = True
        
    if(not samFile or not outTab):
        print(usage)
        sys.exit(0)
        
    with open(samFile, "r") as f, \
        open(outTab, "w") as o:
        for line in f:
            if(re.match("^\@", line)):
                continue
            else:
                coord = coordinate_converter(line, correction)
                if(coord.isValid()):
                    vals = coord.getVals()
                    e = int(vals[1]) + 1
                    o.write("{}\t{}\t{}\t{}\n".format(vals[0], \
                            vals[1], e, vals[2]))
        
    
def coordinate_converter(line, correction):
    """
    Input is a sam file line
    Output is coord class object
    """    
    line.rstrip('\n')
    segs = line.split("\t")
    
    # Get substitution base
    variant = baseSub.search(segs[0]).group(1)
    addition = 0 if variant in oddSubstitutions else 1
    if(not correction):
        addition = 0
    
    # Get raw SNP name without variant base
    name = re.sub("_\[.\/.\]", "", segs[0])
    if(int(segs[1]) & 2048 == 2048):
        # secondary alignment -- skip
        return coord()
    elif(int(segs[1]) & 16 == 16):
        # reverse alignment, keep the 5' alignment position
        return coord(segs[2], int(segs[3]) - addition, name)
    elif(int(segs[1]) & 4 == 4):
        # Unaligned segment
        return coord(segs[2], 0, name)
    else:
        # forward alignment, count the alignment length
        return coord(segs[2], int(segs[3]) + int(cigar_conversion(segs[5])) + addition - 1, name)
    
def cigar_conversion(cigar):
    end = 0
    for m in cigarRe.finditer(cigar):
        count, base = m.groups()
        if(base == 'M' or base == 'I' or base == 'S' or base == '=' or base == 'X'):
            end += int(count)
    return end
    
class coord:    
    def __init__(self, chrom = None, pos = None, name = None):
        self.chrom = chrom
        self.pos = pos
        self.name = name
        
        if(self.chrom is None):
            self.valid = False
        else:
            self.valid = True
        
    def isValid(self):
        return self.valid
    
    def getVals(self):
        return (self.chrom, self.pos, self.name)
        
    
if __name__ == "__main__":
    main(sys.argv[1:])
