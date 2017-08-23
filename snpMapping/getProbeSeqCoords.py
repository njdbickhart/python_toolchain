#!/usr/bin/python3
"""
This script is designed to process Infinium probe alignments to a reference genome
Inputs are a sam file alignment
Output is the tab delimited list of coordinate positions

@author: dbickhart
"""

import sys, getopt
import re

usage = "python3 getProbeSeqCoords.py -s <sam file> -o <output tab file>"
oddSubstitutions = {"[T/A]", "[A/T]", "[G/C]", "[C/G]"}
baseSub = re.compile("(\[.\/.\])")
cigarRe = re.compile("(\d+)([MIDNSHP=X])")

def main(argv):
    samFile = outTab = ''
    try:
        opts, leftover = getopt.getopt(argv, "hs:o:")
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
        
    if(not samFile or not outTab):
        print(usage)
        sys.exit(0)
        
    with open(samFile, "r") as f, \
        open(outTab, "w") as o:
        for line in f:
            if(re.match("^\@", line)):
                continue
            else:
                coord = coordinate_converter(line)
                if(coord.isValid()):
                    e = int(coord.getVals[1]) + 1
                    o.write("{}\t{}\t{}\t{}", coord.getVals[0], \
                            coord.getVals[1], e, coord.getVals[2])
        
    
def coordinate_converter(line):
    """
    Input is a sam file line
    Output is coord class object
    """    
    line.rstrip('\n')
    segs = line.split("\t")
    
    # Get substitution base
    variant = baseSub.search(segs[0]).group(1)
    addition = 0 if variant in oddSubstitutions else 1
    
    if(int(segs[1]) & 2048 == 2048):
        # secondary alignment -- skip
        return coord()
    elif(int(segs[1]) & 16 == 16):
        # reverse alignment, keep the 5' alignment position
        return coord(segs[2], int(segs[3]) - addition, segs[0])
    else:
        # forward alignment, count the alignment length
        return coord(segs[2], int(segs[3]) + cigar_conversion(segs[5]) + addition, segs[0])
    
def cigar_conversion(cigar):
    end = 0
    for m in cigarRe.finditer(cigar):
        count, base = m.groups()
        if(base == 'M' or base == 'I' or base == 'S' or base == '=' or base == 'X'):
            end += count
    return count
    
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
        
    
if __name__ == "main":
    main(sys.argv[1:])