#!/usr/bin/python3
"""
This is a script designed to process a large bam file and generate statistics on read mapping error rates
The goal is to estimate error rates of long reads mapped against a reference genome (useful for Nanopore data)
"""

import sys
import re
import math
import subprocess
#from multiprocessing import Pool
USAGE = "Calculate edit distance error rate for aligned reads.\n python3 <input bam file>"

# Compiled regular expression to ensure capture of MDZ tag info
NMITAG = re.compile("NM:i")
NMIVAL = re.compile("NM:i:(\d+)")

"""
For each read, count the number of NMI tag differences and return a 
ratio against length of read
"""
def MapDiff(L): 
    segs = L.decode("utf-8").split("\t")
    if int(segs[1]) & 2048 == 2048:
        return -1
        
    if len(segs) < 12:
        return -1
    mo = re.match(NMIVAL, segs[11]);
    if mo:
        val = mo.group(1)
    
        seqlen = len(segs[9])
        return int(val) / seqlen
    else:
        return -1
    
"""
Reduce all ratios by taking the average and calculating the stdev
Returns the [avg, stdev]
"""
def ReduceAvgStd(L):
    num = len(L)
    
    c = 0
    high = 0
    low = 1
    for i in L:
        if i > high:
            high = i
        if i < low:
            low = i
        c += i
        
    avg= 0
    if num > 0:
        avg = c / num
    else:
        avg = 0
        
    # Now for the stdev calculation
    ss = 0
    for i in L:
        ss += (i - avg) * (i - avg)
    
    stdev = 0
    if num > 0: 
        stdev = math.sqrt(ss / num)
    else:
        stdev = 0
    
    return [num, high, low, avg, stdev]

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print(USAGE)
        sys.exit(1)
    
    f = subprocess.Popen(['samtools', 'view', sys.argv[1]], \
                        stdout=subprocess.PIPE)
    [num, high, low, avg, stdev] = ReduceAvgStd(list(filter(lambda x: x != -1, map(MapDiff, f.stdout.readlines()))))
    
    # Uncomment the following line if you want to get the list of ratios
    """
    with open("ratio_count.list", "w") as o:
        for l in list(filter(lambda x: x != -1, map(MapDiff, f.stdout.readlines()))):
            o.write("{}\n".format(l))
    """
        
    print("BAM file:{}\tNum reads: {}\tError High: {}\tError Low: {}\tAvg Error: {}\tStdev Error: {}"\
         .format(sys.argv[1], num, high, low, avg, stdev))
    