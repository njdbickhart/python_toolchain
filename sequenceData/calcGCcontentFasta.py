# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:40:22 2017

@author: dbickhart
"""

import argparse
import concurrent.futures
import collections

def gc_task(name, seq):
    count = count_gc(seq)
    return name, count

def count_gc(seq):
    gc = 0
    for i in seq:
        if i == "G" or i == "C" or i == "g" or i == "c":
            gc += 1
    return gc / float(len(seq))

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith('>'):
            line = line.replace('>', '')
            if name: yield(name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield(name, ''.join(seq))

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program to calculate GC content from large fasta entries"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input fasta file containing long (> 100 bp) sequences",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output tab delimited file with fasta names and avg GC content",
                        required=True, type=str
                        )
    parser.add_argument('-t', '--threads',
                        help="Number of threads to use",
                        required=False, type=int, default=1
                        )
    
    return parser.parse_args()


if __name__ == '__main__':
    
    args = parse_user_input()
    fasta = args.fasta
    output = args.output
    
    threads = args.threads
    results = {}
        
    fp = open(fasta, "r")
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for name, seq in read_fasta(fp):
            count = executor.submit(count_gc, seq)
            results[name] = count
    
    fp.close()
    
    sortResults = collections.OrderedDict(sorted(results.items()))
    with open(output, "w") as out:
        for name, count in sortResults.items():
            c = 0.0
            try:
                c = count.result()
            except Exception as ex:
                print("{} generated an exception: {}".format(name, ex))
            out.write("{}\t{}\n".format(name, c))
            
    print("Wrote GC content to output")
