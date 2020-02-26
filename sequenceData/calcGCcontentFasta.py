# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:40:22 2017

@author: dbickhart
Version 2: Added fastq support and read/contig lengths to columns
"""
import sys
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
    
def read_fastq(fp):
    name = fp.readline().rstrip()
    while True:
        seq = ""
        for s in fp:
            if s[0] == '+':
                break
            else:
                seq += s.rstrip()
        qual = ""
        for q in fp:
            if len(qual) > 0 and  q[0] == '@':
                yield name, seq, qual
                name = q.rstrip()
                break
            else:
                qual += q.rstrip()
        else:
            yield name, seq, qual
            return

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program to calculate GC content from large fasta entries"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input fasta file containing long (> 100 bp) sequences",
                        type=str, default=None
                        )
    parser.add_argument('-q', '--fastq',
                        help="Input Fastq file containing long (> 100 bp) sequences",
                        type=str, default=None
                        )
    parser.add_argument('-o', '--output',
                        help="Output tab delimited file with fasta names, seq length and avg GC content, in that order",
                        required=True, type=str
                        )
    parser.add_argument('-t', '--threads',
                        help="Number of threads to use",
                        required=False, type=int, default=1
                        )
    
    return parser, parser.parse_args()


if __name__ == '__main__':
    
    parser, args = parse_user_input()
    if args.fasta == None and args.fastq == None:
        print("At least one file must be specified!")
        print(parser.print_help())
        sys.exit(-1)
    
    output = args.output
    
    threads = args.threads
    results = {}
    lens = {}
        
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        if args.fasta != None:
            with open(args.fasta, 'r') as fp:
                for name, seq in read_fasta(fp):
                    count = executor.submit(count_gc, seq)
                    lens[name] = len(seq)
                    results[name] = count
        elif args.fastq != None:
            with open(args.fastq, 'r') as fp:
                for name, seq, qual in read_fastq(fp):
                    count = executor.submit(count_gc, seq)
                    lens[name] = len(seq)
                    results[name] = count
    
    
    sortResults = collections.OrderedDict(sorted(results.items()))
    with open(output, "w") as out:
        for name, count in sortResults.items():
            c = 0.0
            l = lens[name]
            try:
                c = count.result()
            except Exception as ex:
                print("{} generated an exception: {}".format(name, ex))
            out.write("{}\t{}\t{}\n".format(name, l, c))
            
    print("Wrote GC content to output")
