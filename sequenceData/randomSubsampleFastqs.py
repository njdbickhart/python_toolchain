# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:39:38 2019

@author: dbickhart
"""

import argparse
import gzip
import random
import contextlib

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Subsample pairs of fastq files for approximating a numerical count of reads"
            )
    parser.add_argument('-f', '--first', 
                        help="The first read pair file. Can be specified more than once, but files must be paired in order!",
                        action='append', required=True
                        )
    parser.add_argument('-s', '--second',
                        help="The second read pair file. Can be specified more than once, but must be paired with a \"first\" file!",
                        action='append', required=True
                        )
    parser.add_argument('-o', '--output',
                        help="The output base name for the paired files.",
                        type=str, required=True
                        )
    parser.add_argument('-l', '--lines',
                        help="The number of expected, randomly sampled lines in the output files",
                        type=int, required=True
                        )
    parser.add_argument('-i', '--log',
                        help="The logfile for the edits",
                        type=str, default="subsample.log"
                        )
    
    return parser.parse_args()

def main(args):
    # It sucks, but we need to open the files and count all of the lines. 
    # We can store this information in a log
    log = open(args.log, 'w')
    first = args.first
    second = args.second
    
    # The number of reads should be the same between pairs. We'll assume to save time
    lcount = 0
    for f in first:
        temp = lineCounter(f)
        log.write(f'{temp}\tlines in file\t{f}\n')
        lcount = lcount + temp
        
    # We set the count to a divisor of four to account for the reads
    selection_indicies = random.sample(range(lcount / 4), args.lines)
    selectionSet = set(selection_indicies)
    
    # OK, this is where things get complex. I'm going to keep looping through files
    current_read=0
    with open(args.output + "_R1.fq", 'w') as out1, open(args.output + "_R2.fq", 'w') as out2:
        for f1, f2 in zip(first, second):
            with gzipFile(f1, 'r') as fh1, gzipFile(f2, 'r') as fh2:
                try:
                    f1lines = get_four_lines(fh1)
                    f2lines = get_four_lines(fh2)
                
                    if current_read in selectionSet:
                        out1.write('\n'.join(f1lines) + '\n')
                        out2.write('\n'.join(f2lines) + '\n')
                    
                    current_read += 1
                
                    if current_read % (args.lines / 10):
                        log.write(f'!Currently {current_read / args.lines} through sampling on files: {f1} {f2}\n')
                except EofException as e:
                    log.write(f'!Reached EOF of files {f1} and {f2}\n')
                
    log.write(f'!Finished sampling\n')
    
    log.close()
        
    
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

@contextlib.contextmanager
def gzipFile(filename : str, mode : str = 'r'):
    """
    This is not as refined as my magic number evaluation in java,
    but I couldn't be bothered to write up something more fancy in python!
    sue me :)
    """
    if filename.endswith('\.gz'):
        fh = gzip.open(filename, mode)
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        fh.close()
        
def lineCounter(file : str) -> int:
    with gzipFile(file, "r") as fh:
        return sum(bl.count("\n") for bl in blocks(fh))

def get_four_lines(filehandle):
	"""
	This function pulls four lines at a time and returns a four element
	tuple for parsing
	"""
	
	lines = []
	for i in range(4):
		l = filehandle.readline()
		l = l.rstrip()
		if not l:
			raise EofException("l.rstrip", "Reached end of file!")
		lines.append(l)
	return lines

class EofException(Exception):
	""" Exception for reaching end of file
	
	Attributes:
		expression -- message from error
		message -- explanation of message
	"""
	
	def __init__(self, expression, message):
		self.expression = expression
		self.message = message
    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)