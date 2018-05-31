# -*- coding: utf-8 -*-
"""
Created on Thu May 31 10:16:58 2018

@author: dbickhart
"""

import argparse
import gzip

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A program designed to process SRA fastq files and make them compatible with BWA"
            )
    parser.add_argument('-f', '--forward', 
                        help="forward read file",
                        required=True, type=str
                        )
    parser.add_argument('-r', '--reverse',
                        help="reverse read file",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="base name for the modified output file",
                        required=True, type=str
                        )
    parser.add_argument('-l', '--log',
                        help="Save program stdout to a logfile",
                        required=True, type=str
                        )
    
    return parser.parse_args()

def main(args):
    parse_two_fastq(args.forward, args.reverse, args.output, args.log)

def parse_two_fastq(file1, file2, output, log):
    """
    This function takes two fastq files, processes them four lines
    at a time and removes the line if it failes a filtration check
    """
	
    reads = oneRead = bothRead = 0
	
    with open(file1, "r") as f, open(file2, "r") as r, \
        gzip.open(output + ".1.fq.gz", "wb") as o1, \
        open(log, "w") as log, \
        gzip.open(output + ".2.fq.gz", "wb") as o2:
        log.write("Beginning the parsing of files for: " + output)
        while True:
            try:
                forward = get_four_lines(f)
                revread = get_four_lines(r)
				
                # Remove last two characters of first segment string
                forward[0] = forward[0].split()[:-2]
                revread[0] = revread[0].split()[:-2]
                reads += 1
				
                # Replacing "+" with blank placeholder
                forward[2] = "+"
                revread[2] = "+"
                if(not forward[0] or not revread[0]):
                    raise EofException("not forward or not revread", "Unexpected termination of one fastq file!")
				
                				
                o1.write("".join(map(str, forward)) + "\n")
                o2.write("".join(map(str, revread)) + "\n")
            except EofException as e:
                log.write("Completed {} reads and exited with code: {}".format(reads, str(e)))
                break
    return {"reads" : reads, "oneRead" : oneRead, "bothRead" : bothRead}

def get_four_lines(filehandle):
	"""
	This function pulls four lines at a time and returns a four element
	tuple for parsing
	"""
	
	lines = []
	for i in range(4):
		l = filehandle.readline()
		l.rstrip("\n")
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
        
if __name__ == '__main__':
    args = parse_user_input()
    main(args)