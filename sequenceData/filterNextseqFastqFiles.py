#!/usr/bin/env python3
# This is a script designed to filter and remove "G-base" reads from a pair of fastq files
# The filtration is simple -- if one or both read pairs are comprised entirely of "G-bases" then 
# they are filtered in the output file

import sys, getopt
import re

usage = "filterNextseqFastqFiles.py -f <first read> -r <second read> -o <output base filename>";

def main(argv):
	file1 = file2 = output = ''
	try:
		opts, leftover = getopt.getopt(argv, "hf:r:o:")
	except getopt.GetoptError:
		print("Error parsing arguments!", usage)
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == '-h':
			print(usage)
			sys.exit()
		elif opt == '-f':
			file1 = arg
		elif opt == '-r':
			file2 = arg
		elif opt == '-o':
			output = arg
	
	if (not file1) or (not file2) or (not output):
		print(usage)
		sys.exit()
		
	# Parse both files simultaneously
	results = parse_two_fastq(file1, file2, output)
	
	print("Processed {} reads and identified {} single G-base and {} both G-base artifacts".format(results['reads'], results['oneRead'], results['bothRead']))

def parse_two_fastq(file1, file2, output):
	"""
	This function takes two fastq files, processes them four lines
	at a time and removes the line if it failes a filtration check
	"""
	
	gbasePattern = re.compile("^[GgNn]+$")
	reads = oneRead = bothRead = 0
	
	with open(file1, "r") as f, open(file2, "r") as r, \
			open(output + ".1.fq", "w") as o1, \
			open(output + ".2.fq", "w") as o2:
		while True:
			try:
				forward = get_four_lines(f)
				revread = get_four_lines(r)
				
				reads += 1
				
				if(not forward[0] or not revread[0]):
					break
				
				# If Gbase read pattern is found, then don't print the read
				if gbasePattern.match(forward[1]) or gbasePattern.match(revread[1]):
					if gbasePattern.match(forward[1]) and gbasePattern.match(revread[1]):
						bothRead += 1
					else:
						oneRead += 1
					continue
				
				o1.write("".join(map(str, forward)) + "\n")
				o2.write("".join(map(str, revread)) + "\n")
			except EofException:
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
		
if __name__ == "__main__":
	main(sys.argv[1:])