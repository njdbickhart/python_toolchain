#!/usr/bin/env python3
import sys, getopt
import glob
from itertools import islice

usage = "fastqStats.py -b <folder with fastq files>"

def main(argv):
	basefolder = ''
	try:
		opts, leftover = getopt.getopt(argv, "hb:")
	except getopt.GetoptError:
		print("Error parsing options!", usage)
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print(usage)
			sys.exit()
		elif opt == '-b':
			basefolder = arg
			
	if not basefolder:
		print(usage)
		sys.exit()
		
	# Get a list of fastq files in the directory
	files = []
	for exts in ('*.fq', '*.fastq'):
		files.extend(glob.glob('{}/{}'.format(basefolder, exts)))
		
	print("Total files:\t{}".format(len(files)))
	#files = [glob.glob('{}/{}'.format(basefolder, e)) for e in ['*.fq', '*.fastq']]
	total_reads = total_bases = largest =0
	for f in files:
		counts = parse_fastq(f)
		total_reads += counts[0]
		total_bases += counts[1]
		if counts[2] > largest:
			largest = counts[2]
	
	# Print out final data
	avg = total_bases / total_reads
	print("Total reads:\t{}".format(total_reads))
	print("Total bases:\t{}".format(total_bases))
	print("Avg bases:\t{}".format(avg))
	print("Largest read:\t{}".format(largest))

def parse_fastq(file):
	"""This function returns the read count and base count of a fastq file
	Tuple1: read_count
	Tuple2: base_count
	Tuple3: largest_read
	"""
	read_count = base_count = largest = 0
	with open(file, 'r') as f:
		while True:
			lines = []
			eof = False
			for i in range(4):
				l = f.readline()
				if not l:
					eof = True
					break
				lines.append(l)
			if eof: break
			read_count += 1
			base_count += len(lines[1])
			if base_count > largest:
				largest = base_count
	return (read_count, base_count, largest)
	
if __name__ == "__main__":
	main(sys.argv[1:])