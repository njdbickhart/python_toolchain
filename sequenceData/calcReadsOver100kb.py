import os
import sys
import numpy as np
import contextlib

def fastq_reader_fh(infile):
    name = infile.readline().rstrip()
    while True:
        seq = ""
        for s in infile:
            if s[0] == '+':
                commentp = s.rstrip()
                break
            else:
                seq += s.rstrip()
        qual = ""
        for q in infile:
            if len(qual) > 0 and  q[0] == '@':
                yield name, seq, qual
                name = q.rstrip()
                break
            else:
                qual += q.rstrip()
        else:
            yield name, seq, qual
            return

@contextlib.contextmanager
def smartFile(filename : str, mode : str = 'r'):
    if filename == 'stdin' or filename == 'stdout':
        if filename == 'stdin':
            fh = sys.stdin
        else:
            fh = sys.stdout
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename != 'stdin' and filename != 'stdout':
            fh.close()

if len(sys.argv) < 2:
    print("Usage = python3 calcReadsOver100kb.py <input fastq file or 'stdin' [can use wildcards for bulk processing]>")
    # If the number of input arguments is less than 1, exit the program
    sys.exit()

lengths = list()
# Above is the list that contains the lengths of the reads
lenover100 = list()
# Above is the list that contains the reads > 100kb in size
for filename in sys.argv[1:]:
    with smartFile(filename, 'r') as input:
        for name, seq, qual in fastq_reader_fh(input):
            lengths.append(len(seq))
            # we use the fastq_reader_fh subroutine to get just the sequence
            if len(seq) >= 100000:
                lenover100.append(len(seq))
                # add the length of the sequence (len function) to the lists
                # if it's greater than 100kb, then add it to the special list

totLen = len(lengths)
intLen = len(lenover100)
# These are easy to calculate. Just use the length function to get the sizes of the arrays

# Next, we use the numpy library functions to calculate the median and sums of the lengths
totMedian = np.median(lengths)
totSum = np.sum(lengths)
intSum = np.sum(lenover100)

# Finally, we calculate the X coverage against an approximation of the size of the cattle genome
totXCov = totSum / 2800000000
intXCov = intSum / 2800000000

print("TotalReadNum\tReadsGt100kbNum\tTotalReadLenMedian\tTotalReadBases\tReadsGt100kbBases\tTotalReadXCov(cattle)\tReadsGt100kbXcov(cattle)")
print("\t".join([str(x) for x in [totLen,intLen,totMedian,totSum,intSum,totXCov,intXCov]]))
