import os
import sys
import pysam
import time

def get_chromosomes_names(input):

    # opening the bam file with pysam
    bamfile = pysam.AlignmentFile(input, 'rb')
    # query all the names of  the chromosomes in a list
    list_chromosomes = bamfile.references
    list_length = bamfile.lengths
    bamfile.close()
    return list_chromosomes, list_length


def count_depth(chr_name, size, threshold, input):
    """
    Count the depth of the read. For each genomic coordinate return the
    number of reads
    -----
    Parameters :
        chr : (str) name of the chromosome
        threshold : (int) minimum value to count pileup
    -----
    Returns :
        int : count of pileups above threshold
    """
    bp = 0
    bamfile = pysam.AlignmentFile(input, 'rb')
    for pileupcolumn in bamfile.pileup(chr_name):
        depth = pileupcolumn.nsegments
        if depth >= threshold:
            bp += 1
    bamfile.close()
    return bp


bam = snakemake.input["samples"]
threshold = snakemake.params["threshold"]
print(f'Starting depth estimate for bam: {bam} at threshold {threshold}')
sum = 0

list_chrs, list_sizes = get_chromosomes_names(bam)

print("Found {} chromosomes to count".format(len(list_chrs)))

for chr, size in zip(list_chrs, list_sizes):
    sum += count_depth(chr, size, threshold, bam)
    print(f'Finished with chr: {chr}. {size} {sum}')

print(f'Total bases: {sum}')
with open(snakemake.output["samdepth"], 'w') as final:
    final.write(f'{sum}\n')
