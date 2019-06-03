# This script is designed to run exclusively with a snakemake pipeline
# It removes the generic scaffold names made by spades and makes them sample-unique
# It also removes fasta entries smaller than 1kb
import subprocess as sp
import os
import sys
import re

infile = snakemake.input['asm']
sample = snakemake.params['samp']
ofile = snakemake.output['rasm']

# Check sizes of scaffolds and discard those under 1kb
sp.run(f'module load samtools; samtools faidx {infile}', shell=True)
slens = dict()
with open(infile + '.fai', 'r') as input:
    for l in input:
        l = l.rstrip()
        segs = l.split()
        slens[segs[0]] = int(segs[1])

# Printing to STDOUT to log successful completion
if os.path.isfile(infile + '.fai'):
    print("Successfully indexed scaffolds.")
else:
    print("Failed to index Scaffolds fasta!")
    sys.exit(-1)

count = 0
write = True
with open(infile, 'r') as input, open(ofile, 'w') as output:
    for l in input:
        if l.startswith(">"):
            lsegs = l.split()
            lm = re.match(r">(.+)", lsegs[0])

            # Determine if the contig is > 1000 bp
            sl = slens[lm.group(1)] if lm.group(1) in slens else 1000
            if sl >= 1000:
                write = True
                count += 1
            else:
                write = False

            l = f'>{sample}.scaf.{count}\n'
        if write: output.write(l)

sp.run(f'module load samtools; samtools faidx {ofile}', shell=True)
# Again, another log to confirm completion!
if os.path.isfile(ofile + '.fai'):
    print("Successfully renamed and indexed fasta!")
else:
    print("Job failed at renaming fasta entries!")
    sys.exit(-1)
