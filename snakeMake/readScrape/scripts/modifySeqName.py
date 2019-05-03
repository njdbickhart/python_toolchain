# This script is designed to run exclusively with a snakemake pipeline
# It removes the generic scaffold names made by spades and makes them sample-unique
# It also removes fasta entries smaller than 1kb
import subprocess as sp
import re

infile = snakemake.input['asm']
sample = snakemake.params['samp']
ofile = snakemake.output['rasm']

# Check sizes of scaffolds and discard those under 1kb
subprocess.run(f'module load samtools; samtools faidx {infile}', shell=True)
slens = dict()
with open(infile + '.fai', 'r') as input:
    for l in input:
        l = l.rstrip()
        segs = l.split()
        slens[segs[0]] = int(segs[1])

count = 0
write = True
with open(infile, 'r') as input, open(ofile, 'w') as output:
    for l in input:
        if l.startswith(">"):
            lsegs = l.split()
            lm = re.match(r">(.+)")

            # Determine if the contig is > 1000 bp
            sl = slens[lm.group(1)] if lm.group(1) in slens else 1000
            if sl >= 1000:
                write = True
                count += 1
            else:
                write = False

            l = f'>{sample}.scaf.{count}\n'
        output.write(l) if write

subprocess.run(f'module load samtools; samtools faidx {ofile}', shell=True)
