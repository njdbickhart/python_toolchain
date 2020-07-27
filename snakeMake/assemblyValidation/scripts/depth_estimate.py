import pysam

sum = 0
count = 0

bam = snakemake.input["samples"]
print("Starting depth esimate")
for l in pysam.depth(bam):
    s = l.rstrip().split()
    if int(s[-1]) >= 3:
        sum += 1
    if count % 500000 == 0:
        print(f'At line: {l}')

with open(snakemake.output["samdepth"], 'w') as final:
    final.write(f'{sum}\n')
