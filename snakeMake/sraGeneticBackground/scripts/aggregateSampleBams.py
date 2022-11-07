import sys
import os
import json
import logging
import glob
from collections import defaultdict
from snakemake import shell

usage = f'python3 {sys.argv[0]} config samplename fastqdir reference asmname logfile'
if len(sys.argv) != 7:
    print(usage)
    sys.exit(-1)

# log file
logging.basicConfig(filename=sys.argv[6], encoding='utf-8', level=logging.DEBUG)

# load configfile
config = json.load(open(sys.argv[1], 'r'))

# Get the sample name
sampleToSRA = defaultdict(list)
for k, v in config["samples"].items():
    sampleToSRA[v[1]].append(v[0])
sranames = sampleToSRA.get(sys.argv[2], [])

if len(sranames) == 0:
    logging.error(f'Could not identify any SRA names for sample: {sys.argv[2]}')
    sys.exit(-1)

# Get the SRA files needed for alignment
fastqs = defaultdict(list)
nfiles = 0
for i in sranames:
    tfiles = glob.glob(f'{sys.argv[3]}/{i}*.fastq')
    fastqs[i].append(tfiles)
    nfiles += len(tfiles)

logging.info(f'Found {nfiles} for sample: {sys.argv[2]}')
logging.info(f'Beginning alignment')

os.makedirs(f'mapped/{sys.argv[5]}', exist_ok=True)

mergefiles = []
# Begin alignment
for k, l in fastqs.items():
    tfile = f'mapped/{sys.argv[5]}/{k}.temp.bam'
    cmd = f'bwa mem -t 8 {sys.argv[4]} {l[0]} {l[1]} 2>> {sys.argv[6]}' +
    f'| samtools sort - > {tfile} 2>> {sys.argv[6]}'
    logging.info(f'CMD: {cmd}')
    mergefiles.append(tfile)

    shell(cmd)
    cmd = f'samtools index {tfile}'
    logging.info(f'CMD: {cmd}')

    shell(cmd)

# If merger, then run the following
outfile = f'mapped/{sys.argv[5]}/{sys.argv[2]}.merged.bam'
if len(mergefiles) > 1:
    fstring = " ".join(mergefiles)
    cmd = f'samtools merge -@ 8 {outfile} {fstring}; samtools index {outfile}'
    logging.info(f'CMD: {cmd}')

    shell(cmd)

    # Delete temporary files if finished
    if os.path.exists(outfile):
        for i in mergefiles:
            cmd = f'rm {i}'
            logging.info(f'Deleting {cmd}')
            shell(cmd)
else:
    cmd = f'mv {mergefiles[0]} {outfile}; samtools index {outfile}'
    logging.info(f'Only one file, so changing name\nCMD: {cmd}')

    shell(cmd)



logging.info("Finished!")
