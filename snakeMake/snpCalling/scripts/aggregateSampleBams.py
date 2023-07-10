import sys
import os
import json
import logging
import glob
from collections import defaultdict
from snakemake import shell

usage = f'python3 {sys.argv[0]} config samplename reference logfile'
if len(sys.argv) != 5:
    print(usage)
    sys.exit(-1)

# log file
logging.basicConfig(filename=sys.argv[4], encoding='utf-8', level=logging.DEBUG)

# load configfile
config = json.load(open(sys.argv[1], 'r'))

# Get the sample name
sampleFQs = config["samples"].get(sys.argv[2], [])
if len(sampleFQs) == 0:
    logging.error(f'Could not find any fastqs for sample {sys.argv[2]}!')
    sys.exit(-1)


logging.info(f'Found {len(sampleFQs)} pairs for sample: {sys.argv[2]}')
logging.info(f'Beginning alignment')

os.makedirs(f'mapped/', exist_ok=True)

tempfiles = []
mergefiles = []
# Begin alignment
for k, l in enumerate(sampleFQs):
    tfile = f'mapped/{sys.argv[2]}.{k}.temp.bam'
    sfiles = " ".join(l)
    cmd = f'bwa mem -t 8 -R "@RG\\tID:{sys.argv[2]}.{k}\\tSM:{sys.argv[2]}\\tLB:{k}\\tPL:ILLUMINA" {sys.argv[3]} {sfiles} 2>> {sys.argv[4]} | samtools sort - > {tfile} 2>> {sys.argv[4]}'
    logging.info(f'CMD: {cmd}')
    tempfiles.append(tfile)

    shell(cmd)
    cmd = f'samtools index {tfile}'
    logging.info(f'CMD: {cmd}')

    shell(cmd)
    # Mark Duplicates

    dfile = f'mapped/{sys.argv[2]}.{k}.dedup.bam'
    stats = f'mapped/{sys.argv[2]}.{k}.dedup.stats'
    cmd = f'picard MarkDuplicates I={tfile} O={dfile} M={stats} REMOVE_DUPLICATES=true\n'
    cmd += f'samtools index {dfile}'

    logging.info(f'CMD: {cmd}')

    shell(cmd)
    mergefiles.append(dfile)
    tempfiles.append(dfile)
    tempfiles.append(tfile + '.bai')
    tempfiles.append(dfile + '.bai')

# If merger, then run the following
outfile = f'mapped/{sys.argv[2]}.merged.bam'
if len(mergefiles) > 1:
    fstring = " ".join(mergefiles)
    cmd = f'samtools merge -@ 8 {outfile} {fstring}\n samtools index {outfile} 2>> {sys.argv[4]}'
    logging.info(f'CMD: {cmd}')

    shell(cmd)

    # Delete temporary files if finished
    if os.path.exists(outfile):
        for i in tempfiles:
            cmd = f'rm {i}'
            logging.info(f'Deleting {cmd}')
            shell(cmd)
else:
    cmd = f'mv {mergefiles[0]} {outfile}\n'
    cmd += f'samtools index {outfile}'
    logging.info(f'Only one file, so changing name\nCMD: {cmd}')

    shell(cmd)

# Deleting temp files
for f in tempfiles:
    if os.path.exists(f):
        cmd = f'rm {f}'
        logging.info(f'Deleting {f} {cmd}')
        shell(cmd)
    else:
        logging.info(f'Did not delete {f} as it does not exist')

logging.info("Finished!")
