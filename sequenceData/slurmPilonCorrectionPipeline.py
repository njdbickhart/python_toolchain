# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:41:48 2018

@author: dbickhart
"""

import argparse
import os
import subprocess
from slurmAlignScriptBWA import slurmTools

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A pipeline for running separate pilon tasks on a slurm cluster"
            )
    parser.add_argument('-f', '--frag', 
                        help="Input fragment bam file",
                        required=True, type=str)
    parser.add_argument('-g', '--genome', 
                        help="Input fasta file [must have been used in the alignment!]",
                        required=True, type=str)
    parser.add_argument('-o', '--output', 
                        help="Base output directory",
                        required=True, type=str)
    parser.add_argument('-p', '--partition',
                        help="[Optional] Slurm partition for job run",
                        required=False, type=str, default="general")
    parser.add_argument('-e', '--memory',
                        help="[Optional] Job memory in megabytes",
                        required=False, type=int, default=10000)
   
    parser.add_argument('-t', '--time',
                        help="[Optional] Time limit for jobs",
                        required=False, type=str, default="None")
    parser.add_argument('-q', '--qos',
                        help="[Optional] Slurm partition for job run",
                        required=False, type=str, default="general")
    parser.add_argument('-m', '--meta',
                        help="[Optional] Add multiple targets for metagenome assembly",
                        action='store_true') 
    return parser.parse_args()

def main(args):
    modules = ['pilon/1.23', 'samtools/1.4.1']
    
    curDir = os.getcwd()
    
    chrLens = {}
    # Checking chromosome lengths
    if not os.path.isfile(args.genome + ".fai"):
        print("Did not find index for fasta file ({}). Generating now...".format(args.genome))
        subprocess.check_output('module load samtools/1.4.1 && samtools faidx {}'.format(args.genome), 
                                shell=True)
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        
    
    with open(args.genome + ".fai", "r") as fh:
        for l in fh:
            l = l.rstrip('\n')
            segs = l.split()
            chrLens[segs[0]] = segs[1]
    
    worker = slurmTools(curDir,
                        curDir + "/" + args.output + "/scripts",
                        curDir + "/" + args.output + "/outLog",
                        curDir + "/" + args.output + "/errLog",
                        False,
                        modules,
                        1, 3, args.memory, args.time, args.partition, args.qos)
    
    # Main loop. Separate chr and values and queue up Pilon jobs
    if args.meta :
	# Queue up smaller contig jobs and queue them in duplicate sets according to a job submission limit
        temp = {}
        array = []
        count = 0
        scount = 0
        ilimit = int(round((len(chrLens.items()) / 1000 / 10) + 0.5))
        print('Ilimit = {}'.format(ilimit))
        worker.mem = args.memory
        for k, l in chrLens.items():
            temp[k] = l
            if len(temp.keys()) >= 10:
                targets = ','.join(temp.keys())
                temp = {}
                gmem = int(args.memory / 1000)
                cmd = 'java -Xmx{}G -jar $PILON_HOME/pilon-1.23.jar --genome {} --frags {} --output {}.pilon --outdir {} --fix indels --targets {} --verbose --nostrays'.format(str(gmem), args.genome, args.frag, "meta." + str(count), args.output, targets)
                #worker.createGenericCmd(cmd, "pilon_" + str(count))
                array.append(cmd)
                count += 1
                if count % ilimit == 0 and count != 0:
                    print('Queuing script: {}'.format(scount))
                    worker.createArrayCmd(array, "pilon_" + str(scount))
                    array = []
                    scount += 1
        if array or temp:
            if temp:
                targets = ','.join(temp.keys())
                cmd = 'java -Xmx10G -jar $PILON_HOME/pilon-1.23.jar --genome {} --frags {} --output {}.pilon --outdir {} --fix indels --targets {} --verbose --nostrays'.format(args.genome, args.frag, "meta." + str(count), args.output, targets)
            array.append(cmd)
            worker.createArrayCmd(array, "pilon_" + str(count))
    else:
        for k, l in chrLens.items():
            mem = int((int(l) / 1000) * 1.15)
            mem = 8000 if mem < 8000 else mem
            worker.mem = mem
        
            cmd = 'java -Xmx{}M -jar $PILON_HOME/pilon-1.23.jar --genome {} --frags {} --output {}.pilon --outdir {} --fix indels --targets {} --verbose --nostrays'.format(int(mem), args.genome, args.frag, k, args.output, k)
            worker.createGenericCmd(cmd, "pilon" + k)
    
    print('Queuing jobs...')    
    ids = worker.queueJobs()
    print('Queued a total of {} jobs'.format(len(ids)))

if __name__ == "__main__":
    args = parse_user_input()
    main(args)
