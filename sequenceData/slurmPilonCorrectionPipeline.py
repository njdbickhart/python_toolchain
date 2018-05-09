# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:41:48 2018

@author: dbickhart
"""

import argparse
import os
import subprocess
from slurAlignScriptBWA import slurmTools

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
    
    return parser.parse_args()

def main(args):
    modules = ['pilon/1.22', 'samtools/1.4.1']
    
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
                        1, 2, 8000, -1, args.partition)
    
    # Main loop. Separate chr and values and queue up Pilon jobs
    for k, l in chrLens.items():
        mem = int((l / 1000) * 1.15)
        mem = 8000 if mem < 8000 else mem
        worker.mem(mem)
        
        cmd = 'java -Xmx{}M -jar $PILON_HOME --genome {} --frags {} --output {}.pilon --outdir {} --fix bases --targets {} --verbose --nostrays'.format(mem, args.genome, args.frag, k, args.output, k)
        worker.createGenericCmd(cmd, "pilon" + k)
        
    ids = worker.queueJobs()
    print('Queued a total of {} jobs'.format(len(ids)))

if __name__ == "__main__":
    args = parse_user_input()
    main(args)