#!/usr/bin/python3
"""
Created on Tue May  8 13:32:19 2018

@author: dbickhart
"""

import argparse
import os, sys
import random
import subprocess
import string
from typing import List
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A pipeline for aligning sequence data on a slurm cluster"
            )
    parser.add_argument('-b', '--base', 
                        help="Base output folder name",
                        required=True, type=str
                        )
    parser.add_argument('-t', '--tab',
                        help="Input tab file with fastq file locations and library names",
                        required=True, type=str
                        )
    parser.add_argument('-f', '--fasta',
                        help="Input reference fasta file for alignment",
                        required=True, type=str
                        )
    parser.add_argument('-m', '--merge',
                        help="[Optional] queue alignments and merge scripts",
                        action='store_true'
                        )
    parser.add_argument('-q', '--qos',
                        help="[Optional] name of the slurm qos for jobs",
                        required=False, type=str, default="general"
                        )
    parser.add_argument('-p', '--partition',
                        help="[Optional] name of the slurm partition for jobs",
                        required=False, type=str, default="general"
                        )
    parser.add_argument('-e', '--time',
                        help="[Optional] Time limit on jobs",
                        required=False, type=str, default="None"
                        )
    return parser.parse_args()

def main(args):
    modules = ['bwa/0.7.12', 'samtools/1.4.1']
    
    scriptCount = 0
    slurmWorkers = {} # Each sample gets its own worker
    slurmBams = defaultdict(list) # Each sample gets a list of sorted bams
    
    if not os.path.exists(args.base):
        os.makedirs(args.base)
        
    curDir = os.getcwd()
    if not os.path.isfile(args.fasta):
        if not os.path.isfile(curDir + "/" + args.fasta):
            print("Error finding fasta file location {}!\n".format(args.fasta))
            sys.exit()
    
    # Open tab file and generate scripts for each sample
    with open(args.tab, "r") as infile:
        for l in infile:
            l = l.rstrip()
            segs = l.split()
            
            wkdir = curDir + "/" + args.base + "/" + segs[-1]
            if not segs[-1] in slurmWorkers:
                slurmWorkers[segs[-1]] = slurmTools(
                        wkdir, 
                        wkdir + "/scripts",
                        wkdir + "/outLog",
                        wkdir + "/errLog",
                        False,
                        modules,
                        1, 10, 25000, args.time, args.partition, args.qos)
                
            bname = os.path.basename(segs[0])
            bsegs = bname.split('.')
            fasta = curDir + args.fasta if os.path.isfile(curDir + args.fasta) else args.fasta
            
            uname = bsegs[0] + "." + urlHash()
            
            cmd = "bwa mem -t 8 -M -R '@RG\\tID:{LB}\\tSM:{ID}\\tLB:{LB}' {FA} {seg1} {seg2} | samtools sort -m 2G -o {uname}.sorted.bam -T {uname} -".format(ID=segs[-1], LB=segs[-2], FA=fasta, seg1=segs[0], seg2=segs[1], uname=uname)
            slurmBams[segs[-1]].append(uname + ".sorted.bam")
            
            slurmWorkers[segs[-1]].createGenericCmd(cmd, "bwaAlign")
            scriptCount += 1
    
    # TODO: generate queuing mechanism
    if args.merge:
        for k, worker in slurmWorkers.items():
            jobIds = worker.queueJobs()
            
            print("Sample: {} queued with {} jobs: {}".format(k, len(jobIds),
                  ' '.join(jobIds)))
            wkdir = curDir + "/" + args.base + "/" + k
            merger = slurmTools(
                    wkdir,
                    wkdir + "/scripts",
                    wkdir + "/outLog",
                    wkdir + "/errLog",
                    False,
                    modules,
                    1, 7, 9000, args.time, args.partition, args.qos)
            
            bams = slurmBams[k]
            cmds = []
            
            for b in bams:
                cmds.append('samtools index {}'.format(b))
                
            cmds.append("samtools merge -c -p -@ 6 {}.sorted.merged.bam {}".format(k, ' '.join(bams)))
            cmds.append("samtools index {}.sorted.merged.bam".format(k))
            
            merger.addJobIDs(jobIds)
            merger.createArrayCmd(cmds, "samMerger")
            merger.queueJobs()
            print("Queued jobs with dependencies")

def urlHash():
    #alphabet = list([string.ascii_lowercase, string.digits])
    collections = {}
    for x in range(1, 30):
        length = random.randrange(0, 11) + 10
        key = ''.join(random.choice(list(string.ascii_lowercase)) for x in range(length))
        collections[key] = 1
    return list(collections.keys())[0]
    
class slurmTools:
    def __init__(self, workDir = ".", scriptDir = "./scripts", outDir = "./out",
                 errDir = "./err", useTime = False, modules = [], nodes = 1,
                 tasks = 1, mem = 100, time = "None", partition = "general", qos = "general"):
        self.workDir = workDir
        self.scriptDir = scriptDir
        self.outDir = outDir
        self.errDir = errDir
        self.useTime = useTime
        self.modules = modules
        self.nodes = nodes
        self.tasks = tasks
        self.mem = mem
        self.time = time
        self.partition = partition
        self.qos = qos
        
        # Containers for later processes
        self.scripts = []
        self.dependencies = []
        self.jobIds = []

    def addJobIDs(self, jarray: List[int]) -> None:
        for x in jarray:
            self.dependencies.append(x)

    def queueJobs(self):
        for x in self.scripts:
            stdout = subprocess.getoutput('sbatch {}'.format(x))
            self.jobIds.append(stdout.split()[-1])

        return self.jobIds
        
    def createArrayCmd(self, carray: List[str], sname: str = "script") -> None:
        self._generateFolders()
        
        hashv = self._generateHash()
        sname += "_" + hashv + ".sh"
        
        header = self._generateHeader(sname)
        
        for x in carray:
            header += 'echo "{}"\n{}\n\n'.format(x, x)
        header += "wait\n"
        
        fulldir = self.scriptDir + "/" + sname
        with open(fulldir, "w") as o:
            o.write(header)
            
        self.scripts.append(fulldir)
        
    def createGenericCmd(self, cmd: str, sname: str = "script") -> None:
        self._generateFolders()
        
        hashv = self._generateHash()
        sname += "_" + hashv + ".sh"
        
        header = self._generateHeader(sname)
        header += 'echo "{}"\n'.format(cmd)
        header += '{}\nwait\n'.format(cmd)
        
        fulldir = self.scriptDir + "/" + sname
        with open(fulldir, "w") as o:
            o.write(header)
            
        self.scripts.append(fulldir)
        
    def _generateFolders(self) -> None:
        for x in [self.errDir, self.scriptDir, self.outDir, self.workDir]:
            if not os.path.exists(x):
                os.makedirs(x)
                
    def _generateHash(self) -> str:
        random_set = []
        seen = {}
        
        for x in range(1, 5):
            candidate = int(random.random() * 1185)
            if(candidate in seen):
                candidate = int(random.random() * 1185)
            seen[candidate] = 1
            random_set.append(candidate)
            
        return ''.join(str(x) for x in random_set)
    
    def _generateHeader(self, sName: str) -> str:
        tag = '#SBATCH'
        fields = []
        
        fields.append('#!/usr/bin/bash')
        fields.append('{} --nodes={}'.format(tag, self.nodes))
        fields.append('{} --ntasks-per-node={}'.format(tag, self.tasks))
        fields.append('{} --mem={}'.format(tag, self.mem))
        fields.append('{} --output={}/{}_%j.out'.format(tag, self.outDir, sName))
        fields.append('{} --error={}/{}_%j.err'.format(tag, self.errDir, sName))
        fields.append('{} --workdir={}'.format(tag, self.workDir))
        fields.append('{} --partition={}'.format(tag, self.partition))
        fields.append('{} -q {}'.format(tag, self.qos))

        if self.time != "None":
            fields.append('{} -t {}'.format(tag, self.time))
        
        if self.dependencies:
            depStr = "afterany"
            for x in self.dependencies:
                depStr += ":" + x
            fields.append('{} --dependency={}'.format(tag, depStr))
            
        fields.append("\n")
        
        if self.modules:
            fields.append('module load {}'.format(' '.join(self.modules)))
            
        fields.append("\n")
        return '\n'.join(fields)
            
    
if __name__ == '__main__':
    args = parse_user_input();
    main(args)
