# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:56:15 2018

@author: dbickhart
"""

import argparse
import os
import re
import statsmodels.api as sm
import kneed
import subprocess as sp
from typing import Tuple

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A pipeline for aligning sequence data on a slurm cluster"
            )
    parser.add_argument('-a', '--assembly', 
                        help="Assembly fasta file [full path needed!]",
                        required=True, type=str
                        )
    parser.add_argument('-c', '--contigs',
                        help="Species contig list [full path needed!]",
                        required=True, type=str
                        )
    parser.add_argument('-b', '--bamfile',
                        help="Aligned reads in bam file format [full path needed!]",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="output directory",
                        required=True, type=str
                        )
    parser.add_argument('-d', '--desman',
                        help="Desman home directory",
                        required=True, type=str
                        )
    parser.add_argument('-s', '--strains',
                        help="Desman deviation fit points for strain count determination (usually desman_dic.fits)",
                        required=False, default="desman_dic.fits"
                        )
    
    return parser.parse_args()

def main(args):
    # main routine
    os.chdir(args.output)
    if not os.path.isfile('dfreqstran_df.csv') and not os.path.isfile('dfreqssel_var.csv'):
        print("Missing essential files dfreqstran_df.csv and dfreqssel_var.csv!")
        os.sys.exit(-1)
        
    if not os.path.isfile('species_contigs.fa') and not os.path.isfile('elites.bed'):
        print("Missing essential files elites.bed and species_contigs.fa!")
        os.sys.exit(-1)
        
    # Use LOWESS to smooth data points 
    x = []
    y = []
    with open(args.strains, 'r') as fh:
        next(fh)
        for l in fh:
            segs = re.split(",")
            x.append(int(segs[1]))
            y.append(float(segs[3]))
    
    lowess = sm.nonparametric.lowess(y, x)
    l_x = list(zip(*lowess))[0]
    l_y = list(zip(*lowess))[1]
    
    # Select x value from elbow plot and use that in desman
    k = kneed.KneeLocator(l_x, l_y, invert=True, direction='decreasing')
    strains = k.knee
    
    print(f'Identified {strains} as the optimal strain count from the values')
    with open(args.strains + ".strain.count", 'w') as out:
        out.write(f'{args.output}\t{strains}\n')
    
    # Run desman on the data
    gamma, eta = desmanRun(args.desman, 'dfreqssel_var.csv', 'dfreqstran_df.csv', int(strains))
    
    # Run the strain differentiation 
    strainTigs(args.desman, 'species_contigs.fa', args.assembly, gamma, eta, 
               args.bamfile, 'elites.bed', int(strains))

def desmanRun(desman: str, dfreq_var: str, dfreq_df: str, strains : int) -> Tuple[str, str]:
    cmd = [desman + "/bin/desman", dfreq_var, "-e", dfreq_df, "-o", "cluster_f",
           "-r", "1000", "-i", "100", "-g", strains]
    print(f'CMD: {" ".join(cmd)}')
    sp.run(cmd, shell=False, check=True)
    
    # Moving cluster file
    sp.run("mv cluster/Gamma_star.csv cluster/Eta_star.csv .", shell=True)
    return ["Gamma_star.csv", "Eta_star.csv"]

def strainTigs(desman: str, contigs : str, assembly : str, gamma : str, eta : str,
               bam : str, elite : str, strains : int) -> None:
    print(f'Getting lengths of {contigs}')
    sp.run(f'samtools faidx {contigs}', shell=True)
    
    with open(contigs + ".fai", 'r') as fh:
        for l in fh:
            l = l.rstrip()
            segs = l.split()
            region = f'{segs[0]}:1-{segs[1]}'
            safe = f'{segs[0]}_1_{segs[1]}'
            
            print(f'CMD: samtools mpileup -r {region} -f {assembly} {bam} > s_{safe}.pileup')
            sp.run(f'samtools mpileup -r {region} -f {assembly} {bam} > s_{safe}.pileup', shell=True)
            
    cmd = [desman + "/scripts/pileups_to_freq_table.py", assembly, 's_*.pileup', "contigfreqs.csv"]
    print(f'CMD: {" ".join(cmd)}')
    sp.run(" ".join(cmd), shell=True)
            
    cmd = [desman + "/desman/Variant_Filter.py", "contigfreqs.csv", "-m", "0.0", "-v", "0.03"]
    print(f'CMD: {" ".join(cmd)}')
    sp.run(cmd, check=True)
    
    cmd = [desman + "/scripts/CalcGeneCov.py", "contigfreqs.csv"]
    print(f'CMD: {" ".join(cmd)}')
    sp.run(cmd, stdout=open("contig_cov.csv",'w'))
    
    # Create unique list of contigs with core genes
    ctgs = {}
    with open("elites.bed", 'r') as fh, open("ctgs_core.txt", 'w') as out:
        for l in fh:
            l = l.rstrip()
            segs = l.split()
            ctgs[segs[0]] = 1
    
        for c in ctgs.keys():
            out.write(c + "\n")
            
    cmd = [desman + "/scripts/CalcDelta.py", "contig_cov.csv", "ctgs_core.txt", "cluster_core"]
    print(f'CMD: {" ".join(cmd)}')
    sp.run(cmd, check=True)
    
    # I have little idea what the rest of the pipeline attempts to do, apart from tease out
    # The strains into fastas. Let's stick with the above crop of functions for now
    #cmd = [desman + "/bin/desman", ]
    #print(f'CMD: {" ".join(cmd)}')

if __name__ == "__main__":
    args = parse_user_input()
    main(args)