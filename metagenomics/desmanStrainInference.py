# -*- coding: utf-8 -*-
"""
DESMAN pipeline
Created on Tue Aug 21 15:19:58 2018

@author: dbickhart
"""

import argparse
import os
from typing import Tuple
import subprocess as sp


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
    
    return parser.parse_args()

def main(args):
    # Create output directories
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
        print(f'Created output folder: {args.output}')
        
    os.chdir(args.output)
    
    print("Running FindEliteGenes")
    elites, species = findEliteGenes(args.desman, args.contigs, args.assembly)
    
    print("Running ElitePileups")
    pileup = elitePileups(args.bam, elites, args.assembly)
    
    print("Running CallEliteVariants")
    freq_var, freq_df = callEliteVariants(args.desman, pileup, args.assembly)
    
    print("Running estimateStrainCountDesman")
    fits = estimateStrainCountDesman(args.desman, freq_var, freq_df)
    
    print("Plotting strain replicate data")
    plotDev(args.desman, fits)
        
def findEliteGenes(desman : str, contigs : str, assembly : str) -> Tuple[str, str]:
    cmd = ["python", desman + "/scripts/extract_species_contigs.py", 
           assembly, contigs]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, stdout=open("species_contigs.fa", "w"))
    
    cmd = [desman + "/external/phylosift_v1.0.1/phylosift", "search", "--besthit",
           "--isolate", "species_contigs.fa"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=True)
    
    cmd = [desman + "/external/phylosift_v1.0.1/phylosift", "align", "--besthit",
           "--isolate", "species_contigs.fa"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=True)
    
    cmd = ["python", desman + "/scripts/get_elite_range.py", 
           "PS_temp/species_contigs.fa/blastDir/lookup_ID.1.tbl", 
           "PS_temp/species_contigs.fa/alignDir/DNGNGWU*.codon.updated.1.fasta"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, stdout=open("elites.bed", "w"))
    
    return ["elites.bed", "species_contigs.fa"]

def elitePileups(bam : str, elites : str, assembly : str) -> str:
    cmd = ["samtools", "mpileup", "-l",
           elites, "-f", assembly, bam ]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, stdout=open("elite.pileup", "w"))
    
    return "elite.pileup"

def callEliteVariants(desman : str, pileup : str, assembly : str) -> Tuple[str, str]:
    cmd = ["python", desman + "/scripts/pileups_to_freq_table.py", 
           assembly, pileup, "desmanfreqs.csv"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=True)
    
    cmd = ["python", desman + "/desman/Variant_Filter.py", 
           "desmanfreqs.csv", "-o", "dfreqs", "-p"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=True)
    
    return ["dfreqssel_var.csv", "dfreqstran_df.csv"]

def estimateStrainCountDesman(desman : str, freq_var : str, freq_df : str) -> str:
    fits = []
    for g in range(1,10):
        for repid in range(1,5):
            cmd = [desman + "/bin/desman", freq_var, "-e", freq_df, "-o",
                   f'cluster_{g}_{repid}', "-r", "1000", "-i", "100", "-g", str(g), "-s",
                   str(repid)]
            print(f'Running strain count round {g} and replicate {repid}')
            sp.run(cmd, stdout=open(f'cluster_{g}_{repid}.out', "w"))
            
            sp.run(["cp", "*/fit.txt", f'fit_{g}_{repid}.txt'], shell=True)
            fits.append(f'fit_{g}_{repid}.txt')
            
    with open("desman_dic.fits", 'w') as out:
        out.write("H,G,LP,Dev\n")
        for x in fits:
            with open(x, 'r') as fh:
                next(fh) # skipping header
                for lines in fh:
                    lines = lines.rstrip()
                    out.write(lines + "\n")
                    
    return "desman_dic.fits"

def plotDev(desman : str, alldics : str) -> None:
    cmd = [desman + "/scripts/PlotDev.R", "-l", alldics, "-o", "Dev.pdf"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=True)
    
    print(f'Completed strain testing! Plots are in Dev.pdf')

if __name__ == "__main__":
    args = parse_user_input()
    main(args)