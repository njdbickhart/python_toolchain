#!/software/apps/python_3/gcc/64/3.6.2/bin/python3
#SBATCH --nodes=1
#SBATCH --mem=25000
#SBATCH --ntasks-per-node=8
"""
DESMAN pipeline
Created on Tue Aug 21 15:19:58 2018

@author: dbickhart
"""

import argparse
import os
from typing import Tuple, TextIO, List
import subprocess as sp
import re
import concurrent.futures


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
                        help="Aligned reads in bam file format [full path needed!; Can be specified more than once!]",
                        action="append", default=[]
                        )
    parser.add_argument('-o', '--output',
                        help="output directory",
                        required=True, type=str
                        )
    parser.add_argument('-g', '--genes',
                        help='SCG gene bed file',
                        required =True, type=str
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
    
    #print("Running FindEliteGenes")
    #elites, species = findEliteGenes(args.desman, args.contigs, args.assembly)
    
    sampnames = []
    print("Running ElitePileups")
    for x in args.bamfile:
        samp = elitePileups(x, args.genes, args.assembly)
        sampnames.append(samp)
        print(f'Finished pileup of {samp}')
    
    print("Running CallEliteVariants")
    freq_var, freq_df = callEliteVariants(args.desman, sampnames, args.assembly)
    
    print("Running estimateStrainCountDesman")
    fits = estimateStrainCountDesman(args.desman, freq_var, freq_df)
    
    print("Plotting strain replicate data")
    plotDev(args.desman, fits, args.output)
        
def findEliteGenes(desman : str, contigs : str, assembly : str) -> Tuple[str, str]:
    cmd = ["python3", desman + "/scripts/extract_species_contigs.py", 
           assembly, contigs]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, stdout=open("species_contigs.fa", "w"))
    
    # Clear PS_temp directory if present so that files can be overwritten
    if os.path.isdir('PS_temp'):
        sp.run('rm -r PS_temp', shell=True)
    
    cmd = [desman + "/external/phylosift_v1.0.1/phylosift", "search", "--besthit",
           "--isolate", "species_contigs.fa"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(' '.join(cmd), shell=True)
    
    cmd = [desman + "/external/phylosift_v1.0.1/phylosift", "align", "--besthit",
           "--isolate", "species_contigs.fa"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(' '.join(cmd), shell=True)
    
    cmd = ["python3", desman + "/scripts/get_elite_range.py", 
           "PS_temp/species_contigs.fa/blastDir/lookup_ID.1.tbl", 
           "PS_temp/species_contigs.fa/alignDir/DNGNGWU*.codon.updated.1.fasta"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(' '.join(cmd), shell=True, check = True, stdout=open("elites.bed", "w"))
    
    return ["elites.bed", "species_contigs.fa"]

def getFormatBedLine(bedfh : TextIO) -> str:
    for l in bedfh:
        segs = l.split()
        yield segs[0] + ":" + segs[1] + "-" + segs[2]

def elitePileups(bam : str, elites : str, assembly : str) -> str:
    # Get the sample name from the basename of the bed file
    bfile = os.path.basename(bam)
    bsegs = re.split("\.", bfile)
    change = re.compile('[:-]')
    files = []
    with open(elites, 'r') as bedfh:
        for region in getFormatBedLine(bedfh):
            safe = re.sub(change, '_', region)
            cmd = ["samtools", "mpileup", "-r", 
                   region, "-f", assembly, bam ]
            print(f'Cmd: {" ".join(cmd)}')
            sp.run(cmd, check = True, stdout=open(bsegs[0] + "." + safe + ".pileup", "w"))
            files.append(bsegs[0] + "." + safe + ".pileup")
    
    with open(bsegs[0] + ".pileup", 'w') as out:
        for f in files:
            with open(f, 'r') as fh:
                for l in fh:
                    l = l.rstrip()
                    out.write(l + '\n')
            os.remove(f)
    return bsegs[0] + ".pileup"

def callEliteVariants(desman : str, pileup : List[str], assembly : str) -> Tuple[str, str]:
    cmd = ["python3", desman + "/scripts/pileups_to_freq_table.py", 
           assembly]#, pileup, "desmanfreqs.csv"]
    cmd.extend(pileup)
    cmd.append("desmanfreqs.csv")
                
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=False)
    
    cmd = ["python3", desman + "/desman/Variant_Filter.py", 
           "desmanfreqs.csv", "-o", "dfreqs", "-p", "-m", "1.0", "-f", "25.0", "-c",
           "-sf", "0.50", "-t", "2.5"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=False)
    
    return ["dfreqssel_var.csv", "dfreqstran_df.csv"]

def executeDesman(desman : str, freq_var : str, freq_df : str, out : str, g : str, s : str) -> str:
    cmd = ["python3", desman + "/bin/desman", freq_var, "-e", freq_df, "-o", 
           f'cluster_{g}_{s}', "-r", "1000", "-i", "100", "-g", g, "-s", 
           s]
    print(f'Running strain count round {g} and replicate {s}')
    sp.run(cmd, stdout=open(out + ".out", "w"))
    
    sp.run(f'cp {out}/fit.txt fit_{g}_{s}.txt', shell=True)
    return f'fit_{g}_{s}.txt'

def estimateStrainCountDesman(desman : str, freq_var : str, freq_df : str) -> str:
    fits = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        for g in range(1,8):
            for repid in range(0,10):
                print(f'Queueing: exDesman {desman} {freq_var} {freq_df} cluster_{g}_{repid}') 
                try:
                    temp = executor.submit(executeDesman, desman, freq_var, freq_df, 
                                       f'cluster_{g}_{repid}', str(g), str(repid))
                except NameError as nm:
                    print(nm.args)
                fits.append(temp)
                """
                # Saving this in case my method doesn't work
                cmd = [desman + "/bin/desman", freq_var, "-e", freq_df, "-o",
                       f'cluster_{g}_{repid}', "-r", "1000", "-i", "100", "-g", str(g), "-s", 
                       str(repid)]
                print(f'Running strain count round {g} and replicate {repid}')
                sp.run(cmd, stdout=open(f'cluster_{g}_{repid}.out', "w"))
            
                sp.run(f'cp cluster_{g}_{repid}/fit.txt fit_{g}_{repid}.txt', shell=True)
                fits.append(f'fit_{g}_{repid}.txt')
                """
            
    with open("desman_dic.fits", 'w') as out:
        out.write("H,G,LP,Dev\n")
        for x in fits:
            temp = ""
            try:
                temp = x.result()
            except Exception as ex:
                print(f'{x} process created an exception!')
                print(x.args)
            with open(temp, 'r') as fh:
                for lines in fh:
                    lines = lines.rstrip()
                    segs = re.split(',', lines)
                    segs.pop(0)
                    out.write(','.join(segs) + "\n")
                    
    return "desman_dic.fits"

def plotDev(desman : str, alldics : str, output : str) -> None:
    cmd = [desman + "/scripts/PlotDev.R", "-l", alldics, "-o", output + ".pdf"]
    print(f'Cmd: {" ".join(cmd)}')
    sp.run(cmd, shell=False)
    
    print(f'Completed strain testing! Plots are in Dev.pdf')

if __name__ == "__main__":
    args = parse_user_input()
    main(args)
