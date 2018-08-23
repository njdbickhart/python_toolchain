# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:56:15 2018

@author: dbickhart
"""

import argparse
import os

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
                        help="Estimated count of strains",
                        required=True, type=int
                        )
    
    return parser.parse_args()

def main(args):
    # main routine
    print("hey")
    os.chdir(args.output)
    if not os.path.isfile('dfreqstran_df.csv') and not os.path.isfile('dfreqssel_var.csv'):
        print("Missing essential files dfreqstran_df.csv and dfreqssel_var.csv!")
        os.sys.exit(-1)
        
    if not os.path.isfile('species_contigs.fa') and not os.path.isfile('elites.bed'):
        print("Missing essential files elites.bed and species_contigs.fa!")
        os.sys.exit(-1)
        
    

    

if __name__ == "__main__":
    args = parse_user_input()
    main(args)