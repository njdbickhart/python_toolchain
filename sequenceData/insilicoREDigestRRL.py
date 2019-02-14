# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 10:37:38 2019

@author: dbickhart
"""

import argparse
from Bio import Restriction as res

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Generate restriction enzyme maps for an assembly fasta and select fragments of a specific size range"
            )
    parser.add_argument('-f', '--file', 
                        help="The input sequence fasta file",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--cutoff',
                        help="The mean fragment basepair size to retain after digest",
                        type=int, required=True
                        )
    parser.add_argument('-s', '--stdev',
                        help="The standard deviation (rounded up) around the mean for kept fragments",
                        type=int, required=True
                        )
    parser.add_argument('-e', '--enzyme',
                        help="The name of the restriction enzyme to use in the analysis",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="Text output file name",
                        type=str, required=True
                        )
    
    return parser.parse_args()

