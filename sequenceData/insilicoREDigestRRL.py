# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 10:37:38 2019

@author: dbickhart
"""

import argparse
from Bio import SeqIO, SeqRecord, RestrictionBatch
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

def main(arg):
    print("hey")
    
def fragmentSeq(seq : SeqRecord, rb : RestrictionBatch, ren : str ) -> str:
    sid = seq.id
    slen = len(seq.seq)
    
    print(f'Running insilico digest on {sid} of length {slen}')
    # Running an analysis on this sequence
    analysis = res.Analysis(rb, seq.seq)
    
    # locations of cut sites for this particular restriction enzyme
    resites = analysis.full()[ren]
    
    