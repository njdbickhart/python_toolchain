

import argparse
from typing import Dict, Tuple, List
from collections import defaultdict

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "A script to pull flanking SNP sequence to "
            )
    parser.add_argument('-v', '--vcf', 
                        help="Input vcf file with variant calls",
                        required=True, type=str)
    parser.add_argument('-l', '--list', 
                        help="Input list of animals and haplotype associations",
                        required=True, type=str)
    parser.add_argument('-o', '--output', 
                        help="Base output directory",
                        required=True, type=str)
    parser.add_argument('-s', '--segments',
                        help="Segment coordinates for haplotypes",
                        required=True, type=str)
   
    #parser.add_argument('-m', '--meta',
                        #help="[Optional] Add multiple targets for metagenome assembly",
                        #action='store_true') 
    return parser.parse_args()
