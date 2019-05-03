# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:19:15 2019

@author: dbickhart
"""

import argparse;

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Process long-read alignments and Hi-C links to identify likely viral-host associations in metagenomic assembly data"
            )
    parser.add_argument('-l', '--long_read', 
                        help="Input PAF file for long-read alignments to other contigs",
                        type=str, required=True
                        )
    parser.add_argument('-h', '--hic_links',
                        help="Input Hi-C bipartite link graph",
                        type=str, required=True
                        )
    parser.add_argument('-v', '--viruses',
                        help="Tab delimited list of contigs containing viral sequence and their lengths",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--link_thresh',
                        help="Filter for the number of stdevs of Hi-C links above the average count to be used in viral Hi-C association [2.5]",
                        type=float, default=2.5
                        )
    parser.add_argument('-v', '--overhang',
                        help="Filter for long-read overhang off of viral contigs [150]",
                        type=int, default=150
                        )
    parser.add_argument('-o', '--output',
                        help="Output basename",
                        type=str, required=True
                        )
    
    return parser.parse_args()

def main(args):
    print("Do your stuff, hipster language!")
    
class viralComparison:

    def __init__(self, viralCtgFile : str):
        self.viruses = dict()
        
        # Load contigs and lengths
        with open(viralCtgFile, 'r') as fh:
            for l in fh:
                l = l.rstrip()
                s = l.split()
                self.viruses[s[0]] = s[1]
        
        self.ovlpSizes = list()
        self.readErrors = 0
        
        print(f'Loaded {len(self.viruses)} viral contigs for analysis\n')
        
    def isVirus(self, ctg : str):
        return ctg in self.viruses
    
    def getOverhang(self, rname : str, vctg : str, othresh : int, vstart : int, 
                    vend : int, rstart : int, rend : int, rlen: int, rorient : str,
                    ofh):
        right_end = rlen - rend < othresh
        read_unmapright = self.viruses[vctg] - vend > self.viruses[vctg] - vstart
        
        if rorient == "+":
            # read: ------
            # ctg:  --
            if right_end and read_unmapright:
                ofh.write(f'{rname}\t{rend}\t{rlen}\t{vctg}\n')
                self.ovlpSizes.append(rlen - rend)
            elif not right_end and not read_unmapright:
                # read: ------
                # ctg:      --
                ofh.write(f'{rname}\t0\t{rstart}\t{vctg}\n')
    
if __name__ == "__main__":
    args = parse_user_input()
    main(args)