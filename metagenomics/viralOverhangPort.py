# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:19:15 2019

@author: dbickhart
"""

import argparse
import subprocess
from collections import defaultdict
import numpy as np
import scipy.cluster.vq as sp

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Process long-read alignments and Hi-C links to identify likely viral-host associations in metagenomic assembly data"
            )
    parser.add_argument('-l', '--long_read', 
                        help="Input PAF file for long-read alignments to other contigs",
                        type=str, required=True
                        )
    parser.add_argument('-b', '--blob_tools',
                        help="Input blob tools data table",
                        type=str, required=True
                        )
    parser.add_argument('-s', '--hic_links',
                        help="Input sam/bam file with alignments of Hi-C reads",
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
    parser.add_argument('-a', '--overhang',
                        help="Filter for long-read overhang off of viral contigs [150]",
                        type=int, default=150
                        )
    parser.add_argument('-n', '--noplot',
                        help="[optional flag] Disables plotting of data",
                        action='store_true'
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
        self.ovlpEC = defaultdict(list) # {readname} -> [start, end, vctg]
        self.ovlpASM = defaultdict(dict)
        self.hicASM = defaultdict(dict)        
        
        self.totalCount = 0
        self.finalTable = defaultdict(dict)
        
        self.readErrors = 0
        
        print(f'Loaded {len(self.viruses)} viral contigs for analysis\n')
        
    def isVirus(self, ctg : str):
        return ctg in self.viruses
    
    def printOutFinalTable(self, outtab : str):
        print("Final associations in table {}".format(outtab))
        #TODO: finish the printout
    
    def loadTaxonomy(self, blobtable : str):
        print(f'Loading taxonomic information from file {blobtable}')
        # first, get the set of contig names that we need to process
        contigs = set()
        for v, j in self.finalTable.items():
            contigs.add(v)
            for h, c in self.finalTable.items():
                contigs.add(h)
              
        taxonomy = defaultdict(list) #{contig} -> [kingdom, genus]
        with open(blobtable, 'r') as input:
            # read until first line and get index information
            kingidx = 0
            genusidx = 0
            while(True):
                line = input.readline()
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    line = line.rstrip()
                    segs = line.split()
                    for i in range(len(segs)):
                        if segs[i].startswith("genus"):
                            genusidx = i
                        if segs[i].startswith("superkingdom"):
                            kingidx = i
                    break
            if kingidx == 0 or genusidx == 0:
                print("Could not ID taxonomic table entries! Did you remember to add genus and superkingdom table headers?")
                return
            
            for l in input.readlines():
                l = l.rstrip()
                segs = l.split()
                if segs[0] in contigs:
                    taxonomy[segs[0]] = [segs[kingidx], segs[genusidx]]
        
        # Now assign tax affiliations to the contigs
        taxassign = 0
        for v, j in self.finalTable.items():
            for h, c in self.finalTable.items():
                if v in taxonomy and h in taxonomy:
                    taxassign += 1
                    c.setTaxonomy(taxonomy[v][1], taxonomy[h][0], taxonomy[h][1])
        print(f'Filled in taxonomic entries for {taxassign} associations out of {self.totalCount} possible entries')
                
    def combineTables(self, reads : bool, hic : bool):
        print("Generating final network table")
        if reads:
            for v, j in self.ovlpASM.items():
                for h, c in j.items():
                    self.totalCount += 1
                    self.finalTable[v][h] = c
           
        overlaps = 0
        if hic:
            for v, j in self.hicASM.items():
                for h, c in j.items():
                    if v in self.finalTable:
                        if h in self.finalTable[v]:
                            overlaps += 1
                            self.finalTable[v][h].combine(c)
                            continue
                    self.totalCount += 1
                    self.finalTable[v][h] = c
        
        if reads and hic:
            print("Identified {} overlapping viral-host associations.".format(overlaps))
    
    def generateHiCLinkTable(self, samtools : str, insam : str, outtab : str):
        print("Generating a Hi-C link table from alignment file: {}".format(insam))
        hicLinks = defaultdict(lambda : defaultdict(int))
        selfLinks = defaultdict(float)
        cmd = [samtools, 'view', insam]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
            for l in proc:
                if not l.startswith('@'):
                    segs = l.split()
                    if segs[6] == segs[2]:
                        selfLinks[segs[2]] += 0.5   # half a count to avoid double-counting pairs
                    if segs[6] != segs[2] and self.isVirus(segs[2]):
                        hicLinks[segs[2]][segs[6]] += 1
        
        # Let's try k-means clustering to divide the samples
        # assuming that viral intercontig link noise will be in cluster 1 and signal in cluster 2
        vLCounts = list()
        for v, j in hicLinks.items():
            for h, n in j.items():
                vLCounts.append(float(n))
               
        flattened = array(vLCounts)
        centroids,_ = sp.kmeans(flattened, 2)
        idx,_ = sp.vq(flattened, centroids)
        
        minV1 = min([int(x) for x in flattened[idx==0]])
        minV2 = min([int(x) for x in flattened[idx==1]])
        
        kThresh = minV1 if minV1 > minV2 else minV2
        SNR = np.where(flattened.std(axis=0, ddof=0) == 0, 0, 
                       flattened.mean(0) / flattened.std(axis=0, ddof=0))
        
        print(f'K-means clustering identified a minimum intercontig link count of {kThresh} from {len(vLCounts)} observations with a signal-to-noise ratio of {SNR[0]}')
        
        # Now generate intermediary output tab file and fill data structure
        vassocNum = list() # list of number of host contigs associated with viruses
        with open(outtab, 'w') as out:
            for v, j in hicLinks.items():
                vassocNum.append(len(j))
                for h, n in j.items():
                    if n >= kThresh:
                        out.write(f'{v}\t{h}\t{n}\n')
                        self.hicASM[v][h] = VAssoc(h, v, "HiC", n)
                    
        meanHcontig = np.mean(vassocNum)
        print(f'Found valid Hi-C link associations for {len(hicLinks)} viral contigs out of {len(flattened)} original candidates.')
        print(f'There were an average of {meanHcontig} host contig associations in this dataset')
    
    def samfaidx(self, samtools : str, ctglst, outfasta):
        cmd = [samtools, 'faidx']
        cmd.append(ctglst)
        with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
            for l in proc:
                outfasta.write(l)
    
    def realignECOverhangs(self, asmCtgFasta : str, ecReads : str, minimap : str,
                           minimapOpts =['-x', 'map-pb'], outfasta: str, outfile : str,
                           samtools : str, rcountThresh = 3):
        print("Generating overhangs in {} fasta file".format(outfasta));
        # Create the overhanging read fasta for alignment
        with open(outfasta, 'w') as fasta:
            container = list()
            for k, f in self.ovlpEC.items():
                container.append(f'{k}:{f[0]}-{f[1]}')
                if len(container) >= 100:
                    self.samfaidx(samtools, container, fasta)
                    container = list()
            
            if len(container) > 0:
                self.samfaidx(samtools, container, fasta)
                container = list()
        
        # Now, map the reads and filter the results
        print("Aligning overhangs to full genome fasta file")
        redundancies = 0
        overlaps = dict() # temp container for EC read associations
        with subprocess.Popen([minimap, minimapOpts, asmCtgFasta, outfasta], stdout=subprocess.PIPE) as proc, open(outfile, 'w') as out:
            for l in proc:
                l = l.rstrip()
                segs = l.split()
                
                # Get original alignments
                vir = self.ovlpEC[segs[0]]
                if segs[0] not in overlaps and not self.isVirus(segs[5]):
                    # Note: this preferentially prints out the first alignment encountered
                    overlaps[segs[0]] = VAssoc(segs[5], vir[2])
                    out.write(f'{segs[0]}\t{segs[5]}\t{vir[2]}')
                else:
                    overlaps[segs[0]].setRedundant()
                    # Ignore subsequent alignments
                    redundancies += 1
        
        successes = len({s.vctg for k, s in overlaps.items()})                    
        print(f'Found successful associations for {successes} out of {len(self.viruses)} viral contigs and {redundancies} ambiguously aligned reads')
        
        # Loading into the final container
        readcounts = defaultdict(lambda : defaultdict(int))
        for k, s in overlaps.items():
            # Skip over any reads that may have been ambiguous alignments 
            if not s.redund:
                readcounts[s.vctg][s.hostctg] += 1
                
        for v, j in readcounts.items():
            for h, c in j.items():
                if c >= rcountThresh:
                    # The above selects only associations that have at least X read alignments
                    self.ovlpASM[v][h] = VAssoc(h, v, "Read", c)
    
    def alignECReads(self, viralCtgFasta : str, ecReads : str, minimap : str, 
                     minimapOpts = ['-x', 'map-pb'], algLen = 500, oThresh = 200,
                     outfile : str):
        cmd = [minimap]
        cmd.extend(minimapOpts)
        cmd.extend([viralCtgFasta, ecReads])
        
        with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc, open(outfile, 'w') as out:
            for l in proc:
                l = l.rstrip()
                segs = l.split()
                
                if int(segs[9]) > algLen and (int(segs[7]) < oThresh or int(segs[6]) - int(segs[8]) < oThresh):
                    self.getOverhang(segs[0], segs[5], oThresh, int(segs[7]), 
                                     int(segs[8]), int(segs[2]), int(segs[3]), 
                                     int(segs[1]), segs[4], out)
        count = len(self.ovlpSizes)
        mean = np.mean(self.ovlpSizes)
        viruscounts = len({s[2] for s in self.ovlpEC})
        print(f'Identified {count} overlapping reads with a mean length of {mean} for {viruscounts} unique viral contigs')
    
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
                self.ovlpEC[rname] = [rend, rlen, vctg]
            elif not right_end and not read_unmapright:
                # read: ------
                # ctg:      --
                ofh.write(f'{rname}\t0\t{rstart}\t{vctg}\n')
                self.ovlpSizes.append(rstart)
                self.ovlpEC[rname] = [1, rstart, vctg]
            else:
                self.readErrors += 1
        else:
            # read: -----
            # ctg:  --
            if right_end and not read_unmapright:
                ofh.write(f'{rname}\t0\t{rstart}\t{vctg}\n')
                self.ovlpSizes.append(rstart)
                self.ovlpEC[rname] = [1, rstart, vctg]
            elif not right_end and not read_unmapright:
                # read: ------
                # ctg:     ---
                ofh.write(f'{rname}\t{rend}\t{rlen}\t{vctg}\n')
                self.ovlpSizes.append(rlen - rend)
                self.ovlpEC[rname] = [rend, rlen, vctg]
            else:
                self.readErrors += 1

class VAssoc:
    def __init__(self, hostctg : str, vctg : str, cat = "Read", count = 0):
        self.hostctg = hostctg
        self.vctg = vctg
        self.redund = False
        self.category = cat
        self.count = 0
        # Taxonomic placeholders
        self.vTax = "N/A"
        self.hostKing = "N/A"
        self.hostGenus = "N/A"
        
        # used only if combined
        self.complex = defaultdict(int)
    
    def setRedundant(self):
        self.redund = True      
        
    def setTaxonomy(self, vTax : str, hKing : str, hGen : str):
        self.vTax = vTax
        self.hostKing = hKing
        self.hostGenus = hGen
    
    def combine(self, other : VAssoc):
        self.complex[self.category] = self.count
        self.complex[other.category] = other.count
        self.category = "Both"
        
    def getEvidence(self) -> str:
        if self.category == "Both":
            return ';'.join([f'{k}:{v}' for k, v in self.complex.items()])
        else:
            return f'{self.category}:{self.count}'
        
    def getListOutput(self) -> list:
        return [self.vctg, self.hostctg, self.vTax, self.hostKing, 
                self.hostGenus, self.getEvidence()]
  
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
