# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:43:36 2023

@author: Derek.Bickhart
"""

import re
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
from collections import defaultdict

revdict = {'T' : 'A', 'G' : 'C', 'C' : 'G'}

class Marker:
    
    def __init__(self, snpID, Ref, Alt, Chr, Pos, population, reverse = False):
        # Defining attributes for later filling
        self.snpID = snpID
        self.Ref = Ref
        self.Alt = Alt
        self.Chr = Chr
        self.Pos = Pos
        self.population = population
        self.reverse = reverse        
                
        if reverse:
            self.Ref = revdict.get(self.Ref, self.Ref)
            self.Alt = revdict.get(self.Alt, self.Alt)
        
class Genotype:
    
    def __init__(self, marker, allele, reverse = False):
        # Defining attributes for later filling
        self.marker = marker
        self.variantState = self._determineVariantState(allele, reverse)
        
    def  _determineVariantState(self, allele, reverse=False):

        if allele == "0" or allele == "1" or allele == "2" or allele == "5":
            if allele == '5':
                return '.'
            else:
                return allele
        elif re.search(r'^.\/.', allele):
            asegs = re.split(r'[:\/]', allele)
            var = 0
            for x in range(2):
                if asegs[x] == '0':
                    var += 0 
                elif asegs[x] == '1':
                    var += 1
                else:
                    var = '' # Multivariant sites in VCF files will not be handled right now
            return str(var)
        else:
            asegs = allele.split()
            var = 0
            for a in asegs:
                if a == self.marker.Ref:
                    var += 1 if reverse else 0
                elif a == self.marker.Alt:
                    var += 0 if reverse else 1
            
            #var = "".join(asegs)
            #if reverse:
            #    var = "".join(revdict.get(base, base) for base in var)
            return var
        
    def isConcordant(self, genotype):
        return self.variantState == genotype.variantState
    
class GenotypeFactory:
    
    def __init__(self, markerInfoFile, skiplist = []):
        self.markerInfoFile = markerInfoFile
        self.MarkerList = list()
        self.MarkerByID = dict()
        # Data structure to identify genotypes from VCF files
        self.MarkerByCoord = dict()
        
        self.skip = set(skiplist)
        self.dsNameList = list()
        self.individualSetA = set()
        self.individualSetB = set()
        # Keys: dataset -> animal -> SNPID
        self.GenotypeKeys = defaultdict(lambda: defaultdict(dict))
        
    # Assumes that markers are in a "bim" file format with the forward strand alleles in the last two columns
    def loadMarkers(self, population='test'):
        loaded = 0
        with open(self.markerInfoFile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                loaded += 1
                rev = False
                if len(s) == 7:
                    rev = True if s[6] == '1' else False
                m = Marker(s[1], s[4], s[5], s[0], s[3], population, reverse=rev)
                self.MarkerList.append(m)
                self.MarkerByID[s[1]] = m
                self.MarkerByCoord[f'{s[0]}:{s[3]}'] = m
                
        print(f'Generated {loaded} marker site classes')
        
    # Convert a BED file to a PED file:
    # plink --bfile template/LW_HyporGeno2 --recode --tab --out LW_test --extract test_markers.list
    # PED files are tab delimited but have space delmited genotypes
    def loadPED(self, pedfile, dsName, runReverse = False):
        self.dsNameList.append(dsName)
        with open(pedfile, 'r') as input:
            icount = 0
            msum = 0
            firstelement = True
            fixmarker = False
            for l in input:
                s = l.rstrip().split("\t")
                if s[0] in self.skip:
                    continue
                self.individualSetA.add(s[0])
                
                allele = ''
                icount += 1
                msum += len(s) - 6
                for i in range(6, len(s)):
                    if firstelement:
                        if re.search(r'\d+ \d+', s[i]):
                            print(s[i])
                            fixmarker = True
                        firstelement = False
                    if fixmarker:
                        tsegs = s[i].split()
                        tallele = 0
                        for j in tsegs:
                            if int(j) == 2:
                                tallele += 1
                        allele = str(tallele)
                    else:
                        allele = str(s[i])
                    marker = self.MarkerList[i - 6]
                    rev = True if runReverse and marker.reverse else False
                    g = Genotype(marker, allele, reverse=rev)
                    self.GenotypeKeys[dsName][s[0]][self.MarkerList[i - 6].snpID] = g
            print(f'PED: Out of {icount} individuals, had an average of {msum / icount} markers per individual')
            
    # Load VCF entries assuming non phased encodings
    def loadVCF(self, vcffile, dsName, runReverse = False):
        self.dsNameList.append(dsName)
        icount = 0
        msum = 0
        skipCols = set()
        with open(vcffile, 'r') as input:
            # Animal IDs start at the 10th (9 - 0-base) column
            header = list()
            for l in input:
                if l.startswith('##'):
                    continue
                elif l.startswith('#'):
                    header = l.rstrip().split()
                    icount = len(header) - 10
                    for i, anim in enumerate(header[9:]):
                        if anim in self.skip:
                            skipCols.add(i + 9)
                        else:
                            self.individualSetB.add(anim)
                    continue
                s = l.rstrip().split()
                coord = f'{s[0]}:{s[1]}'
                if coord in self.MarkerByCoord:
                    msum += 1
                    m = self.MarkerByCoord[coord]
                    rev = True if runReverse and m.reverse else False
                    for i in range(9, len(header)):
                        if i in skipCols:
                            continue
                        anim = header[i]
                        allele = s[i]
                        g = Genotype(m, allele, reverse=rev)
                        self.GenotypeKeys[dsName][anim][m.snpID] = g
                elif s[2] in self.MarkerByID:
                    msum += 1
                    m = self.MarkerByID[s[2]]
                    rev = True if runReverse and m.reverse else False
                    for i in range(9, len(header)):
                        if i in skipCols:
                            continue
                        anim = header[i]
                        allele = s[i]
                        g = Genotype(m, allele, reverse=rev)
                        self.GenotypeKeys[dsName][anim][m.snpID] = g
        print(f'VCF: Out of {icount} individuals, had an average of {msum / icount} markers per individual')
     
    # Returns: Contingency Table, By Marker DF, by Animal DF
    def compareTwoDatasets(self):
        # Contingency table starter
        contingency = defaultdict(list)
        # marker, Identical call count, total calls, Dataset1 p freq, Dataset1 q freq, Dataset2 p freq, Dataset2 q freq
        bymarker = defaultdict(list)
        # Animal, Dataset1 homozygotes proportion, Dataset2 homozygotes proportion, 
        byanimal = defaultdict(list)
        
        # Prepare the contingency table and homozygote proportion                    
        for anim in self.individualSetB:
            aHomCount = defaultdict(int)
            totCount = defaultdict(int)
            for m in self.MarkerList:                
                contingency['Animal'].append(anim)
                contingency['SNPID'].append(m.snpID)
                for ds in self.dsNameList:
                    gtype = self.GenotypeKeys[ds][anim][m.snpID].variantState
                    contingency[ds].append(gtype)
                    if gtype == "0" or gtype == "2":
                        aHomCount[ds] += 1
                    else:
                        aHomCount[ds] += 0
                    totCount[ds] += 1
            for d, v in aHomCount.items():
                byanimal[f'{d}_hom'].append(v / totCount[d])
        
        byanimal['Animal'].extend([x for x in self.individualSetB])
        
        # Create contingency table
        dss = list(self.dsNameList)
        cDF = pd.DataFrame(contingency)
        ctDF = pd.crosstab(cDF[dss[0]], cDF[dss[1]], normalize=True)
        
        # Create the byanimal dataframe
        byADF = pd.DataFrame(byanimal)
        
        # Create marker level stats
        
        #byM = cDF[['SNPID'] + dss].groupby(['SNPID'], as_index=False).agg({'identical' : lambda x: })
        
        return ctDF, cDF, byADF