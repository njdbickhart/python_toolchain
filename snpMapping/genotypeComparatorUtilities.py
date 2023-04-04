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

class Marker:
    
    def __init__(self, snpID, Ref, Alt, Chr, Pos, population):
        # Defining attributes for later filling
        self.snpID = snpID
        self.Ref = Ref
        self.Alt = Alt
        self.Chr = Chr
        self.Pos = Pos
        self.population = population
        
class Genotype:
    
    def __init__(self, marker, allele):
        # Defining attributes for later filling
        self.marker = marker
        self.variantState = self._determineVariantState(allele)
        
    def  _determineVariantState(self, allele):
        if allele == "0" or allele == "1" or allele == "2" or allele == "5":
            if allele == '5':
                return '.'
            else:
                return allele
        elif re.search(r'^\d{1}\/\d{1}\:', allele):
            asegs = re.split(r'[:\/]', allele)
            var = 0
            for x in range(2):
                if asegs[x] == '0':
                    var += 0 
                elif asegs[x] == '1':
                    var += 1
                else:
                    var = '.' # Multivariant sites in VCF files will not be handled right now
            return str(var)
        
    def isConcordant(self, genotype):
        return self.variantState == genotype.variantState
    
class GenotypeFactory:
    
    def __init__(self, markerInfoFile):
        self.markerInfoFile = markerInfoFile
        self.MarkerList = list()
        self.MarkerByID = dict()
        # Data structure to identify genotypes from VCF files
        self.MarkerByCoord = dict()
        
        self.dsNameList = list()
        self.individualSet = set()
        # Keys: dataset -> animal -> SNPID
        self.GenotypeKeys = defaultdict(lambda: defaultdict(dict))
        
    # Assumes that markers are in a "bim" file format with the forward strand alleles in the last two columns
    def loadMarkers(self, population='test'):
        loaded = 0
        with open(self.markerInfoFile, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                loaded += 1
                m = Marker(s[1], s[4], s[5], s[0], s[3], population)
                self.MarkerList.append(m)
                self.MarkerByID[s[1]] = m
                self.MarkerByCoord[f'{s[0]}:{s[3]}'] = m
                
        print(f'Generated {loaded} marker site classes')
        
    # Convert a BED file to a PED file:
    # plink --bfile template/LW_HyporGeno2 --recode --tab --out LW_test --extract test_markers.list
    # PED files are tab delimited but have space delmited genotypes
    def loadPED(self, pedfile, dsName):
        self.dsNameList.append(dsName)
        with open(pedfile, 'r') as input:
            icount = 0
            msum = 0
            for l in input:
                s = l.rstrip().split("\t")
                self.individualSet.add(s[0])
                firstelement = True
                fixmarker = False
                allele = ''
                icount += 1
                msum += len(s) - 6
                for i in range(6, len(s)):
                    if firstelement:
                        if re.search(r'\d \d', s[i]):
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
                    g = Genotype(self.MarkerList[i - 6], allele)
                    self.GenotypeKeys[dsName][s[0]][self.MarkerList[i - 6].snpID] = g
            print(f'PED: Out of {icount} individuals, had an average of {msum / icount} markers per individual')
            
    # Load VCF entries assuming non phased encodings
    def loadVCF(self, vcffile, dsName):
        self.dsNameList.append(dsName)
        icount = 0
        msum = 0
        with open(vcffile, 'r') as input:
            # Animal IDs start at the 10th (9 - 0-base) column
            header = list()
            for l in input:
                if l.startswith('##'):
                    continue
                elif l.startswith('#'):
                    header = l.rstrip().split()
                    icount = len(header) - 10
                    for anim in header[9:]:
                        self.individualSet.add(anim)
                    continue
                s = l.rstrip().split()
                coord = f'{s[0]}:{s[1]}'
                if coord in self.MarkerByCoord:
                    msum += 1
                    m = self.MarkerByCoord[coord]
                    for i in range(9, len(header)):
                        anim = header[i]
                        allele = s[i]
                        g = Genotype(m, allele)
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
        for anim in self.individualSet:
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
                    totCount[ds] += 1
            for d, v in aHomCount.items():
                byanimal[f'{d}_hom'].append(v / totCount[d])
        
        byanimal['Animal'].extend([x for x in self.individualSet])
        
        # Create contingency table
        dss = list(self.dsNameList)
        cDF = pd.DataFrame(contingency)
        ctDF = pd.crosstab(cDF[dss[0]], cDF[dss[1]], normalize=True)
        
        # Create the byanimal dataframe
        byADF = pd.DataFrame(byanimal)
        
        # Create marker level stats
        
        #byM = cDF[['SNPID'] + dss].groupby(['SNPID'], as_index=False).agg({'identical' : lambda x: })
        
        return ctDF, cDF, byADF