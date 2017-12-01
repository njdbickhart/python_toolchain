#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Variant calling helper utilities

"""

class haplotype:
    
    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        
class individual:
    
    def __init__(self, regid):
        self.regid = regid
        self.haplotypes = []
        self.haplookup = set()
        self.vcfs = []
        
    def loadHap(self, hap):
        self.haplotypes.append(hap)
        self.haplookup.add(hap)
        
    def hasHap(self, hap):
        return hap in self.haplookup