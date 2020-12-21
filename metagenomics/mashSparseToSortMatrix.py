# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 14:05:37 2020

@author: derek.bickhart-adm
"""

import argparse
import numpy as np
import os

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Take pairwise Mash distance scores and order into a sorted matrix"
            )
    parser.add_argument('-r', '--rlist', 
                        help="Input list of bins for the reference set",
                        required=True, type=str
                        )
    parser.add_argument('-q', '--qlist',
                        help="Input list of bins for the query set",
                        required=True, type=str
                        )
    parser.add_argument('-d', '--distance',
                        help="Input tab delimited list of distance scores for bins",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are (basename).output and (basename).output.scored",
                        required=True, type=str,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    worker = sparseToFullMatrix()
    
    worker.createIdx(args.rlist, args.qlist)
    
    with open(args.distance, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            s[0] = os.path.basename(s[0])
            s[1] = os.path.basename(s[1])
            worker.updateMat(s[0], s[1], 1.0 - float(s[2]))
            
    worker.sortDiagonal()
    
    worker.printOutMat(args.output)
    
    
class sparseToFullMatrix:
    
    def __init__(self):
        self.rlistIdx = dict()
        self.qlistIdx = dict()
        
        self.idtoRefID = dict()
        self.idtoQID = dict()
        
        self.rlistnum = 0
        self.qlistnum = 0
        
        self.mat = None
        
        
    def updateMat(self, refid, qid, value):
        self.mat[self.rlistIdx[refid]][self.qlistIdx[qid]] = value
        
    def createIdx(self, reflist, qlist):
        with open(reflist, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                self.rlistIdx[s[0]] = self.rlistnum
                self.rlistnum += 1
                
        with open(qlist, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                self.qlistIdx[s[0]] = self.qlistnum
                self.qlistnum += 1
                
        self.idtoRefID = {v : k for k, v in self.rlistIdx.items()}
        self.idtoQID = {v : k for k, v in self.qlistIdx.items()}
        self.mat = self.getEmptyMat()
                
        
    def getEmptyMat(self):
        return [[0 for col in range(self.qlistnum)] for row in range(self.rlistnum)]
    

    def sortDiagonal(self):
        
        self.mat = np.array(self.mat)
        
        indicies = np.argsort(np.diag(self.mat))
        self.mat = self.mat[indicies,:]
        
        nid = 0
        for i in indicies:
            self.rlistIdx[self.idtoRefID[i]] = nid
            nid += 1
            
        self.idtoRefID = {v : k for k, v in self.rlistIdx.items()}
        
            
    def printOutMat(self, outfile):
        with open(outfile, 'w') as out:
            out.write("\t" + "\t".join([j for j in self.idtoQID.keys()]) + "\n")
            for i in range(self.rlistnum):
                out.write(self.idtoRefID + "\t" + "\t".join(self.mat[i]) + "\n")
        

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
