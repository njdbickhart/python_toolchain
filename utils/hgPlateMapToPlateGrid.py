# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:04:41 2025

@author: Derek.Bickhart
"""

import argparse
import string

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script to convert platemap files into plate grids"
            )
    parser.add_argument('-f', '--file', 
                        help="Input platemap file with 1. platename, 4. row, 5. column, 6. animalid",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name. Tab delimited",
                        required=True, type=str,
                        )
    parser.add_argument('-r', '--row',
                        help="Row column number 0-based [3]",
                        default=3, type=int,
                        )
    parser.add_argument('-c', '--col',
                        help="Plate column number 0-based [4]",
                        default=4, type=int,
                        )
    parser.add_argument('-a', '--an',
                        help="Animal ID column number 0-based [5]",
                        default=5, type=int,
                        )
    parser.add_argument('-e', '--empty',
                        help="Empty well string",
                        default="NTC", type=str,
                        )
    return parser.parse_args(), parser

# This variable helps to sort a linear numerical sample number (1-96) into a plate well format (A1-H8)
hexCon = {'A' : 0, 'B' : 12, 'C' : 24, 'D' : 36, 'E' : 48, 'F' : 60, 'G' : 72, 'H' : 84}

def transform_well(well):
    hexAdd = hexCon.get(well[0], -96)
    return (hexAdd + int(well[1:]))

class plate:
    
    def __init__(self, name):
        self.name = name
        
        # Integer number -> sampleID
        self.samples = dict() 
        
    def loadSample(self, well, sampleID):
        if well in self.samples:
            print(f'Warning! Found duplicate sample in well {well} of plate {self.name}')
        self.samples[well] = sampleID
    
    def generateRow(self, start, empty):
        outstr = ''
        for i in range(start, start + 12):
            outstr += '\t' + self.samples.get(i, empty)
        return outstr
    
    
class plateManager:
    
    def __init__(self, empty):
        # plate name -> plate object
        self.plates = dict()
        self.empty = empty
        
    def loadRow(self, plateID, well, sampleID):
        if plateID not in self.plates:
            self.plates[plateID] = plate(plateID)
            
        self.plates[plateID].loadSample(well, sampleID)
        
    def createPlateGrids(self, outfile):
        with open(outfile, 'w') as out:
            for k, v in self.plates.items():
                out.write(f'Plate\t{k}\n')
                out.write('\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\n')
                rowstart = 1
                for rowlabel in [chr(i) for i in range(ord('A'), ord('H') + 1)]:
                    out.write(rowlabel + v.generateRow(rowstart, self.empty) + '\n')
                    rowstart += 12
                out.write('\n')
                
        

def main(args, parser):
    worker = plateManager(args.empty)
    with open(args.file, 'r') as input:
        # clear header
        input.readline()
        for l in input:
            s = l.rstrip().split()
            if len(s) < 4:
                continue
            worker.loadRow(s[0], transform_well(s[args.row] + s[args.col]), s[args.an])

    print("Finished loading plates")
    
    worker.createPlateGrids(args.output)
    
    print("Finished printing plategrid to file")
        

if __name__ == "__main__":
    (args, parser) = arg_parse()
    main(args, parser)