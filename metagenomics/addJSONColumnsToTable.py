# -*- coding: utf-8 -*-
"""
This is a script designed to take a json list of new columns to add to an existing data table
Created on Tue Jul 31 10:47:21 2018

@author: dbickhart
"""

import argparse
import json as js
import sys
from typing import Dict, List
import os.path


def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A selective tab file grep analog script"
            )
    parser.add_argument('-j', '--json', 
                        help="Input json table with files and information for loading",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Table output file name",
                        required=True, type=str,
                        )
    parser.add_argument('-t', '--table',
                        help="Input table file for processing",
                        required=True, type=str
                        )
    return parser.parse_args()

def main(args):
    # Make sure that the input JSON file is there
    if not os.path.isfile(args.json):
        print("Error reading json file! exiting...")
        sys.exit(-1)
    
    # Parse json and load files
    fileDict = jsonParser(args.json)
    
    keyOrder = fileDict.keys()
    tab = '\t'
    # Load the table file and add the additional columns in the order determined by the dictionary
    with open(args.table, 'r') as fh, open(args.output, 'w') as out:
        # Deal with the header
        head = fh.readline()
        head = head.rstrip('\n')
        
        tempCols = []
        for k in keyOrder:
            tempCols.extend(fileDict[k].colLabels)
        
        head = f'{head}\t{tab.join(tempCols)}\n'        
        out.write(head)
        
        # Now deal with the rest of the table
        for l in fh:
            l = l.rstrip('\n')
            segs = l.split("\t")
            newCols = generateNewColumns(fileDict, segs[0], keyOrder)
            out.write(f'{l}\t{tab.join(newCols)}\n')
    
    print("Done with file")
            
def generateNewColumns(fileDict : Dict, contig : str, keyOrder : List) -> List:
    tempCols = []
    for k in keyOrder:
        tempCols.extend(fileDict[k].getValues(contig))
    
    return tempCols
    
def jsonParser(json : str) -> Dict:
    with open(json, 'r') as fh:
        data = js.load(fh)
    
    # Check the number of files
    fileDict = {}
    print("Identified {} JSON objects. Processing...".format(len(data)))
    for key, value in data.items():
        print(f'Loading data from {key}...')
        fileDict[key] = jsonData(value)
        fileDict[key].processFile()
        
    print("All detected json entries loaded...")
    return fileDict

class jsonData:
    
    def __init__(self, jsonDict : Dict) -> None:
        self.file = jsonDict["File"]
        self.default = jsonDict["Default"]
        self.numColumns = int(jsonDict["NumCol"])
        self.colLabels = jsonDict["ColLabels"]
        self.data = {}
        
    def processFile(self) -> None:
        with open(self.file, "r") as fh:
            for l in fh:
                l = l.rstrip('\n')
                segs = l.split()
                self.data[segs[0]] = segs[1:]
                
    def getValues(self, contig : str) -> List:
        if contig not in self.data:
            # Return empty values for the number of columns needed
            return [self.default for x in range(self.numColumns)]
        else:
            return self.data[contig]
                
    
if __name__ == "__main__":
    args = arg_parse()
    main(args)
