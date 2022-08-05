# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 19:14:31 2022

@author: derek.bickhart-adm
"""

import sys
import argparse
import pandas as pd
import re

NAABCODE = re.compile(r'(\d+)(\D{2})(\d+)')

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A script designed to validate and padd NAAB codes"
            )
    parser.add_argument('-i', '--input', 
                        help="Input excel table with NAAB codes to validate",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name. Please add the xlsx extension!",
                        required=True, type=str,
                        )
    parser.add_argument('-c', '--column',
                        help="Column with NAAB codes to validate. Starts from Zero",
                        required=True, type=int,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    df = None
    try:
        df = pd.read_excel(args.input, header=0)
        #print(df.describe())
    except:
        print(f'Could not open excel file: {args.input}!')
        sys.exit(-1)
    
    #replacements = 0
    try:    
        df.iloc[:, args.column] = df.iloc[:, args.column].apply(nameConversion)
        #for i in range(1, len(df.iloc[:,args.column])):
        #    first, second, third = nameConversion(df.iloc[i, args.column])
        #if first == None:
        #    print(f'{df.iloc[i, args.column]} {first}')
        #first = int(first)
        #third = int(third)
        #code = f'{first:03}{second}{third:05}'
        #if df.iloc[i, args.column] != code:
        #    replacements += 1
        #df.iloc[i, args.column] = code
    except:
        print(f'Could not access column {args.column} in file {args.input}. Does the spreadsheet have that many columns?')
        sys.exit(-1)
        
    try:
        df.to_excel(args.output, index=False)
    except:
        print(f'Could not write output to file {args.output}! Please check the file name or path.')
        sys.exit(-1)
        
    print(f'Finished!')
        
    
def nameConversion(name):
    name = "".join(name.split())
    match = re.match(NAABCODE, name)
    if match:
        first, second, third = match.group(1, 2, 3)
        first = int(first)
        third = int(third)
        return f'{first:03}{second}{third:05}'
    
if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
