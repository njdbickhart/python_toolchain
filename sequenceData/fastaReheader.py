# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:24:09 2021

@author: derek.bickhart-adm
"""

import argparse
import re
import sys
import os

# regular expression to extract header sequence
FAHEADER = re.compile(r'>(\S+)')
# regular expression to add newlines to fasta sequence to comply with fasta file format guidelines
FALINES = re.compile(r'(.{60})')
# regular expression to remove all potential confounding delimiters from the fasta name
DELIMITERS = re.compile(r'[:-|,_]')


def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to process a fasta file. If run with only the fasta file name, will create a table file that the user can edit to rename the fasta\n" +
            "To generate initial table: python fastaReheader -f file.fasta\n" +
            "To generate reheadered fasta: python fastaReheader -f file.fasta -t file.fasta.table -o reheadered.fasta\n"
            )
    parser.add_argument('-f', '--fasta', 
                        help="Input reference fasta file. If run as the only option will create a table file [filename + .table] suitable for option -t",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output fasta file name. Only used with the -t option",
                        required=False, type=str, default="NO"
                        )
    parser.add_argument('-t', '--table',
                        help="Input table file. Used to convert fasta file entry names",
                        required=False, type=str, default="NO"
                        )
    parser.add_argument('-d', '--debug',
                        help="Provide verbose output to debug code",
                        action="store_true"
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    tableMode = True if args.table == "NO" and args.output == "NO" else False 
    parseMode = True if args.table != "NO" and args.output != "NO" else False
    
    if not tableMode and not parseMode:
        print("Error with input arguments! Could not determine which mode to use! Please select one of the options displayed in the help menu")
        parser.print_help()
        sys.exit(-1)
        
    if args.debug:
        print("Debug selected. Starting script")
        
    if tableMode:
        if os.path.exists(args.fasta + ".table"):
            print("Started to generate a table file but found an older table file that would have been overwritten. Do you want to overwrite the file?")
            x = input()
            if x[0].lower() != 'y':
                print("Stopping the program!")
                sys.exit(-1)
            else:
                print("OK, overwriting file...")
                
        tableGeneration(args.fasta, args.fasta + ".table", args.debug)
        print("Done! Please edit the table file and use it in the next iteration of the program to edit your fasta")
    elif parseMode:
        if args.debug:
            print("Starting to convert fasta file...")
            
        replaceHeader(args.fasta, args.table, args.output, args.debug)
    else:
        print("You shouldn't be here!")
        
    if args.debug:
        print("That's all folks!")

def tableGeneration(filename, outputname, debug):
    with open(filename, 'r') as input, open(outputname, 'w') as output:
        for name, seq in fasta_reader_fh(input):
            if debug:
                print(f'{name}\t{len(seq)}')
            (nname, nsubs) = re.subn(DELIMITERS, '', name)
            if debug:
                print(f'Changed fasta header name {name} to {nname} and replaced {nsubs} delimiters')
            output.write(f'{name}\t{nname}\n')
            
def replaceHeader(filename, tablefile, outputname, debug):
    
    converter = dict()
    with open(tablefile, 'r') as input:
        for l in input:
            s = l.rstrip().split()
            converter[s[0]] = s[1]
            
    if debug:
        print("Loaded conversion table file")
        
    with open(filename, 'r') as input, open(outputname, 'w') as output:
        for name, seq in fasta_reader_fh(input):
            if name not in converter:
                print(f'ERROR! Could not find {name} in the tablefile: {tablefile}!')
                print("Exiting...")
                sys.exit(-1)
                
            nname = converter[name]
            (seq, nsubs) = re.subn(FALINES, r'\1\n', seq)
            if debug:
                print(f'For {name} -> {nname} added {nsubs} newlines to sequence')
                
            output.write(f'>{nname}\n{seq}\n')
    
    print(f'File successfully converted. Output is: {outputname}')
    
def fasta_reader_fh(infile):
    name = infile.readline().rstrip()
    name = re.match(FAHEADER, name).group(1)
    while True:
        seq = ""
        for s in infile:
            if s[0] == '>':
                yield name, seq 
                name = s.rstrip()
                re.match(FAHEADER, name).group(1)
                break
            else:
                seq += s.rstrip()        
        else:
            yield name, seq
            return
  
if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
