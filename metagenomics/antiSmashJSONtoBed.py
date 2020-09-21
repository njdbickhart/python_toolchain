# -*- coding: utf-8 -*-
"""
This is a script designed to take a json from Antismash and turn it into a bed file
Created on Tue Jul 31 10:47:21 2018

@author: dbickhart
"""

import argparse
import json as js
import sys
import re
import os.path

TAGS = re.compile(r'[<>]')
LOC = re.compile(r'\[(\d+):(\d+)\]\((.{1})\)'

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A conversion script for antiSmash output"
            )
    parser.add_argument('-j', '--json',
                        help="Input json table with files and information for loading",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output bed file",
                        required=True, type=str,
                        )
    return parser.parse_args()

def main(args):
    fileDict = jsonParser(args.json)

    with open(args.output, 'w') as bed:
        for r in fileDict["records"]:
            if "qualifiers" not in r:
                continue

            locus = r["qualifiers"]["locus_tag"][0].split('_')[0]
            location = r["location"]

            locmatch = re.match(LOC, location)
            start = locmatch.group(1)
            end = locmatch.group(2)
            orient = locmatch.group(3)

            type = r["type"]

            bed.write(f'{locus}\t{start}\t{end}\t{orient}\t{type}\n')      

if __name__ == "__main__":
    args = arg_parse()
    main(args)
