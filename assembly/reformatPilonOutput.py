# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:01:45 2020

@author: derek.bickhart-adm
"""

import argparse
import re
from os import path
import io
import subprocess as sp
import glob
import sys
import concurrent.futures

isChr = re.compile(r'[Cc]hr(.{1,2})_?')
isNum = re.compile(r'(\d+)')
pTag = re.compile(r'_pilon')

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Resolve pilon fastas into karyotypic order"
            )
    parser.add_argument('-f', '--folder',
                        help="The base folder where the pilon fastas are contained",
                        required=False, type=str, default=None
                        )
    parser.add_argument('-b', '--file',
                        help="OR the fasta file to edit",
                        required=False, type=str, default=None
                        )
    parser.add_argument('-o', '--output',
                        help="Output file full name",
                        required=True, type=str,
                        )
    parser.add_argument('-t', '--threads',
                        help="Max threads to use",
                        default=1, type=int,
                        )
    return parser.parse_args(), parser

def main(args, parser):
    if args.folder == None and args.file == None:
        parser.print_help()
        sys.exit()

    # Get a list of fasta files in the directory
    files = []
    if args.file == None:
        for exts in ('*.fa', '*.fasta'):
            files.extend(glob.glob('{}/{}'.format(args.folder, exts)))
    else:
        files.append(args.file)

    # Now samtools faidx them all and tabulate all chr names
    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        for f in files:
            results[f] = executor.submit(samtoolsFaiTask, f)

    # Now organize results for sequential samtools faidx processing
    completed = {}
    chrlist = []
    for k, v in results.items():
        for c in v.result():
            chrlist.append(c)
            completed[c] = k

    # Sort the list by karyotype and then process
    chrlist.sort(key=lambda x: karyotypeSort(x))
    with open(args.output, 'w') as out:
        for c in chrlist:
            f = completed[c]
            print(f'Working on chr: {c} from file {f}')
            with sp.Popen(f'samtools faidx {f} {c}', shell=True, stdout=sp.PIPE, bufsize=1, universal_newlines=True) as sf:
                h = sf.stdout.readline().rstrip()
                h = re.sub(pTag, '', h)
                out.write(f'{h}\n')
                for l in sf.stdout:
                    out.write(l)


def samtoolsFaiTask(file):
    sp.run(f'samtools faidx {file}', shell=True)
    chrs = []
    if path.isfile(f'{file}.fai'):
        with open(f'{file}.fai', 'r') as input:
            for l in input:
                s = l.rstrip().split()
                chrs.append(s[0])
    return chrs

def karyotypeSort(key):
    c = re.match(isChr, key)
    d = re.match(isNum, key)
    if c:
        v = c.group(1)
        if not re.match(isNum, v):
            if v == 'X':
                v = 100
            elif v == 'Y':
                v = 101
            else:
                v = 102
        return int(v)
    elif d:
        v = d.group(1)
        return int(v)
    return 103



if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
