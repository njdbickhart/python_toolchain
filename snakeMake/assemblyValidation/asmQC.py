#!/usr/bin/env python3

import argparse
import os

version = "0.0.1"

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Compare multiple genome assemblies with short-read alignment data" + version
            )
    parser.add_argument('-a', '--assembly',
                        help="An assembly fasta file to compare. Provide the full path to the file",
                        action="append", default=[]
                        )
    parser.add_argument('-n', '--name',
                        help="A name for the assembly file. Add a name in the order in which the -a option is specified.",
                        action="append", default=[]
                        )
    parser.add_argument('-b', '--busco',
                        help="The name of the BUSCO database to use for assessment",
                        type=str, required=True
                        )
    parser.add_argument('-f', '--fastq',
                        help="Paired short read fastqs separated by commas.",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--sample',
                        help="A sample name for the fastqs provided. Add in the order in which the fastq files are specified",
                        action="append", default=[]
                        )

    return parser.parse_args(), parser


if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
