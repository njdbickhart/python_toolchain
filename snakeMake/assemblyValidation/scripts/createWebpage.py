#!/usr/bin/env python3

import argparse
import os


def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Compare multiple genome assemblies with short-read alignment data" + version
            )
    parser.add_argument('-f', '--final',
                        help="The location of the final folder of the pipeline",
                        type=str, required=True
                        )
    parser.add_argument('-o', '--output',
                        help="output file basename ",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--combos',
                        help="Expected assembly combinations",
                        action="append", default=[]
                        )

    return parser.parse_args(), parser


if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
